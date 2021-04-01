program uiuj_largesmall
use dnsdata
use ifport  ! this library is intel compiler specific; only used to create a directory (makedirqq)
            ! please comment this and substitute the call to makedirqq if using gfortran or other
implicit none

! stuff for parsing arguments

character(len=32) :: cmd_in_buf ! needed to parse arguments
character(len=32) :: arg

integer(C_INT) :: nfmin,nfmax,dnf,nftot
integer :: nmin,nmax,dn
integer :: nmin_cm = 0, nmax_cm = 0, dn_cm = 0

logical :: custom_mean = .FALSE.

! counters

integer :: ii ! counter initially used for parsing arguments
integer :: ix, iy, iz, iv ! for spatial directions and components
integer :: i, j, k, c, irs, ilasm

! global stuff

complex(C_DOUBLE_COMPLEX), allocatable :: pressure(:,:,:)
complex(C_DOUBLE_COMPLEX), allocatable :: dertemp_in(:), dertemp_out(:)
complex(C_DOUBLE_COMPLEX), allocatable :: Vgrad(:,:,:,:,:)

integer, parameter :: file_vel = 883, file_press = 884
character(len=40) :: istring, foldername, currfname

real(C_DOUBLE), allocatable :: mean(:,:), uiujprofiles(:,:,:,:)
real(C_DOUBLE) :: m_grad(3,3) ! m_grad(i,j) = dUi/dxj

! shortcut parameters
integer, parameter :: var = 1, prod = 2, psdiss = 3, ttrsp = 4, tcross = 5, vdiff = 6, pstrain = 7, ptrsp = 8, PHIttrsp = 9, PHIvdiff = 10, PHIptrsp = 11
integer, parameter :: large = 1, small = 2
logical :: ignore

! large-small stuff

integer :: z_threshold

! MPI stuff

TYPE(MPI_Datatype) :: press_read_type, press_field_type, vel_read_type, vel_field_type
TYPE(MPI_Datatype) :: mean_write_type, mean_inmem_type, uiuj_write_type, uiuj_inmem_type
TYPE(MPI_File) :: fh
integer :: ierror
integer(MPI_OFFSET_KIND) :: offset


! Program begins here
!----------------------------------------------------------------------------------------------------------------------------



    ! read arguments from command line
    call parse_args()
    if (custom_mean) then
        ! apply settings
        nfmin = nmin_cm; nfmax = nmax_cm; dnf = dn_cm
    end if
    nftot = ((nfmax - nfmin) / dnf) + 1
    
    ! Init MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call read_dnsin()
    call largesmall_setup()
    ! set npy to total number of processes
    ! this effectively DEACTIVATES PARALLELISATION IN X/Z
    ! and makes life easier (without significant losses in perf.)
    npy = nproc
    ! notice that this overrides value from dns.in
    call init_MPI(nx+1,nz,ny,nxd+1,nzd)
    call init_memory(.FALSE.) ! false flag avoids allocation of RHS related stuff!
    call init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB)
    call setup_derivatives()

    ! allocate stuff
    call init_uiuj_mpitypes()
    allocate(pressure(ny0-2:nyN+2,-nz:nz,nx0:nxN))
    allocate(Vgrad(ny0-2:nyN+2,-nz:nz,nx0:nxN,1:3,1:3)) ! iy, iz, ix, iv (veloctiy), ider (direction of derivative)
    allocate(dertemp_in(ny0-2:nyN+2))
    allocate(dertemp_out(ny0-2:nyN+2))
    allocate(mean(1:7, ny0-2:nyN+2))
    allocate(uiujprofiles(1:6, 1:11, ny0-2:nyN+2, 1:2))

    !---------------------------------------------------!
    !----------------- COMPUTE AVERAGE -----------------!
    !---------------------------------------------------!

    if (has_terminal) print *, "Computing average..."

    mean = 0
    do ii = nfmin, nfmax, dnf ! loop over files
        write(istring,*) ii
        open(file="Dati.cart."//TRIM(ADJUSTL(istring))//".out", unit=file_vel, status="old", access="stream", action="read")
        open(file="pField"//TRIM(ADJUSTL(istring))//".fld", unit=file_press, status="old", access="stream", action="read")
            do iy = ny0-2, nyN+2
                ! cumulate mean data
                mean(1,iy) = mean(1,iy) + real(fieldmap(file_vel, iy, 0, 0, 1)) ! U
                mean(2,iy) = mean(2,iy) + real(fieldmap(file_vel, iy, 0, 0, 3)) ! W
                mean(7,iy) = mean(7,iy) + real(pressmap(file_press, iy, 0, 0) ) ! P
            end do
        close(file_vel)
        close(file_press)
    end do
    mean = mean / nftot ! divide by number of files

    ! derivate stuff
    CALL REALderiv(mean(1,:), mean(3,:))  ! U --> Uy
    CALL REALderiv2(mean(1,:), mean(5,:)) ! U --> Uyy
    CALL REALderiv(mean(2,:), mean(4,:))  ! W --> Wy
    CALL REALderiv2(mean(2,:), mean(6,:)) ! W --> Wyy

    ! revert to desired indices for calculation of tke
    nfmin = nmin; nfmax = nmax; dnf = dn
    nftot = ((nfmax - nfmin) / dnf) + 1

    !---------------------------------------------------!
    !--------------- COMPUTE TKE BUDGET ----------------!
    !---------------------------------------------------!

    uiujprofiles = 0;

    do ii = nfmin, nfmax, dnf ! loop over files

        ! read velocity, pressure
        write(istring,*) ii
        currfname = trim("Dati.cart."//TRIM(ADJUSTL(istring))//".out")
        call read_vel(currfname, V)
        currfname = trim("pField"//TRIM(ADJUSTL(istring)))//".fld"
        call read_press(currfname, pressure)

        ! remove average
        V(:,0,0,1) = V(:,0,0,1) - cmplx(mean(1,:))
        V(:,0,0,3) = V(:,0,0,3) - cmplx(mean(2,:))
        pressure(:,0,0) = pressure(:,0,0) - cmplx(mean(7,:))

        ! compute derivatives
        call get_gradient_33(V,Vgrad)

#define u(cmp) V(iy,iz,ix,cmp)
#define gu(cmp,dd) Vgrad(iy,iz,ix,cmp,dd)
#define uiuj(trm) uiujprofiles(irs,trm,iy,ilasm)
#define p pressure(iy,iz,ix)

        ! parseval theorem method for var, prod, pstrain, psdiss, PHIptrsp
        do iy = ny0-2, nyN+2
            do ix = nx0, nxN
                c = 2; if (ix == 0) c = 1 ! multiplier for doubling points in x direction
                do iz = -nz, nz
                    ilasm = small; if (is_large(ix,iz)) ilasm = large ! determine if mode is large or small
                    do irs = 1, 6
                        
                        ! prepare stuff
                        call get_indexes(irs, i, j)
                        call get_mean_grad(iy)
                        
                        ! calculate statistics
                        uiuj(var) = uiuj(var) + c*cprod(u(i), u(j))
                        uiuj(pstrain) = uiuj(pstrain) + c*cprod(p, gu(i,j)+gu(j,i) )
                        do k = 1,3 ! prod, psdiss need sum over k
                            uiuj(prod) = uiuj(prod) - c * ( cprod(u(i),u(k))*m_grad(j,k) + cprod(u(j),u(k))*m_grad(i,k) )
                            uiuj(psdiss) = uiuj(psdiss) - 2 * ni * c * cprod(gu(i,k), gu(j,k))
                        end do
                        select case(irs) ! PHIptrsp: only for specific terms
                            case(2)
                                uiuj(PHIptrsp) = uiuj(PHIptrsp) - c * ( cprod(p,u(2)) + cprod(u(2),p) )
                            case(4)
                                uiuj(PHIptrsp) = uiuj(PHIptrsp) - c * cprod(p,u(1))
                            case(5)
                                uiuj(PHIptrsp) = uiuj(PHIptrsp) - c * cprod(p,u(3))
                        end select

                    end do
                end do
            end do
        end do

        ! antitransform

        ! calculate ttrsp and tcross

    end do

    ! divide by nftot to obtain average
    uiujprofiles = uiujprofiles / nftot

    ! derivate fluxes

!    ! read file
!    CALL read_vel(fname,V)
!
!    ! for each y plane
!    do iy = miny,maxy ! skips halo cells: only non-duplicated data
!        
!        ! skip ghost cells
!        if (iy < 0) cycle
!        if (iy > ny) cycle
!        
!        ! do fourier transform
!        VVdz(1:nz+1,1:nxB,1:3,1)=V(iy,0:nz,nx0:nxN,1:3);         VVdz(nz+2:nzd-nz,1:nxB,1:3,1)=0;
!        VVdz(nzd+1-nz:nzd,1:nxB,1:3,1)=V(iy,-nz:-1,nx0:nxN,1:3); 
!        DO iV=1,3
!            CALL IFT(VVdz(1:nzd,1:nxB,iV,1))  
!            CALL MPI_Alltoall(VVdz(:,:,iV,1), 1, Mdz, VVdx(:,:,iV,1), 1, Mdx, MPI_COMM_X)
!            VVdx(nx+2:nxd+1,1:nzB,iV,1)=0;    CALL RFT(VVdx(1:nxd+1,1:nzB,iV,1),rVVdx(1:2*nxd+2,1:nzB,iV,1))
!        END DO
!        ! rVVdx is now containing the antitransform
!        ! it has size 2*(nxd+1),nzB,6,6
!        ! but you only need (1:(2*nx+1),1:(2*nz+1),1:3,1)
!
!        ! TODO your stuff here
!
!    end do

    !---------------------------------------------------!
    !-----------------      WRITE      -----------------!
    !---------------------------------------------------!

    ! create folder if not existing
    if (has_terminal) then
        if (custom_mean) then
            foldername = "cm_largesmall"
        else
            foldername = "largesmall"
        end if
        ignore = makedirqq(foldername)
    end if
    
    ! write to disk
    if (has_terminal) print *, "Saving to disk..."
    CALL MPI_File_open(MPI_COMM_WORLD, "uiuj_largesmall.bin", IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
        
        ! write header
        if (has_terminal) CALL MPI_file_write(fh, [nfmin, nfmax, dnf, nftot], 4, MPI_INTEGER, MPI_STATUS_IGNORE)

        ! write mean data
        offset = 4 * sizeof(nfmin)
        CALL MPI_File_set_view(fh, offset, MPI_DOUBLE_PRECISION, mean_write_type, 'native', MPI_INFO_NULL)
        CALL MPI_File_write_all(fh, mean, 1, mean_inmem_type, MPI_STATUS_IGNORE)

        ! write uiuj data
        offset = offset + sizeof(uiujprofiles)
        CALL MPI_File_set_view(fh, offset, MPI_DOUBLE_PRECISION, uiuj_write_type, 'native', MPI_INFO_NULL)
        CALL MPI_File_write_all(fh, uiujprofiles, 1, uiuj_inmem_type, MPI_STATUS_IGNORE)

    call MPI_File_close(fh)

    !---------------------------------------------------!
    !-----------------     FINALISE    -----------------!
    !---------------------------------------------------!

    ! realease memory
    CALL free_fft()
    CALL free_memory(.FALSE.) 
    CALL MPI_Finalize()



contains !-------------------------------------------------------------------------------------------------------------------



    ! initialise largesmall stuff
    subroutine largesmall_setup()
        open(15, file='largesmall_settings.in')
            read(15, *) z_threshold
        close(15)
        if (has_terminal) then
            write(*,"(A,I5)") "   z_threshold =", z_threshold
        end if
    end subroutine largesmall_setup



    ! returns true if a given mode is part of large scale field
    function is_large(xx,zz) result(islarge)
    integer, intent(in) :: xx, zz
    logical :: islarge 
        islarge = .FALSE.
        if (abs(beta0*zz) <= z_threshold) then
            islarge = .TRUE.
        end if
    end function is_large



    ! define and commit MPI filetypes
    subroutine init_uiuj_mpitypes()

        ! define type for reading pressure from disk
        CALL MPI_Type_create_subarray(3, [ny+3, 2*nz+1, nx+1], [nyN-ny0+5, 2*nz+1, nxB], [ny0-1,0,nx0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, press_read_type, ierror)
        CALL MPI_Type_commit(press_read_type, ierror)
        ! type describing pressure in memory
        CALL MPI_Type_create_subarray(3, [nyN-ny0+5, 2*nz+1, nxB], [nyN-ny0+5, 2*nz+1, nxB], [0,0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, press_field_type, ierror)
        CALL MPI_Type_commit(press_field_type, ierror)

        ! define type for reading velocity from disk
        CALL MPI_Type_create_subarray(4, [ny+3, 2*nz+1, nx+1, 3], [nyN-ny0+5, 2*nz+1, nxB, 3], [ny0-1,0,nx0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_read_type, ierror)
        CALL MPI_Type_commit(vel_read_type, ierror)
        ! type describing velocity in memory
        CALL MPI_Type_create_subarray(4, [nyN-ny0+5, 2*nz+1, nxB, 3], [nyN-ny0+5, 2*nz+1, nxB, 3], [0,0,0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_field_type, ierror)
        CALL MPI_Type_commit(vel_field_type, ierror)

        ! define type for writing mean on disk
        CALL MPI_Type_create_subarray(2, [7, ny+1], [7, min(ny,maxy)-max(0,miny)+1], [0,max(0,miny)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mean_write_type, ierror)
        CALL MPI_Type_commit(mean_write_type, ierror)
        CALL MPI_Type_create_subarray(2, [7, nyN-ny0+5], [7, min(ny,maxy)-max(0,miny)+1], [0,max(0,miny)-(ny0-2)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mean_inmem_type, ierror)
        CALL MPI_Type_commit(mean_inmem_type, ierror)

        ! define type for writing uiujprofiles on disk
        CALL MPI_Type_create_subarray(3, [6, 11, ny+1, 2], [6, 11, min(ny,maxy)-max(0,miny)+1, 2], [0,0,max(0,miny),0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, uiuj_write_type, ierror)
        CALL MPI_Type_commit(uiuj_write_type, ierror)
        CALL MPI_Type_create_subarray(3, [6, 11, nyN-ny0+5, 2], [6, 11, min(ny,maxy)-max(0,miny)+1, 2], [0,0,max(0,miny)-(ny0-2),0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, uiuj_inmem_type, ierror)
        CALL MPI_Type_commit(uiuj_inmem_type, ierror)

    end subroutine



    ! access pressure file on unit as a memmap
    function pressmap(unit, iy_zero, iz_zero, ix_zero) result(element)
        integer, intent(in) :: unit, ix_zero, iy_zero, iz_zero
        integer(C_SIZE_T) :: ix, iy, iz
        complex(C_DOUBLE_COMPLEX) :: element
        integer(C_SIZE_T) :: position, el_idx ! position
        
        ! calculate indices starting from 1
        ix = ix_zero + 1
        iz = iz_zero + nz + 1
        iy = iy_zero + 2

        ! calculate position
        el_idx = (ix-1_C_SIZE_T)*(2_C_SIZE_T*nz+1_C_SIZE_T)*(ny+3_C_SIZE_T) + (iz-1_C_SIZE_T)*(ny+3_C_SIZE_T) + iy
        position = 1_C_SIZE_T + (el_idx - 1) * 2*C_DOUBLE_COMPLEX

        ! read element
        read(unit,pos=position) element

    end function pressmap



    subroutine read_press(filename,R)
        complex(C_DOUBLE_COMPLEX), intent(inout) :: R(ny0-2:nyN+2,-nz:nz,nx0:nxN)
        character(len=40), intent(IN) :: filename
        TYPE(MPI_File) :: fh

        call MPI_file_open(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh)
        call MPI_file_set_view(fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE_COMPLEX, press_read_type, 'native', MPI_INFO_NULL)
        call MPI_file_read_all(fh, R, 1, press_field_type, MPI_STATUS_IGNORE)
        call MPI_file_close(fh)

    end subroutine read_press



    subroutine read_vel(filename,R)
        complex(C_DOUBLE_COMPLEX), intent(inout) :: R(ny0-2:nyN+2,-nz:nz,nx0:nxN,3)
        character(len=40), intent(IN) :: filename
        INTEGER(MPI_OFFSET_KIND) :: disp = 3*C_INT + 7*C_DOUBLE
        TYPE(MPI_File) :: fh

        if (has_terminal) print *, "Reading from file "//filename

        call MPI_file_open(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh)
        call MPI_file_set_view(fh, disp, MPI_DOUBLE_COMPLEX, vel_read_type, 'native', MPI_INFO_NULL)
        call MPI_file_read_all(fh, R, 1, vel_field_type, MPI_STATUS_IGNORE)
        call MPI_file_close(fh)

    end subroutine read_vel



    subroutine get_gradient_33(R, grad)
        complex(C_DOUBLE_COMPLEX), intent(in) :: R(ny0-2:nyN+2,-nz:nz,nx0:nxN,3)
        complex(C_DOUBLE_COMPLEX), intent(out) :: grad(ny0-2:nyN+2,-nz:nz,nx0:nxN,3,3)

        ! y derivatives
        do iv = 1,3
            do ix = nx0, nxN
                do iz = -nz, nz
                    ! stuff
                    call COMPLEXderiv(R(:,iz,ix,iv), grad(:,iz,ix,iv,2))
                    grad(:,iz,ix,iv,1) = R(:,iz,ix,iv) * alfa0 * ix
                    grad(:,iz,ix,iv,3) = R(:,iz,ix,iv) * beta0 * iz
                end do
            end do
        end do
    end subroutine get_gradient_33



    subroutine get_indexes(irs, i, j)
    integer, intent(in) :: irs
    integer, intent (out) :: i, j
        select case(irs)
            case(1)
                i=1; j=1
            case(2)
                i=2; j=2
            case(3)
                i=3; j=3
            case(4)
                i=1; j=2
            case(5)
                i=2; j=3
            case(6)
                i=1; j=3
        end select
    end subroutine get_indexes



    subroutine get_mean_grad(iy)
    integer, intent(in) :: iy
        m_grad = 0
        m_grad(1,2) = mean(3,iy)
        m_grad(3,2) = mean(4,iy)
    end subroutine



    function cprod(a, b) result(r)
    complex(C_DOUBLE_COMPLEX), intent(in) :: a, b
    real(C_DOUBLE) :: r
        r = real(conjg(a)*b)
    end function cprod



!----------------------------------------------------------------------------------------------------------------------------
! less useful stuff here
!----------------------------------------------------------------------------------------------------------------------------



    subroutine print_help()
        print *, "Calculates TKE budget; sharp Fourier filtering is used to decompose the fluctuation field into large and small components."
        print *, "Statistics are calculated on files ranging from index nfmin to nfmax with step dn. Usage:"
        print *, ""
        print *, "   mpirun [mpi args] uiuj_largesmall [-h] nfmin nfmax dn [--custom_mean nmin_m nmax_m dn_m]"
        print *, ""
        print *, "If the flag --custom_mean is passed, the mean field is calculated on fields (nmin_m nmax_m dn_m); the remaining statistics are still calculated on (nfmin nfmax dn)."
        print *, ""
        print *, "This program is meant to be used on plane channels."
        print *, ""
        print *, "Results are output to uiuj.bin. Use uiuj2ascii to get the results in a human readable format."
        print *, ""
        print *, "Mean TKE budget terms are calculated as:"
        print *, "INST    --> dK/dt"
        print *, "CONV    --> Ui*dK/dxi"
        print *, "PROD    --> -<uiuj>dUj/dxi"
        print *, "DISS*   --> nu<(duj/dxi + dui/dxj)*duj/dxi>"
        print *, "TDIFF   --> -0.5*d/dxi<ui*uj*uj>"
        print *, "PDIFF   --> -d/dxi<ui*p>"
        print *, "VDIFF1  --> nu*d2K/dxi2"
        print *, "VDIFF2* --> nu*d2/dxjdxi<ui*uj>"
        print *, "*-terms can be summed into the PDISS=nu*<duj/dxi*duj/dxi>"
        print *, ""
        print *, "which in a statistically stationary and fully-developed turbulent"
        print *, "channel flow with spanwise wall oscillations reduces to"
        print *, "PROD  --> -<uv>dU/dy-<vw>dW/dy         [this is computed after the fields loop]"
        print *, "PDISS --> nu*<dui/dxj*dui/dxj>"
        print *, "TDIFF --> -0.5*d/dy(<vuu>+<vvv>+<vww>)"
        print *, "PDIFF --> -d/dy<vp>"
        print *, "VDIFF --> nu*d2K/dy2"
        print *, ""
        print *, "The MKE buget equation, in a statistically stationary"
        print *, "and fully-developed turbulent channel flow, reduces to"
        print *, "pump   --> -dP/dx*U = tau_w*U"
        print *, "TPROD  --> <uv>dU/dy+<vw>dW/dy         [TKE production, here a sink]"
        print *, "ttrsp  --> -d(<uv>U)/dy-d(<vw>W)/dy"
        print *, "vdiff  --> ni*d(U*Uy)/dy+d(W*Wy)/dy"
        print *, "dissU  --> ni*dUdy^2"
        print *, "dissW  --> ni*dWdy^2"
        stop
    end subroutine print_help



    ! read arguments
    subroutine parse_args()
        ii = 1
        if (command_argument_count() < 3) then ! handle exception: no input
            print *, 'ERROR: please provide one input file as command line argument.'
            stop
        end if
        do while (ii <= command_argument_count()) ! parse optional arguments
            call get_command_argument(ii, arg)
            select case (arg)
                case ('-h', '--help') ! call help
                    call print_help()
                    stop
                case ('-c', '--custom_mean') ! specify undersampling
                    custom_mean = .TRUE.
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) nmin_cm
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) nmax_cm
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) dn_cm
                case default
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) nmin
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) nmax
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) dn
                    ! apply settings
                    nfmin = nmin; nfmax = nmax; dnf = dn
            end select
            ii = ii + 1
        end do
    end subroutine parse_args



end program



!----------------------------------------------------------------------------------------------------------------------------



!--------------------------------------------------!
!-----------------     LEGEND     -----------------!
!--------------------------------------------------!


! mean(iv, iy)      where iv        U   W   Uy  Wy  Uyy Wyy P
!                                   1   2   3   4   5   6   7


! irs ----> index used for Reynolds stress
!   uu  vv  ww  uv  vw  uw
!   1   2   3   4   5   6

! uiujprofiles(irs, iterm, iy, largesmall)      where irs defined before, while iterm:
!   var prod    psdiss      ttrsp   tcross  vdiff   pstrain     ptrsp   PHIttrsp    PHIvdiff    PHIptrsp
!   1   2       3           4       5       6       7           8       9           10          11


