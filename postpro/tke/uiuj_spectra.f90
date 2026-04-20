
!--------------------------------------------------!
!-----------------     LEGEND     -----------------!
!--------------------------------------------------!

! mean(iv, iy)      where iv        U   W   Uy  Wy  Uyy Wyy P
!                                   1   2   3   4   5   6   7

! irs ----> index used for Reynolds stress
!   uu  vv  ww  uv  vw  uw
!   1   2   3   4   5   6

! uiujspectra(irs, iterm, iz, iy)      where irs defined before, while iterm:
!   var prod    psdiss      ttrsp   vdiff   pstrain     ptrsp   PHIttrsp    PHIvdiff    PHIptrsp
!   1   2       3           4       5       6           7       8           9           10      

! *VVd*(zf(iz), xf(ix), term, 1)
! : WARNING :  iz, ix indices here is scrambled as required by FFT -> use zf(iz), xf(ix)
! Index term has 3 possible values = 1:3
! - term=1 -> dui / dxk
! - term=2 -> duj / dxk
! - term=3 -> uk
! - term=4 -> uk * ( dui / dxk )
! - term=5 -> uk * ( duj / dxk )
! See parameters i_ft, j_ft, k_ft, ik_ft, jk_ft

! : WARNING : triadic interactions are still not implemented
! Array for triadic interactions: convs
! convs(itri,-nz:nz,0:nz,iy)
! itri = 1,4
!   uu  vv  ww  uv
!   1   2   3   4



#include "header.h"

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

logical :: custom_mean = .FALSE., timestwo = .FALSE.

! counters

integer :: ii ! counter initially used for parsing arguments
integer :: ix, iy, iz, iv ! for spatial directions and components
integer :: i, j, k, c, irs, cntr

! global stuff

complex(C_DOUBLE_COMPLEX), allocatable :: pressure(:,:,:)
complex(C_DOUBLE_COMPLEX), allocatable :: dertemp_in(:), dertemp_out(:)
complex(C_DOUBLE_COMPLEX), allocatable :: Vgrad(:,:,:,:,:)

integer, parameter :: file_vel = 883, file_press = 884
character(len=40) :: istring, currfname

real(C_DOUBLE), allocatable :: mean(:,:), uiujspectra(:,:,:,:), uiujprofiles(:,:,:), convs(:,:,:,:)
real(C_DOUBLE) :: m_grad(3,3) ! m_grad(i,j) = dUi/dxj

character(40), parameter :: nfofile="uiuj_spectra.nfo" 

! shortcut parameters

integer, parameter :: var = 1, prod = 2, psdiss = 3, ttrsp = 4, vdiff = 5, pstrain = 6, ptrsp = 7, PHIttrsp = 8, PHIvdiff = 9, PHIptrsp = 10
integer, parameter :: i_ft = 1, j_ft = 2, k_ft = 3, ik_ft = 4, jk_ft = 5

! MPI stuff

TYPE(MPI_Datatype) :: press_read_type, press_field_type
TYPE(MPI_Datatype) :: mean_write_type, mean_inmem_type, spectra_write_type, spectra_inmem_type, profiles_write_type, profiles_inmem_type, convs_write_type, convs_inmem_type
TYPE(MPI_File) :: fh
integer :: ierror
integer(MPI_OFFSET_KIND) :: offset

! Program begins here
!----------------------------------------------------------------------------------------------------------------------------

    ! Init MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

    ! read arguments from command line
    call parse_args()
    call check_in_files()
    if (custom_mean) then
        ! apply settings
        nfmin = nmin_cm; nfmax = nmax_cm; dnf = dn_cm
    end if
    nftot = ((nfmax - nfmin) / dnf) + 1
    
    ! setup
    call read_dnsin(.TRUE.) ! passing .TRUE. cause dns.in is in parent folder wrt cwd
    ! set npy to total number of processes
    ! this effectively DEACTIVATES PARALLELISATION IN X/Z
    ! and makes life easier (without significant losses in perf.)
    npy = nproc
    ! notice that this overrides value from dns.in
    call init_MPI(nx+1,nz,ny,nxd+1,nzd)
    call init_memory(.FALSE.) ! false flag avoids allocation of RHS related stuff!
    call init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB,.FALSE.,[5,1]) ! size of fft array is specified with optional flag

    call setup_derivatives()

    ! allocate stuff
    call init_uiuj_mpitypes()
    allocate(pressure(ny0-2:nyN+2,-nz:nz,nx0:nxN))
    allocate(Vgrad(ny0-2:nyN+2,-nz:nz,nx0:nxN,1:3,1:3)) ! iy, iz, ix, iv (veloctiy), ider (direction of derivative)
    allocate(dertemp_in(ny0-2:nyN+2))
    allocate(dertemp_out(ny0-2:nyN+2))
    allocate(mean(1:7, ny0-2:nyN+2))
    allocate(uiujspectra(1:6, 1:10, -nz:nz, ny0-2:nyN+2))
    allocate(uiujprofiles(1:6, 1:10, ny0-2:nyN+2))
    allocate(convs(1:4,-nz:nz,0:nz,ny0-2:nyN+2)); convs=0

    !---------------------------------------------------!
    !----------------- COMPUTE AVERAGE -----------------!
    !---------------------------------------------------!

    if (has_terminal) print *, "Computing average..."

    mean = 0
    do ii = nfmin, nfmax, dnf ! loop over files
        write(istring,*) ii
        open(file="../Dati.cart."//TRIM(ADJUSTL(istring))//".out", unit=file_vel, status="old", access="stream", action="read")
        open(file="../pField"//TRIM(ADJUSTL(istring))//".fld", unit=file_press, status="old", access="stream", action="read")
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
    call REALderiv(mean(1,:), mean(3,:))  ! U --> Uy
    call REALderiv2(mean(1,:), mean(5,:)) ! U --> Uyy
    call REALderiv(mean(2,:), mean(4,:))  ! W --> Wy
    call REALderiv2(mean(2,:), mean(6,:)) ! W --> Wyy

    ! revert to desired indices for calculation of tke
    nfmin = nmin; nfmax = nmax; dnf = dn
    nftot = ((nfmax - nfmin) / dnf) + 1

    !---------------------------------------------------!
    !--------------- COMPUTE TKE BUDGET ----------------!
    !---------------------------------------------------!

    uiujspectra = 0;

    do ii = nfmin, nfmax, dnf ! loop over files

        ! read velocity, pressure
        write(istring,*) ii
        currfname = trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
        call read_restart_file(currfname, V)
        currfname = trim("../pField"//TRIM(ADJUSTL(istring)))//".fld"
        call read_press(currfname, pressure)

        ! remove average
        V(:,0,0,1) = V(:,0,0,1) - cmplx(mean(1,:))
        V(:,0,0,3) = V(:,0,0,3) - cmplx(mean(2,:))
        pressure(:,0,0) = pressure(:,0,0) - cmplx(mean(7,:))

        ! compute derivatives
        call get_gradient_33(V,Vgrad)

#define u(cmp) V(iy,iz,ix,cmp)
#define gu(cmp,dd) Vgrad(iy,iz,ix,cmp,dd)
#define uiuj(trm) uiujspectra(irs,trm,iz,iy)
#define p pressure(iy,iz,ix)

        ! parseval theorem method for var, prod, pstrain, psdiss, PHIptrsp
        do iy = ny0-2, nyN+2
            do ix = nx0, nxN
                c = 2; if (ix == 0) c = 1 ! multiplier for doubling points in x direction
                do iz = -nz, nz
                    do irs = 1, 6
                        
                        ! prepare stuff
                        call get_indexes(irs, i, j)
                        call get_mean_grad(iy)
                        
                        ! calculate statistics
                        uiuj(var) = uiuj(var) + c*cprod( u(i), u(j) )
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

#define xf(a) a-nx0+1

            ! calculate ttrsp
            do irs = 1, 6

                call get_indexes(irs, i, j)

                ! loop over k
                do k = 1,3
                    
                    ! first compute transofrm of ui, uj
                    VVdz(:,:,i_ft:k_ft,1) = 0
                    ! prepare Fourier transform by copying in correct order
                    do ix = nx0,nxN
                        do iz = -nz,nz
                            ! copy arrays
                            VVdz(zf(iz), xf(ix), i_ft, 1) = gu(i,k)
                            VVdz(zf(iz), xf(ix), j_ft, 1) = gu(j,k)
                            VVdz(zf(iz), xf(ix), k_ft, 1) = u(k)
                        end do
                    end do
                    ! antitransform
                    do cntr = i_ft,k_ft
                        call IFT(VVdz(1:nzd,1:nxB,cntr,1))
                        call zTOx(VVdz(:,:,cntr,1), VVdx(:,:,cntr,1))
                        !call MPI_Alltoall(VVdz(:,:,cntr,1), 1, Mdz, VVdx(:,:,cntr,1), 1, Mdx, MPI_COMM_X) ! you could just copy without MPI since you deactivated xz parallelisation but ok
                        VVdx(nx+2:nxd+1,1:nzB,cntr,1)=0
                        call RFT(VVdx(1:nxd+1,1:nzB,cntr,1),rVVdx(1:2*nxd+2,1:nzB,cntr,1))
                    end do

                    ! compute product in real space
                    rVVdx(:,:,ik_ft,1) = rVVdx(:,:, i_ft, 1) * rVVdx(:,:, k_ft, 1) ! uk*( dui / dxk )
                    rVVdx(:,:,jk_ft,1) = rVVdx(:,:, j_ft, 1) * rVVdx(:,:, k_ft, 1) ! uk*( dui / dxk )

                    ! transform back
                    do cntr = ik_ft, jk_ft
                        call HFT(rVVdx(1:2*nxd+2,1:nzB,cntr,1), VVdx(1:nxd+1,1:nzB,cntr,1)); 
                        call xTOz(VVdx(:,:,cntr,1), VVdz(:,:,cntr,1))
                        !call MPI_Alltoall(VVdx(:,:,cntr,1), 1, Mdx, VVdz(:,:,cntr,1), 1, Mdz, MPI_COMM_X) ! you could just copy without MPI since you deactivated xz parallelisation but ok
                        call FFT(VVdz(1:nzd,1:nxB,cntr,1));
                        VVdz(1:nzd,1:nxB,cntr,1) = VVdz(1:nzd,1:nxB,cntr,1) * factor
                    end do

                    ! use Parseval's theorem to calculate statistics
                    do ix = nx0, nxN
                        c = 2; if (ix == 0) c = 1 ! multiplier for doubling points in x direction
                        do iz = -nz,nz
                            uiuj(ttrsp) = uiuj(ttrsp) - c * cprod( VVdz(zf(iz), xf(ix), ik_ft, 1), u(j) ) - c * cprod( VVdz(zf(iz), xf(ix), jk_ft, 1), u(i) )
                        end do
                    end do

                end do
                
            end do

        end do

    end do

    ! divide by nftot to obtain average
    uiujspectra = uiujspectra / nftot

    ! COMPUTE DERIVATIVES OF FLUX TERMS
    !----------------------------------

    do iz = -nz,nz
        do irs = 1,6

            ! var --> PHIvdiff
            call REALderiv(uiujspectra(irs, var, iz, :), uiujspectra(irs, PHIvdiff, iz, :))
            uiujspectra(irs, PHIvdiff, iz, :) = uiujspectra(irs, PHIvdiff, iz, :) * ni

            ! var --> vdiff
            call REALderiv2(uiujspectra(irs, var, iz, :), uiujspectra(irs, vdiff, iz, :))
            uiujspectra(irs, vdiff, iz, :) = uiujspectra(irs, vdiff, iz, :) * ni

            ! PHIptrsp --> ptrsp
            call REALderiv(uiujspectra(irs, PHIptrsp, iz, :), uiujspectra(irs, ptrsp, iz, :))

        end do
    end do

    ! COMPUTE PROFILES FROM SPECTRA
    !------------------------------

    ! if requested in input, multiply spectral budget times two
    if (timestwo) uiujspectra(:,var,:,:) = uiujspectra(:,var,:,:) * 2

    uiujprofiles = 0
    do iz = -nz,nz
        uiujprofiles(:,:,:) = uiujprofiles(:,:,:) + uiujspectra(:,:,iz,:)
    end do

    !---------------------------------------------------!
    !-----------------      WRITE      -----------------!
    !---------------------------------------------------!

    if (has_terminal) print *, "Saving to disk..."

    ! write metadata
    currfname = nfofile
    if (has_terminal) then
        open(15, file=currfname)
            write(15,*) "This is uiuj_spectra.f90."
            if (timestwo) write(15,*) "Correcting results by factor 2, to match older versions."
            write(15,*) ""
            write(15,*) "nmin", nmin 
            write(15,*) "nmax", nmax
            write(15,*) "dn", dn
            write(15,*) "nftot", nftot
            write(15,*) ""
            if (custom_mean) then
                write(15,*) "nmin_cm", nmin_cm
                write(15,*) "nmax_cm", nmax_cm
                write(15,*) "dn_cm", dn_cm
            end if
        close(15)
    end if

    ! write to disk
    currfname = "uiuj_spectra.bin"
    call MPI_File_open(MPI_COMM_WORLD, trim(currfname), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
        
        ! write mean data
        offset = 0
        CALL MPI_File_set_view(fh, offset, MPI_DOUBLE_PRECISION, mean_write_type, 'native', MPI_INFO_NULL)
        CALL MPI_File_write_all(fh, mean, 1, mean_inmem_type, MPI_STATUS_IGNORE)

        ! write uiuj spectral data
        offset = offset + (ny+3)*7*sizeof(mean(1,1)) ! DON'T DO SIZEOF(mean)! Program is PARALLEL!!!
        CALL MPI_File_set_view(fh, offset, MPI_DOUBLE_PRECISION, spectra_write_type, 'native', MPI_INFO_NULL)
        CALL MPI_File_write_all(fh, uiujspectra, 1, spectra_inmem_type, MPI_STATUS_IGNORE)

        ! write uiuj profiles data
        offset = offset + ( int(ny+3,MPI_OFFSET_KIND) * int((2*nz)+1,MPI_OFFSET_KIND) * 60_MPI_OFFSET_KIND * sizeof(uiujspectra(1,1,1,1)) ) ! DON'T DO SIZEOF(uiujspectra)! Program is PARALLEL!!!
        CALL MPI_File_set_view(fh, offset, MPI_DOUBLE_PRECISION, profiles_write_type, 'native', MPI_INFO_NULL)
        CALL MPI_File_write_all(fh, uiujprofiles, 1, profiles_inmem_type, MPI_STATUS_IGNORE)

        ! write convs data - NOT IMPLEMENTED YET!
        !offset = offset + ( int(ny+3,MPI_OFFSET_KIND) * 60_MPI_OFFSET_KIND * sizeof(uiujprofiles(1,1,1)) ) ! DON'T DO SIZEOF(uiujprofiles)! Program is PARALLEL!!!
        !CALL MPI_File_set_view(fh, offset, MPI_DOUBLE_PRECISION, convs_write_type, 'native', MPI_INFO_NULL)
        !CALL MPI_File_write_all(fh, convs, 1, convs_inmem_type, MPI_STATUS_IGNORE)

    call MPI_File_close(fh)

    if (has_terminal) call watermark()

    ! be polite and say goodbye
    if (has_terminal) print *, "Goodbye man!"

    !---------------------------------------------------!
    !-----------------     FINALISE    -----------------!
    !---------------------------------------------------!

    ! realease memory
    CALL free_fft()
    CALL free_memory(.FALSE.) 
    CALL MPI_Finalize()



contains !-------------------------------------------------------------------------------------------------------------------



    ! define and commit MPI filetypes
    subroutine init_uiuj_mpitypes()

        ! define type for reading pressure from disk
        CALL MPI_Type_create_subarray(3, [ny+3, 2*nz+1, nx+1], [nyN-ny0+5, 2*nz+1, nxB], [ny0-1,0,nx0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, press_read_type, ierror)
        CALL MPI_Type_commit(press_read_type, ierror)
        ! type describing pressure in memory
        CALL MPI_Type_create_subarray(3, [nyN-ny0+5, 2*nz+1, nxB], [nyN-ny0+5, 2*nz+1, nxB], [0,0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, press_field_type, ierror)
        CALL MPI_Type_commit(press_field_type, ierror)

        ! define type for writing mean on disk
        CALL MPI_Type_create_subarray(2, [7, ny+3], [7, maxy-miny+1], [0,miny+1], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mean_write_type, ierror)
        CALL MPI_Type_commit(mean_write_type, ierror)
        CALL MPI_Type_create_subarray(2, [7, nyN-ny0+5], [7, maxy-miny+1], [0,miny-(ny0-2)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mean_inmem_type, ierror)
        CALL MPI_Type_commit(mean_inmem_type, ierror)

        ! define type for writing uiujspectra on disk
        CALL MPI_Type_create_subarray(4, [6, 10, (2*nz)+1, ny+3], [6, 10, (2*nz)+1, maxy-miny+1], [0,0,0,miny+1], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, spectra_write_type, ierror)
        CALL MPI_Type_commit(spectra_write_type, ierror)
        CALL MPI_Type_create_subarray(4, [6, 10, (2*nz)+1, nyN-ny0+5], [6, 10, (2*nz)+1, maxy-miny+1], [0,0,0,miny-(ny0-2)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, spectra_inmem_type, ierror)
        CALL MPI_Type_commit(spectra_inmem_type, ierror)

        ! define type for writing uiujprofiles on disk
        CALL MPI_Type_create_subarray(3, [6, 10, ny+3], [6, 10, maxy-miny+1], [0,0,miny+1], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, profiles_write_type, ierror)
        CALL MPI_Type_commit(profiles_write_type, ierror)
        CALL MPI_Type_create_subarray(3, [6, 10, nyN-ny0+5], [6, 10,maxy-miny+1], [0,0,miny-(ny0-2)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, profiles_inmem_type, ierror)
        CALL MPI_Type_commit(profiles_inmem_type, ierror)

        ! define type for writing convs on disk
        CALL MPI_Type_create_subarray(4, [4, (2*nz)+1, nz+1, ny+3], [4, (2*nz)+1, nz+1, maxy-miny+1], [0,0,0,miny+1], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, convs_write_type, ierror)
        CALL MPI_Type_commit(convs_write_type, ierror)
        CALL MPI_Type_create_subarray(4, [4, (2*nz)+1, nz+1, nyN-ny0+5], [4, (2*nz)+1, nz+1, maxy-miny+1], [0,0,0,miny-(ny0-2)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, convs_inmem_type, ierror)
        CALL MPI_Type_commit(convs_inmem_type, ierror)

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



    subroutine get_gradient_33(R, grad)
        complex(C_DOUBLE_COMPLEX), intent(in) :: R(ny0-2:nyN+2,-nz:nz,nx0:nxN,3)
        complex(C_DOUBLE_COMPLEX), intent(out) :: grad(ny0-2:nyN+2,-nz:nz,nx0:nxN,3,3)

        ! y derivatives
        do iv = 1,3
            do ix = nx0, nxN
                do iz = -nz, nz
                    call COMPLEXderiv(R(:,iz,ix,iv), grad(:,iz,ix,iv,2))
                    call LeftLU5divStep2(D0mat, grad(:,iz,ix,iv,2))
                    grad(:,iz,ix,iv,1) = R(:,iz,ix,iv) * dcmplx(0, alfa0*ix)
                    grad(:,iz,ix,iv,3) = R(:,iz,ix,iv) * dcmplx(0, beta0*iz)
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
        r = dreal(dconjg(a)*b)
    end function cprod



    function zf(iz) result(iz_fft)
    integer, intent(in) :: iz
    integer :: iz_fft
        if (iz < 0) then
            iz_fft = nzd + 1 + iz
        else
            iz_fft = iz + 1
        end if
    end function zf



!----------------------------------------------------------------------------------------------------------------------------
! less useful stuff here
!----------------------------------------------------------------------------------------------------------------------------



    subroutine print_help()
        if (iproc == 0) then
            print *, "Calculates spectral TKE budget along the spanwise direction."
            print *, "Statistics are calculated on files ranging from index nfmin to nfmax with step dn. Usage:"
            print *, ""
            print *, "   mpirun [mpi args] uiuj_spectra [-h] nfmin nfmax dn [--custom_mean nmin_m nmax_m dn_m]"
            print *, ""
            print *, "If the flag --custom_mean is passed, the mean field is calculated on fields (nmin_m nmax_m dn_m); the remaining statistics are still calculated on (nfmin nfmax dn)."
            print *, ""
            print *, "This program is meant to be used on plane channels."
            print *, ""
            print *, "Results are output to uiuj_spectra.bin. Input (dns.in) is read from parent folder."
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
        end if
        call MPI_Finalize()
        stop
    end subroutine print_help



    ! read arguments
    subroutine parse_args()
        ii = 1
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
                case ('-two', '--timestwo') ! multiply uiuj budget by two
                ! to match output of old version
                    timestwo = .TRUE.
                case default
                    if (command_argument_count() < 3) then ! handle exception: no input
                        print *, 'ERROR: please provide one input file as command line argument.'
                        stop
                    end if
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



    subroutine check_in_files()
    logical :: inputs_missing = .FALSE.

        if (iproc == 0) then

            print *
            print *, "Checking input files..."

            currfname = "../dns.in"
            call check_da_file(inputs_missing)

            ! check fields
            if (.NOT. custom_mean) then
                do ii = nfmin, nfmax, dnf ! loop over files
                    ! check velocity
                    write(istring,*) ii; currfname = trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
                    call check_da_file(inputs_missing)
                    currfname = trim("../pField"//TRIM(ADJUSTL(istring))//".fld")
                    call check_da_file(inputs_missing)
                end do
            else ! if custom mean is used
                do ii = min(nmin,nmin_cm), max(nmax, nmax_cm), min(dn, dn_cm) ! loop over files
                    ! check velocity
                    write(istring,*) ii; currfname = trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
                    call check_da_file(inputs_missing)
                    currfname = trim("../pField"//TRIM(ADJUSTL(istring))//".fld")
                    call check_da_file(inputs_missing)
                end do
            end if

            ! return error if anything is missing
            if (inputs_missing) then
                print *
                print *, "ERROR: SOME INPUT FILES ARE MISSING"
                print *
                call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
            end if

        end if

        call MPI_BARRIER(MPI_COMM_WORLD)

    end subroutine
    


    subroutine check_da_file(inputs_missing)
    logical, intent(inout) :: inputs_missing
    logical :: file_exists
    ! name of file to check needs to be written to currfname
        inquire(file=trim(currfname), EXIST=file_exists)
        if (.NOT. file_exists) then
            print *, "Missing ", trim(currfname)
            inputs_missing = .TRUE.
        end if
    end subroutine



    subroutine watermark()
    integer :: datetime(8)
    character(40) :: daystring, monthstring, yearstring, hstr, mstr, sstr

        call date_and_time(values=datetime)
        write(daystring,'(I0.2)') datetime(3)
        write(monthstring,'(I0.2)') datetime(2)
        write(yearstring,'(I0.4)') datetime(1)
        write(hstr,'(I0.2)') datetime(5)
        write(mstr,'(I0.2)') datetime(6)
        write(sstr,'(I0.2)') datetime(7)
        open(15, file=nfofile,access='append')
            write(15,*)
            write(15,*) "EXECUTION COMPLETED ON ", trim(daystring), "/", trim(monthstring), "/", trim(yearstring), " at ", trim(hstr), ":", trim(mstr), ":", trim(sstr) 
        close(15)

    end subroutine



end program
