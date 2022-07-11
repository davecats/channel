! Scalar SCAL was never assigned any value
! one cool thing would be to parallelise the file format as well
! check out: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

! also, it might be nice to allow for mutltiple fields at once



program out2bin
! convert out file to VTK (XML) format for paraview
! syntax:
! mpirun -np 1 out2vtk Dati.cart.xx.out

use dnsdata
implicit none

character(len=32) :: cmd_in_buf ! needed to parse arguments
logical :: fluct_only = .FALSE., get_large = .FALSE., get_small = .FALSE., undersample = .FALSE.
real*4 :: xcent_undersample = 1
character(len=32) :: arg

integer :: iv, ix, iz, iy
integer :: ierror, ii=1, imin, imax
integer :: ndim, nxtot, nytot, nztot, nxloc, nyloc, nzloc, z_threshold
integer, parameter :: large = 1, small = 2, default = 0
real*8, dimension(:,:,:,:), allocatable :: vec ! watch out for this type! must be the same as 

integer :: y_start = 0

type(MPI_Datatype) :: wtype_3d, type_towrite ! , wtype_scalar

    call parse_args()
    
    ! Init MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call read_dnsin()

    !-----------!
    ! ATTENTION !
    !-----------!
    ! set npy to total number of processes
    ! this effectively DEACTIVATES PARALLELISATION IN Z
    ! and makes life easier (without significant losses in perf.)
    npy = nproc
    ! notice that this overrides value from dns.in
    !-----------------------------------------------------------------!
    ! NEVER RUN THIS PROGRAM WITH X/Z PARALLELISATION (it won't work) !
    !-----------------------------------------------------------------!

    call init_MPI(nx+1,nz,ny,nxd+1,nzd)
    call init_memory(.FALSE.)
    call init_fft(VVdz,VVdx,rVVdx,nxB,nxB,2*nz+1,2*nz+1,.TRUE.,[3,1])
    ! LAST FLAG OF init_fft (.TRUE.) SPECIFIES THAT REAL FFT TRANSFORM HAS AN ODD LOGICAL SIZE
    ! Notice that in the call to init_fft the values of nxd and nzd have been replaced by nxB and 2*nz+1.
    ! Since the parallelisation in xz is deactivated, nxB=nx+1.
    ! Instead, nxd and nzd are the expanded number of points that are needed for dealiasing;
    ! here, there is no problem with aliasing, so no need to increase the spatial-resolution before
    ! antitransforming. Therefore, the arguments corresponding to the increased number of points for dealiasing
    ! are replaced with the actual, regular number of points.
    ! Also, since nzB is by definition the same as nzd, its value has been replaced with 2*nz+1.

    ! WARNING: NEVER TRUST nz0, nzN, nzB: these are all dependend on nzd!
    ! Instead, nx0, nxN, nxB can be used (as they are defined on nx+1).



    !------------!
    ! PARAMETERS !
    !------------!

    ndim  = 3  ! number of spatial dimension

    ! GLOBAL QUANTITIES, which is, overall sizes of stuff
    nxtot = nint(xcent_undersample * ((2*nx) + 1)) ! no. points
    nytot = (ny + 1) ! no. points
    nztot = nint(xcent_undersample * (2*nz + 1)) ! no. points

    ! things that are local to each process
    ! number of points in each direction
    nxloc = nxtot
    nyloc = min(ny,maxy)-max(0,miny)+1 ! these refer to UNIQUE data! It does not consider halo cells, nor ghost cells
    nzloc = nztot

    !-------------------!
    ! END OF PARAMETERS !
    !-------------------!



    ! allocate stuff
    allocate(vec(ndim, 1:nxtot, max(0,miny):min(ny,maxy), 1:nztot))    ! vectorial field



    !---------------------------!
    ! MPI FILETYPES FOR WRITING !
    !---------------------------!

    ! each process calculates where its y chunk starts (zero based: MPI needs it so)
    y_start = max(0,miny)

    ! data owned by each process that needs to be written
    CALL MPI_Type_create_subarray(3, [nxloc, nyloc, nzloc], [nxloc, nyloc, nzloc], [0, 0, 0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, type_towrite)
    CALL MPI_Type_commit(type_towrite, ierror)
    
    ! XYZ AND VEC - for writing (setting view)
    ! 0) number of dimensions of array you want to write
    ! 1) size along each dimension of the WHOLE array ON DISK
    ! 2) size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    ! 3) starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(4, [3, nxtot, nytot, nztot], [3, nxloc, nyloc, nzloc], [0, 0, y_start, 0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, wtype_3d)
    !                             0)            1)                        2)                    3)
    CALL MPI_Type_commit(wtype_3d, ierror)

    !----------------------------------!
    ! END OF MPI FILETYPES FOR WRITING !
    !----------------------------------!

    ! loop over files
    do ii=imin,imax
        call convert_outfile(default)
        if (get_large) call convert_outfile(large)
        if (get_small) call convert_outfile(small)
    end do

    call write_metadata()

    ! realease memory
    CALL free_fft()
    CALL free_memory(.FALSE.)
    CALL MPI_Finalize()



contains !-------------------------------------------------------------------------------------------------------------------



    subroutine convert_outfile(ls)

        integer, intent(in) :: ls

        ! read file
        call get_command_argument(ii, cmd_in_buf)
        read(cmd_in_buf, *) fname
        CALL read_restart_file(fname,V)

        ! remove average if necessary
        if (fluct_only) then
            V(:,0,0,:) = 0
        end if

        ! filter (if necessary)
        if (.NOT. ls == default) call filter_vfield(ls)

        ! for each y plane
        do iy = miny,maxy ! skips halo cells: only non-duplicated data

            ! skip ghost cells
            if (iy < 0) cycle
            if (iy > ny) cycle

            ! do fourier transform
            ! WARNING: this differs from normal program, since there is no extension of the Fourier domain
            ! (or, no increase of spatial resolution to avoid dealiasing)
            ! i.e., no need to set anything zo zero in the arrays
            VVdz(1:nz+1,:,1:3,1)=V(iy,0:nz,:,1:3);
            VVdz(nz+2:2*nz+1,:,1:3,1)=V(iy,-nz:-1,:,1:3); 
            do iV=1,3
                call IFT(VVdz(1:2*nz+1,1:nxB,iV,1)) ! first transform in z

                ! Transpose for next Fourier transform
                ! WARNING: IF YOU WANT TO USE XZ PARALLELISATION, you need to use a MPI_Alltoall
                ! call similar to the one used in channel to do this matrix transposition.
                ! Thing is, you need to redefine the MPI types npx and npz.
                do ix = 1, nxB
                    do iz = 1, 2*nz + 1
                        VVdx(ix,iz,iV,1) = VVdz(iz,ix,iV,1)
                    end do
                end do

                call RFT(VVdx(:,:,iV,1),rVVdx(:,:,iV,1)) ! second transform in x

            end do
            ! rVVdx is now containing the antitransform
            ! it has size 2*nx+1,2*nz+1,6,6

            ! convert velocity vector for each process
            do iV = 1,3
                do iz=1,nztot
                    do ix=1,nxtot
                        call xz_interpolation(iV,ix,iy,iz)
                    end do
                end do
            end do

        end do

        call write_brutal_bin(fname,ls)

    end subroutine



    subroutine write_brutal_bin(fname,ls) 

        character(len = 40) :: fname
        integer(C_SIZE_T), parameter :: disp = 0
        integer, intent(in) :: ls
        character(len = 40) :: interstring = ""

        integer   :: dotpos

        type(MPI_File) :: fh
        TYPE(MPI_Status) :: status

        ! adapt filename
        dotpos = scan(trim(fname),".", BACK= .true.)
        if (ls == large) interstring = "large."
        if (ls == small) interstring = "small."
        if ( dotpos > 0 ) fname = fname(1:dotpos)//TRIM(interstring)//"bin"

        if (has_terminal) print *, "Writing "//TRIM(fname)

        CALL MPI_File_open(MPI_COMM_WORLD, TRIM(fname), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
            CALL MPI_File_set_view(fh, disp, MPI_DOUBLE_PRECISION, wtype_3d, 'native', MPI_INFO_NULL)
            CALL MPI_File_write_all(fh, vec, 3, type_towrite, status)
        CALL MPI_File_close(fh)

    end subroutine



    subroutine xz_interpolation(ii,xx,iy,zz) ! this interpolates the field on a smaller grid if undersampling is requested
    integer, intent(in) :: xx, zz, iy, ii
    real*8 :: xprj, zprj
    integer :: xc, xf, zc, zf
    real*8 :: u_cc, u_cf, u_fc, u_ff ! shortcuts for values; first letter of pedix refers to x, second to z
    real*8 :: w_cc, w_cf, w_fc, w_ff, w_xc, w_zc, w_xf, w_zf ! weights for values above
    
        ! first off, project indeces so that maximum range is (2*nx+1), (2*nz+1)
        call undersampled_to_fullindex(xx,zz, xprj,zprj)

        ! find 4 nearest points
        xc = ceiling(xprj); zc = ceiling(zprj)
        xf = floor(xprj); zf = floor(zprj)

        ! find weights in x
        if (xc == xf) then
            w_xc = 0.5; w_xf = 0.5
        else
            w_xc = abs(real(xf) - xprj) / abs(xc-xf) ! the xf there is correct! the closer xprj to xc, the bigger
            w_xf = abs(real(xc) - xprj) / abs(xc-xf) ! the xc there is correct! the closer xprj to xf, the bigger
        end if

        ! find weights in z
        if (zc == zf) then
            w_zc = 0.5; w_zf = 0.5
        else
            w_zc = abs(real(zf) - zprj) / abs(zc-zf) ! the zf there is correct! the closer zprj to zc, the bigger
            w_zf = abs(real(zc) - zprj) / abs(zc-zf) ! the zc there is correct! the closer zprj to zf, the bigger
        end if

        ! shortcuts for values and weights
        u_cc = rVVdx(xc,zc,ii,1);   w_cc = w_xc * w_zc
        u_ff = rVVdx(xf,zf,ii,1);   w_ff = w_xf * w_zf
        u_cf = rVVdx(xc,zf,ii,1);   w_cf = w_xc * w_zf
        u_fc = rVVdx(xf,zc,ii,1);   w_fc = w_xf * w_zc

        ! do interpolation
        vec(ii,ix,iy,iz) = u_cc*w_cc + u_ff*w_ff + u_fc*w_fc + u_cf*w_cf

    end subroutine



    subroutine undersampled_to_fullindex(xx,zz, xprj,zprj)
    ! Takes indeces xx, zz from the undersampled domain and returns their index-position in the full domain.
    ! The full domain is what is actually stored in memory (rVVdx); so, xprj and zprj *could* be indeces of
    ! the rVVdx array, if they were integers. Instead, xx and zz are indeces of the array vec that is being
    ! written to memory.
    ! Notice that xx=1,nxtot and zz=1,nztot,
    ! whereas xprj=1,2*nx+1 and zprj=1,2*nz+1.
    ! However, xprj and zprj are REAL NUMBERS, meaning that (xx,zz) can correspond to a position that is
    ! inbetween two points of the array rVVdx.
    ! The logic here is that:
    ! - xx=1 gets mapped to xprj=1 (as this is x=0 in space, i.e. the first extreme of the periodic domain)
    ! - xx=nxtot+1 gets mapped to xprj=2*nx+2 (as this would be the other extreme of the periodic domain)
    ! Same applies for z.

    integer, intent(in) :: xx, zz
    real*8, intent(out) :: xprj,zprj

        xprj = real(2*nx+1)/real(nxtot) * (xx - 1.0) + 1
        zprj = real(2*nz+1)/real(nztot) * (zz - 1.0) + 1

    end subroutine undersampled_to_fullindex



    subroutine write_metadata()
    real*8 :: xp, zp, dx, dz

        ! get resolution in x and z directions
        call undersampled_to_fullindex(2,2, dx,dz)
        call undersampled_to_fullindex(1,1, xp,zp)

        dx = dx - xp; dx = abs(dx); dx = dx * 2 * PI / alfa0 / (2*nx + 1)
        dz = dz - zp; dz = abs(dz); dz = dz * 2 * PI / beta0 / (2*nz + 1)

        if (has_terminal) then
            open(1977, file="spatial_grid.nfo")
                write(1977,*) nxtot, nztot, dx, dz
                write(1977,*) ny, ny+1, a, ymin, ymax
                write(1977,*) "# * ~ - ~ * # LEGEND # * ~ - ~ * #"
                write(1977,*) "nxtot, nztot, dx, dz"
                write(1977,*) "ny (see dns.in), nytot, a, ymin, ymax, nztot, dx, dz"
            close(1977)
        end if

    end subroutine



    subroutine largesmall_setup()
        open(15, file='largesmall_settings.in')
            read(15, *) z_threshold
        close(15)
        if (has_terminal) then
            write(*,"(A,I5)") "   z_threshold =", z_threshold
        end if
    end subroutine largesmall_setup



    subroutine filter_vfield(ls)
    integer, intent(in) :: ls

        do iz = -nz,nz
            if ( (abs(beta0*iz) <= z_threshold .AND. ls == small) .OR. (abs(beta0*iz) > z_threshold .AND. ls == large ) ) then
                V(:,iz,:,:) = 0
            end if
        end do

    end subroutine
    




!-------------------------------------------------------------------------------------------------------------------
! less useful stuff here
!-------------------------------------------------------------------------------------------------------------------



    subroutine print_help()
        print *, "Converts a given .out file from channel into a pure binary velocity array."
        print *, "The binary array is Fourier-antitransformed with respect to the .out file,"
        print *, "so that it represents a velocity field in space (and not in the Fourier domain)."
        print *, "Syntax:"
        print *, "   out2bin [-h] [-f] [-l -s -ls -sl] [-u undersample_xcent] file.name"
        print *, "If flag '-f' (or '--fluctuation') is passed, the (spatial) average is"
        print *, "subtracted from the field before converting."
        print *, "If flags -l, -s or so are passed, filtering thresholds are read from largesmall.in"
        print *, "and the desired large or small scale field is returned."
        print *, "Flag -u is used to undersample the field. The the number n_new of points in the x and z"
        print *, "directions of the output files is given by n_new = usample_xcent*n_old, where n_old is"
        print *, "the number of points of the input."
    end subroutine



    subroutine parse_args()
        
        ! read arguments
        imin=command_argument_count(); imax=0
        if (command_argument_count() < 1) then ! handle exception: no input
            print *, 'ERROR: please provide one input file as command line argument.'
            stop
        end if
        do while (ii <= command_argument_count()) ! parse optional arguments
            call get_command_argument(ii, arg)
            select case (arg)
                case ('-f', '--fluctuation') ! flag to only have fluctuations
                    fluct_only = .TRUE.
                case ('-l', '--large')
                    get_large = .TRUE.
                    call largesmall_setup()
                case ('-s', '--small')
                    get_small = .TRUE.
                    call largesmall_setup()
                case ('-ls', '-sl', '--largesmall', '--smalllarge')
                    get_small = .TRUE.
                    get_large = .TRUE.
                    call largesmall_setup()
                case ('-u', '--undersample') ! specify undersampling
                    undersample = .TRUE.
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) xcent_undersample
                case ('-h', '--help') ! call help
                    call print_help()
                    stop
                case default
                    if (ii < imin) imin=ii
                    if (ii > imax) imax=ii
            end select
            ii = ii + 1
        end do

    end subroutine



end program






