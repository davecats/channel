program get_fake_psd

use dnsdata
implicit none

character(len=32) :: cmd_in_buf ! needed to parse arguments
character(len=32) :: arg

integer :: iv, ix, iz, iy, y_start
integer :: ierror, ii, nmin, nmax, dn
integer :: ndim, nxtot, nytot, nztot, nxloc, nyloc, nzloc
real(C_DOUBLE), dimension(:), allocatable :: mean
real(C_DOUBLE), allocatable :: fake_psd(:,:,:,:)

character(len = 40) :: istring, istring2, filename, barename
type(MPI_Datatype) :: wtype_3d, type_towrite 



!---------------------!
! PROGRAM BEGINS HERE !
!---------------------!

    call parse_args()
    
    ! Init MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call read_dnsin(.TRUE.) ! passing .TRUE. cause dns.in is in parent folder wrt cwd

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
    CALL init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB)
    call init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB,.FALSE.,[3,1])
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

    ndim  = 3 ! number of spatial dimension
    ! 1, 2, 3     ->  (u, v, w)
    ! 4, 5, 6, 7  ->  u'v', u'u', v'v', w'w'
    ! 8, 9, 10    ->  dudy, sqrt(dudy), dudy^2
    
    ! GLOBAL QUANTITIES, which is, overall sizes of stuff
    nxtot = nx+1 ! no. points
    nytot = ny+1 ! no. points
    nztot = 2*nz+1 ! no. points

    ! things that are local to each process
    ! number of points in each direction
    nxloc = nxtot
    nyloc = min(ny,maxy)-max(0,miny)+1 ! these refer to UNIQUE data! It does not consider halo cells, nor ghost cells
    nzloc = nztot

    !-------------------!
    ! END OF PARAMETERS !
    !-------------------!

    ! allocate stuff
    allocate(fake_psd(-nz:nz,nx0:nxN,ny0-2:nyN+2,1:3)) ! iz, ix, iy, iv (veloctiy)
    allocate(mean(ny0-2:nyN+2)) ! contains mean velocity in streamwise direction
    ! notice that mean includes ghost and halo cells!
    ! this makes it easy to get the wall-normal derivative with built in routines

    !---------------------------!
    ! MPI FILETYPES FOR WRITING !
    !---------------------------!

    ! each process calculates where its y chunk starts (zero based: MPI needs it so)
    y_start = max(0,miny)

    ! data owned by each process that needs to be written
    ! notice that this is "wrong"! but if you multiply by 3, you get the correct number of floats that you need to write.
    ! see the writing routine later.
    ! also keep in mind that there are no "holes" (halo/ghost cells) in cond_field
    CALL MPI_Type_create_subarray(3, [nzloc, nxloc, nyloc], [nzloc, nxloc, nyloc], [0, 0, 0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, type_towrite)
    CALL MPI_Type_commit(type_towrite, ierror)
    
    ! cond_avg - for writing (setting view)
    ! 0) number of dimensions of array you want to write
    ! 1) size along each dimension of the WHOLE array ON DISK
    ! 2) size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    ! 3) starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(4, [nztot, nxtot, nytot, ndim], [nzloc, nxloc, nyloc, ndim], [0, 0, y_start, 0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, wtype_3d)
    !                             0)            1)                        2)                    3)
    CALL MPI_Type_commit(wtype_3d, ierror)

    !----------------------------------!
    ! END OF MPI FILETYPES FOR WRITING !
    !----------------------------------!

    ! get the mean - fast!
    mean = 0
    do ii=nmin,nmax ! loop over files
        write(istring,*) ii
        filename=trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
        open(file=filename, unit=999, status="old", access="stream", action="read")
            do iy = ny0-2, nyN+2
                mean(iy) = mean(iy) + real(fieldmap(999, iy, 0, 0, 1)) ! U
            end do
        close(999)
    end do
    mean = mean / (nmax-nmin+1) ! divide by number of files

    ! empty fake_psd
    fake_psd = 0

    ! loop over files
    do ii=nmin,nmax

        ! read file
        write(istring,*) ii
        filename=trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
        call read_restart_file(filename, V)
        ! remove mean
        V(:,0,0,1) = V(:,0,0,1) - cmplx(mean(:))

        ! conditionally average current file; store in temp_fileavg
        do iy = miny, maxy ! loop for all y positions
            ! skip ghost cells
            if (iy < 0) cycle
            if (iy > ny) cycle
            ! fourier transform current y position
            ! result on rVVdx
            call fourier_antitransform(iy) 
            ! get square of velocity
            ! also, factor of Fourier transform is included now
            rVVdx(:,:,:,1) = rVVdx(:,:,:,1)*rVVdx(:,:,:,1) * factor
            ! transform back to Fourier domain
            call fourier_transform(iy)
            do iV=1,3
                do iz=-nz,nz
                    do ix=0,nx
                        fake_psd(iz,ix,iy,iV) = fake_psd(iz,ix,iy,iV) + dreal(dconjg(V(iy,iz,ix,iV))*V(iy,iz,ix,iV))
                    end do
                end do
            end do
        end do
        
    end do

    ! write to disk
    call write_out()

    ! be polite and say goodbye
    if (has_terminal) print *, "Done."

    ! realease memory
    CALL free_fft()
    CALL free_memory(.FALSE.)
    CALL MPI_Finalize()



contains !-------------------------------------------------------------------------------------------------------------------



    function convert_idx(i_read, n, center_v) result(i_write)
        integer, intent(in) :: i_read, n, center_v
        integer  :: i_write
        i_write = MODULO((3*n+1 + i_read - center_v),(2*n+1)) + 1
    end function



    subroutine fourier_antitransform(iy)

        integer, intent(in) :: iy

            ! do fourier transform
            ! WARNING: this differs from normal program, since there is no extension of the Fourier domain
            ! (or, no increase of spatial resolution to avoid dealiasing)
            ! i.e., no need to set anything zo zero in the arrays

            ! get velocity
            VVdz(nz+2:nzd-nz,1:nxB,1:3,1)=0
            VVdz(1:nz+1,1:nxB,1:3,1)=V(iy,0:nz,nx0:nxN,1:3)
            VVdz(nzd+1-nz:nzd,1:nxB,1:3,1)=V(iy,-nz:-1,nx0:nxN,1:3)
            
            do iV=1,3
                call IFT(VVdz(1:nzd,1:nxB,iV,1)) ! first transform in z

                ! Transpose for next Fourier transform
                ! WARNING: IF YOU WANT TO USE XZ PARALLELISATION, you need to use a MPI_Alltoall
                ! call similar to the one used in channel to do this matrix transposition.
                ! Thing is, you need to redefine the MPI types npx and npz.
                call MPI_Alltoall(VVdz(:,:,iV,1), 1, Mdz, VVdx(:,:,iV,1), 1, Mdx, MPI_COMM_X)

                call RFT(VVdx(:,:,iV,1),rVVdx(:,:,iV,1)) ! second transform in x

            end do
            ! rVVdx is now containing the antitransform
            ! its size is specified in init_fft

    end subroutine



    subroutine fourier_transform(iy)

        integer, intent(in) :: iy

        ! do fourier transform
        ! WARNING: this differs from normal program, since there is no extension of the Fourier domain
        ! (or, no increase of spatial resolution to avoid dealiasing)
        ! i.e., no need to set anything zo zero in the arrays
        do iV=1,3
            call HFT(rVVdx(1:2*nxd+2,1:nzB,iV,1),VVdx(1:nxd+1,1:nzB,iV,1)); 
            call MPI_Alltoall(VVdx(:,:,iV,1), 1, Mdx, VVdz(:,:,iV,1), 1, Mdz, MPI_COMM_X)
            call FFT(VVdz(1:nzd,1:nxB,iV,1));
        end do

        ! get velocity
        V(iy,0:nz,nx0:nxN,1:3) = VVdz(1:nz+1,1:nxB,1:3,1)
        V(iy,-nz:-1,nx0:nxN,1:3) = VVdz(nzd+1-nz:nzd,1:nxB,1:3,1)

    end subroutine



    subroutine write_out()

        integer(C_SIZE_T), parameter :: disp = 0

        type(MPI_File) :: fh
        TYPE(MPI_Status) :: status

        if (has_terminal) print *, "Writing fake psd..."
        barename = "fake_psd."

        write(istring,*) nmin
        write(istring2,*) nmax
        filename=trim(TRIM(ADJUSTL(barename))//TRIM(ADJUSTL(istring))//"_"//TRIM(ADJUSTL(istring2))//".bin")

        CALL MPI_File_open(MPI_COMM_WORLD, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
            CALL MPI_File_set_view(fh, disp, MPI_DOUBLE_PRECISION, wtype_3d, 'native', MPI_INFO_NULL)
            CALL MPI_File_write_all(fh, fake_psd, ndim, type_towrite, status)
        CALL MPI_File_close(fh)

        if (iproc == 0) then
            filename=trim(TRIM(ADJUSTL(barename))//TRIM(ADJUSTL(istring))//"_"//TRIM(ADJUSTL(istring2))//".nfo")
            open(15, file=filename)
                write(15,*) "This is fake_psd.f90."
                write(15,*) ""
                write(15,*) "nmin", nmin 
                write(15,*) "nmax", nmax
                write(15,*) "dn", dn
            close(15)
        end if

    end subroutine

    



!-------------------------------------------------------------------------------------------------------------------
! less useful stuff here
!-------------------------------------------------------------------------------------------------------------------



    subroutine print_help()
        print *, "Calculates 1D and 2D 'fake psd': <|FT(u_i u_i)|>"
        print *, "in the 1D case, FT is done only in the z direction and averaging is performed"
        print *, "also in the x direction. Otherwise, FT is done in the x,z directions and the"
        print *, "result is time-averaged."
        print *, "Syntax:"
        print *, "   fake_psd [-h] nmin nmax dn"
    end subroutine



    subroutine parse_args()
        ii = 1
        ! read arguments
        if (command_argument_count() < 1) then ! handle exception: no input
            print *, 'ERROR: please provide one input file as command line argument.'
            stop
        end if
        do while (ii <= command_argument_count()) ! parse optional arguments
            call get_command_argument(ii, arg)
            select case (arg)
                case ('-h', '--help') ! call help
                    call print_help()
                    stop
                case default
                    if (command_argument_count() < 3) then ! handle exception: no input
                        print *, 'ERROR: 5 positional arguments are needed.'
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
            end select
            ii = ii + 1
        end do

    end subroutine





end program






