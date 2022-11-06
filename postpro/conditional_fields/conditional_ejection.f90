program conditional_ejection

use dnsdata
implicit none

character(len=32) :: cmd_in_buf ! needed to parse arguments
character(len=32) :: arg

integer :: iv, ix, iz, iy, iy_ref, iproc_ref
logical :: has_ref
logical, dimension(:), allocatable :: gather_ref
logical, dimension(:,:), allocatable :: mask
integer :: ierror, ii, nmin, nmax, dn
integer :: ndim, nxtot, nytot, nztot, nxloc, nyloc, nzloc
real(C_DOUBLE), dimension(:,:,:), allocatable :: cond_avg, temp_fileavg, temp_cumul
integer :: it_cumul
integer, parameter :: max_it_cumul = 10 ! try decreasing this value if you get NaN in the output bin file!
integer*8 :: no_samples = 0

real(C_DOUBLE) :: y_ref, u_thr
character(len = 40) :: istring, istring2, filename

integer :: y_start = 0

type(MPI_Datatype) :: wtype_3d, mask_type, type_towrite ! , wtype_scalar

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
    nxtot = (2*nx) + 1 ! no. points
    nytot = (ny + 1) ! no. points
    nztot = (2*nz) + 1 ! no. points

    ! things that are local to each process
    ! number of points in each direction
    nxloc = nxtot
    nyloc = min(ny,maxy)-max(0,miny)+1 ! these refer to UNIQUE data! It does not consider halo cells, nor ghost cells
    nzloc = nztot

    !-------------------!
    ! END OF PARAMETERS !
    !-------------------!

    ! allocate stuff
    allocate(cond_avg(ndim, max(0,miny):min(ny,maxy), 1:nztot)) ! output - conditional average
    allocate(temp_fileavg(ndim, max(0,miny):min(ny,maxy), 1:nztot)) ! conditionally averaged field on a single velocity file
    allocate(temp_cumul(ndim, max(0,miny):min(ny,maxy), 1:nztot)) ! conditionally averaged field: intermediate cumulations
    allocate(mask(1:nxtot, 1:nztot)) ! tells all processes where the ejections are
    allocate(gather_ref(0:(nproc-1))) ! used to determine who has y_ref

    !---------------------------!
    ! MPI FILETYPES FOR WRITING !
    !---------------------------!

    ! each process calculates where its y chunk starts (zero based: MPI needs it so)
    y_start = max(0,miny)

    ! data owned by each process that needs to be written
    ! notice that this is "wrong"! but if you multiply by 3, you get the correct number of floats that you need to write.
    ! see the writing routine later.
    ! also keep in mind that there are no "holes" (halo/ghost cells) in cond_field
    CALL MPI_Type_create_subarray(2, [nyloc, nzloc], [nyloc, nzloc], [0, 0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, type_towrite)
    CALL MPI_Type_commit(type_towrite, ierror)

    ! describes mask for positions of ejections that gets broadcasted
    CALL MPI_Type_create_subarray(2, [nxloc, nzloc], [nxloc, nzloc], [0, 0], MPI_ORDER_FORTRAN, MPI_LOGICAL, mask_type)
    CALL MPI_Type_commit(mask_type, ierror)
    
    ! cond_avg - for writing (setting view)
    ! 0) number of dimensions of array you want to write
    ! 1) size along each dimension of the WHOLE array ON DISK
    ! 2) size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    ! 3) starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(3, [3, nytot, nztot], [3, nyloc, nzloc], [0, y_start, 0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, wtype_3d)
    !                             0)            1)                        2)                    3)
    CALL MPI_Type_commit(wtype_3d, ierror)

    !----------------------------------!
    ! END OF MPI FILETYPES FOR WRITING !
    !----------------------------------!

    ! determine which process has y_ref
    has_ref = .FALSE.
    call get_nearest_y(y_ref, iy_ref, has_ref)
    if (has_ref) then
        print *
        print *, "Found reference y position"
        print *, "iy_ref", iy_ref
        print *, "effective y_ref", y(iy_ref)
    end if

    ! let other processes know who has y_ref
    call MPI_Gather(has_ref,1,MPI_LOGICAL,gather_ref,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
    if (iproc==0) then
        do ii=0,(nproc-1)
            if (gather_ref(ii)) iproc_ref = ii
        end do
    end if
    call MPI_Bcast(iproc_ref,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

    ! empty cond_avg, it_cumul, temp_cumul
    cond_avg = 0
    it_cumul = 0
    temp_cumul = 0

    ! loop over files
    do ii=nmin,nmax

        ! read file
        write(istring,*) ii
        filename=trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
        call read_restart_file(filename, V)

        ! get mask and no_samples
        if (has_ref) then
            mask = .FALSE.
            no_samples = 0
            call fourier_antitransform(iy_ref)
            do iz=1,nztot
                do ix=1,nxtot
                    if (rVVdx(ix,iz,1,1) < u_thr) then
                        mask(ix,iz) = .TRUE.
                        no_samples = no_samples + 1
                    end if
                end do
            end do
        end if

        ! broadcast mask and no_samples for this file
        call MPI_Bcast(no_samples,1,MPI_INTEGER8,iproc_ref,MPI_COMM_WORLD,ierror)
        call MPI_Bcast(mask,1,mask_type,iproc_ref,MPI_COMM_WORLD,ierror)

        ! reset temp_fileavg
        temp_fileavg = 0

        ! conditionally average current file; store in temp_fileavg
        do iy = miny, maxy ! loop for all y positions
            ! skip ghost cells
            if (iy < 0) cycle
            if (iy > ny) cycle
            ! fourier transform current y position
            call fourier_antitransform(iy)
            ! go through each point and see if there is an ejection there
            do iz=1,nztot
                do ix=1,nxtot
                    if (mask(ix,iz)) then
                        call cumulate_on_temp(iz) ! if there is an ejection, use ejection as center and cumulate
                        ! notice that this subroutine divides by no_samples already
                    end if
                end do
            end do
        end do

        ! cumulate average of single file into actual array with average
        temp_cumul = temp_cumul + temp_fileavg
        it_cumul = it_cumul + 1
        
        ! if you reach max_it_cumul: cumulate temp_cumul into cond_avg
        if (it_cumul == 10) then
            temp_cumul = temp_cumul / (nmax - nmin + 1)
            cond_avg = cond_avg + temp_cumul
            ! reset the arrays
            temp_cumul = 0
            it_cumul = 0
        end if
        
    end do

    ! last cumulation (if you didn't reach max_it_cumul on the last round)
    ! if you did reach max_it_cumul, temp_cumul will actually be zero, so this section does nothing (as desired)
    temp_cumul = temp_cumul / (nmax - nmin + 1)
    cond_avg = cond_avg + temp_cumul

    ! write to disk
    call write_cond_field()

    ! be polite and say goodbye
    if (has_terminal) print *, "Done."

    ! realease memory
    CALL free_fft()
    CALL free_memory(.FALSE.)
    CALL MPI_Finalize()



contains !-------------------------------------------------------------------------------------------------------------------



    subroutine cumulate_on_temp(cz)
        integer, intent(in) :: cz
        integer :: cc, zz

        do cc = 1,3
            do zz=1,nztot
                temp_fileavg(cc,iy,convert_idx(zz,nz,cz)) = temp_fileavg(cc,iy,convert_idx(zz,nz,cz)) + rVVdx(ix,zz,cc,1)/no_samples
            end do
        end do

    end subroutine



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

    end subroutine



    subroutine write_cond_field() 

        integer(C_SIZE_T), parameter :: disp = 0

        type(MPI_File) :: fh
        TYPE(MPI_Status) :: status

        if (has_terminal) print *, "Writing 2D conditional field..."

        write(istring,*) nmin
        write(istring2,*) nmax
        filename=trim("cond_field_ejection."//TRIM(ADJUSTL(istring))//"_"//TRIM(ADJUSTL(istring2))//".bin")

        CALL MPI_File_open(MPI_COMM_WORLD, filename, IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
            CALL MPI_File_set_view(fh, disp, MPI_DOUBLE_PRECISION, wtype_3d, 'native', MPI_INFO_NULL)
            CALL MPI_File_write_all(fh, cond_avg, 3, type_towrite, status)
        CALL MPI_File_close(fh)

        if (has_ref) then
            filename=trim("cond_field_ejection."//TRIM(ADJUSTL(istring))//"_"//TRIM(ADJUSTL(istring2))//".nfo")
            open(15, file=filename)
                write(15,*) "This is conditional_ejection.f90."
                write(15,*) ""
                write(15,*) "nmin", nmin 
                write(15,*) "nmax", nmax
                write(15,*) "dn", dn
                write(15,*) "" 
                write(15,*) "y_ref", y_ref
                write(15,*) "effective y_ref", y(iy_ref)
                write(15,*) "u_thr", u_thr
            close(15)
        end if

    end subroutine



    subroutine get_nearest_y(reqy, iyy, has_yy)
        integer, intent(out) :: iyy
        logical, intent(out) :: has_yy
        real(C_DOUBLE), intent(in) :: reqy
        real(C_DOUBLE) :: cdiff
        integer :: ii
    
            iyy = ny0-2; cdiff = abs(reqy - y(iyy))
            do ii = ny0-2, nyN+2
                if (abs(reqy - y(ii)) < cdiff) then
                    iyy = ii
                    cdiff = abs(reqy - y(iyy))
                end if
            end do
    
            if ( (iyy >= miny) .AND. (iyy <= maxy) ) then
                has_yy = .TRUE.
            end if
            
    end subroutine

    




!-------------------------------------------------------------------------------------------------------------------
! less useful stuff here
!-------------------------------------------------------------------------------------------------------------------



    subroutine print_help()
        print *, "Calculates 2D (spanwise - wall-normal) conditional field centered around"
        print *, "a sweeping event."
        print *, "Syntax:"
        print *, "   conditional_ejection [-h] nmin nmax dn y_ref u_thr"
        print *, "The sweeping event is defined as u(x,y_ref,z,t) < u_thr, where u is the"
        print *, "streamwise velocity, x,z describe the stream- and span-wise positions and"
        print *, "t represents the time. The field is effectively conditionally averaged in"
        print *, "x and t."
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
                    if (command_argument_count() < 5) then ! handle exception: no input
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
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) y_ref
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) u_thr
            end select
            ii = ii + 1
        end do

    end subroutine



end program






