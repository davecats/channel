! Scalar SCAL was never assigned any value
! one cool thing would be to parallelise the file format as well
! check out: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

! also, it might be nice to allow for mutltiple fields at once



program out2vtk
! convert out file to VTK (XML) format for paraview
! syntax:
! mpirun -np 1 out2vtk Dati.cart.xx.out

use dnsdata
implicit none

character(len=32) :: cmd_in_buf ! needed to parse arguments

character(len=40) :: in_fname, out_fname ! input and output files

integer :: n_nx, n_nz
real*8 :: na0, nb0

real*8 :: cap_x, cap_z

integer :: ix, iz, iy
integer :: ierror, ii=1
complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:), allocatable :: vec ! watch out for this type! must be the same as 

type(MPI_Datatype) :: type_w2disk, type_towrite ! , wtype_scalar



    ! read arguments
    if (command_argument_count() /= 6 .AND. command_argument_count() /= 1) then ! handle exception: no input
        print *, 'ERROR: please provide arguments as specified by help screen (xz_resize -h).'
        stop
    end if
    call get_command_argument(ii, cmd_in_buf)
    select case (cmd_in_buf)
        case ('-h', '--help') ! call help
            call print_help()
            stop
        case default
            read(cmd_in_buf, *) in_fname
            ii = ii+1
            call get_command_argument(ii, cmd_in_buf)
            read(cmd_in_buf, *) out_fname
            ii = ii+1
            call get_command_argument(ii, cmd_in_buf)
            read(cmd_in_buf, *) n_nx
            ii = ii+1
            call get_command_argument(ii, cmd_in_buf)
            read(cmd_in_buf, *) n_nz
            ii = ii+1
            call get_command_argument(ii, cmd_in_buf)
            read(cmd_in_buf, *) na0
            ii = ii+1
            call get_command_argument(ii, cmd_in_buf)
            read(cmd_in_buf, *) nb0
    end select  
    
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

    ! allocate stuff
    allocate(vec(ny0-2:nyN+2,-n_nz:n_nz,0:n_nx,1:3)); vec=0 ! this mimicks V from dnsdata, just using new parameters

    ! determine bounds for interpolation
    cap_x = min(n_nx, floor(alfa0 * nx / na0))
    cap_z = min(n_nz, floor(beta0 + nz / nb0))



    !---------------------------!
    ! MPI FILETYPES FOR WRITING !
    !---------------------------!

    ! data owned by each process that needs to be written
    CALL MPI_Type_create_subarray(4, [nyN-ny0+5, 2*n_nz+1, n_nx+1, 3], [maxy-miny+1, 2*n_nz+1, n_nx+1, 3], [miny-(ny0-2), 0, 0, 0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, type_towrite)
    CALL MPI_Type_commit(type_towrite, ierror)
    
    ! vec - for writing (setting view)
    CALL MPI_Type_create_subarray(4, [ny+3,2*n_nz+1,n_nx+1,3], [maxy-miny+1, 2*n_nz+1, n_nx+1, 3], [miny+1, 0, 0, 0], MPI_ORDER_FORTRAN, MPI_REAL, type_w2disk)
    CALL MPI_Type_commit(type_w2disk, ierror)

    !----------------------------------!
    ! END OF MPI FILETYPES FOR WRITING !
    !----------------------------------!



    ! read file
    call read_restart_file(in_fname,V)

    ! interpolate
    do ii = 1,3
        do ix=0,cap_x
            do iz=-cap_z,cap_z
                do iy = miny,maxy ! skips halo cells: only non-duplicated data
                    call xz_interpolation(iy,iz,ix,ii)
                end do
            end do
        end do
    end do
    
    if (has_terminal) print *, "Saving to disk..."
    call write_resized(out_fname,vec)
    if (has_terminal) print *, "... done!"

    ! realease memory
    CALL free_fft()
    CALL free_memory(.FALSE.) 
    CALL MPI_Finalize()



contains !-------------------------------------------------------------------------------------------------------------------



    SUBROUTINE write_resized(filename,R)
    complex(C_DOUBLE_COMPLEX), intent(in) :: R(ny0-2:nyN+2,-n_nz:n_nz,0:n_nx,1:3)
    character(len=40), intent(in) :: filename
    ! mpi stuff
    TYPE(MPI_File) :: fh
    INTEGER(MPI_OFFSET_KIND) :: disp 
    TYPE(MPI_Status) :: status
    
        ! open file
        CALL MPI_File_open(MPI_COMM_WORLD, TRIM(filename), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh) 

            ! write header
            IF (has_terminal) THEN ! only one process does this
              CALL MPI_file_write(fh, [n_nx,ny,n_nz], 3, MPI_INTEGER, status)
              CALL MPI_file_write(fh, [na0,nb0,ni,a,ymin,ymax,time], 7, MPI_DOUBLE_PRECISION, status)
            END IF

            ! set view to subarray
            disp = 3*C_INT + 7*C_DOUBLE ! offset to skip header
            CALL MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, type_w2disk, 'native', MPI_INFO_NULL)

            ! finally write field
            CALL MPI_File_write_all(fh, R, 1, type_towrite, status)

        call MPI_File_close(fh)

    END SUBROUTINE write_resized



    subroutine xz_interpolation(iy,zz,xx,ii) ! this interpolates the field on a smaller grid if undersampling is requested
    integer, intent(in) :: xx, zz, iy, ii
    real*8 :: xprj, zprj
    integer :: xc, xf, zc, zf
    complex(C_DOUBLE_COMPLEX) :: u_cc, u_cf, u_fc, u_ff ! shortcuts for values; first letter of pedix refers to x, second to z
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
        u_cc = V(iy,zc,xc,ii);   w_cc = w_xc * w_zc
        u_ff = V(iy,zf,xf,ii);   w_ff = w_xf * w_zf
        u_cf = V(iy,zf,xc,ii);   w_cf = w_xc * w_zf
        u_fc = V(iy,zc,xf,ii);   w_fc = w_xf * w_zc

        ! do interpolation
        vec(iy,zz,xx,ii) = u_cc*w_cc + u_ff*w_ff + u_fc*w_fc + u_cf*w_cf

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

        xprj = xx*na0/alfa0
        zprj = zz*nb0/beta0

    end subroutine undersampled_to_fullindex



!-------------------------------------------------------------------------------------------------------------------
! less useful stuff here
!-------------------------------------------------------------------------------------------------------------------



    subroutine print_help()
        print *, "Resizes a given .out file from channel.f90 to match provided mesh"
        print *, "details in the x and z directions. Details on the old mesh are read"
        print *, "from dns.in."
        print *, "   xz_resize [-h]  old_file.input  name_of_new_file.output"
        print *, "                   new_nx  new_nz  new_alfa0  new_beta0"
        print *, "The new number of points and the Fourier resolutions in x and z are"
        print *, "passed as arguments."
    end subroutine



end program






