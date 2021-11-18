program resize
use iso_c_binding, only: C_INT, C_DOUBLE, C_DOUBLE_COMPLEX
implicit none

! input parsing
integer :: ii=1, imin, imax, io
character(len=32) :: arg
character(len=40) :: fname
integer*1, dimension(68) :: header


! old files
complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:), allocatable :: old_v
integer :: old_x, old_y, old_z

! new files
complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:), allocatable :: v
integer :: nx, ny, nz

! looping
integer :: ub_x, ub_z
integer :: ic, ix, iz, iy


    call parse_args()

    ! setup
    call get_command_argument(imin, arg)
    read(arg, *) fname
    call get_old_size(fname)
    call get_new_size()
    if (.NOT. ny==old_y) then
        print *, "ERROR: dimensions in y need to be the same."
        print *, "Old y dimension:", old_y
        print *, "Required y dimension:", ny
        stop
    end if
    ub_x = MIN(nx,old_x); ub_z = min(nz, old_z)

    ! loop over files
    do ii=imin,imax
        
        ! read file
        call get_command_argument(ii, arg)
        read(arg, *) fname
        call read_old_field(fname)

        ! copy into new array
        v = 0
        do ic = 1,3
            do ix = 0,ub_x
                do iz = (-ub_z),ub_z
                    do iy = (-1), ny+1
                        v(iy,iz,ix,ic) = old_v(iy,iz,ix,ic)
                    end do
                end do
            end do
        end do

        ! write new file
        fname = trim(fname)//'.new'
        open(unit=100,file=trim(fname),access="stream",action="write",iostat=io)
            write(100) header, v
        close(100)

    end do

contains



subroutine print_help()
    print *
    print *, "resize [-h] file_1 [file_2 ... file_n]"
    print *
    print *, "Resizes input file to desired dimensions, which are read from dns.in."
    print *, 'If input is "file.name", output is "file.name.resized".'
    print *, 'If the desired dimensions are larger than the old ones, the field is'
    print *, 'zero-padded on high wavenumbers. Otherwise, high wavenumbers are truncated.'
    print *, 'IMPORTANT: the desired number of points in y has to be the same as the old one.'
    print *, 'No interpolation is done. In no direction.'
end subroutine



subroutine parse_args()
    imin=command_argument_count(); imax=0
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
                if (ii < imin) imin=ii
                if (ii > imax) imax=ii
        end select
        ii = ii + 1
    end do
end subroutine parse_args



subroutine get_old_size(filename)
character(len=40) :: filename

    open(unit=100,file=trim(filename),access="stream",status="old",action="read",iostat=io)
        read(100) old_x, old_y, old_z
    close(100)

    ! read header
    open(unit=100,file=trim(filename),access="stream",status="old",action="read",iostat=io)
        read(100) header
    close(100)
    
    allocate(old_v(-1:(old_y+1),(-old_z):old_z,0:old_x,1:3))

end subroutine



subroutine read_old_field(filename)
character(len=40) :: filename

    open(unit=100,file=trim(filename),access="stream",status="old",action="read",iostat=io)
        read(100, pos=69), old_v
    close(100)

end subroutine



subroutine get_new_size()
character(len=40), parameter :: filename = 'dns.in'
    open(100, file=trim(filename))
        read(100, *) nx, ny, nz;
    close(100)

    allocate(v(-1:(ny+1), (-nz):nz, 0:nx, 1:3))

end subroutine


end program