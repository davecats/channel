program resize
use iso_c_binding, only: C_INT, C_DOUBLE, C_DOUBLE_COMPLEX
implicit none

! input parsing
integer :: ii=1, imin, imax, io
character(len=32) :: arg
character(len=40) :: fname


! old files
complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:), allocatable :: read_old_v, old_v
integer :: old_x, old_y, old_z, nxzp, nzzp
real(C_DOUBLE) :: old_alfa0, old_beta0

! new files
complex(C_DOUBLE_COMPLEX), dimension(:,:,:,:), allocatable :: v
integer :: nx, ny, nz
real(C_DOUBLE) :: alfa0,beta0,ni,a,ymin,ymax,deltat,cflmax,time,dt_field,dt_save,t_max,gamma,meanpx,meanpz,meanflowx,meanflowz,u0,uN
logical :: time_from_restart, CPI
integer :: CPI_type, nstep, npy

! looping
integer :: ub_x, ub_z
integer :: ic, ix, iz, iy
integer :: ixr, ixl, izt, izb
real(C_DOUBLE) :: wr, wl, wt, wb


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

    ! allocate a zero-padded array for the old field, so that domain length in Fourier matches new domain length
    call get_weights(nz, beta0, old_beta0, izb, nzzp, wb, wt)
    call get_weights(nx, alfa0, old_alfa0, ixl, nxzp, wl, wr)
    allocate(old_v(-1:(ny+1), (-nzzp):(nzzp), 0:nxzp, 1:3))
    old_v = 0
    ub_x = min(nxzp,old_x); ub_z = min(nzzp, old_z)

    ! loop over files
    do ii=imin,imax
        
        ! read file
        call get_command_argument(ii, arg)
        read(arg, *) fname
        call read_old_field(fname)
        old_v(:,-ub_z:ub_z,0:ub_x,:) = read_old_v(:,-ub_z:ub_z,0:ub_x,:)

        ! copy into new array
        do ic = 1,3
            do ix = 0,nx
                call get_weights(ix, alfa0, old_alfa0, ixl, ixr, wl, wr)
                do iz = (-nz),nz
                    call get_weights(iz, beta0, old_beta0, izb, izt, wb, wt)
                    do iy = (-1), ny+1
                        v(iy,iz,ix,ic) = wb*wl*old_v(iy,izb,ixl,ic) + wb*wr*old_v(iy,izb,ixr,ic) + wt*wl*old_v(iy,izt,ixl,ic) + wt*wr*old_v(iy,izt,ixr,ic)
                    end do
                end do
            end do
        end do

        ! write new file
        fname = trim(fname)//'.new'
        open(unit=100,file=trim(fname),access="stream",action="write",iostat=io)
            write(100) nx,ny,nz
            write(100) alfa0,beta0,ni,a,ymin,ymax,time
            write(100) v
        close(100)

    end do

contains



subroutine print_help()
    print *
    print *, "resize [-h] file_1 [file_2 ... file_n]"
    print *
    print *, "Interpolates the velocity fiel of the input file(s) to a new x-z mesh;"
    print *, "properties of the new mesh are read from dns.in, while properties of the "
    print *, "old mesh are read from the input binary file. New and old ny must match."
    print *, 'If input is "file.name", output is "file.name.resized".'
    print *, 'If the desired dimensions are larger than the old ones, the field is'
    print *, 'zero-padded on high wavenumbers. Otherwise, high wavenumbers are truncated.'
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
        read(100) old_alfa0, old_beta0
    close(100)
   
    allocate(read_old_v(-1:(old_y+1),(-old_z):old_z,0:old_x,1:3))

end subroutine



subroutine read_old_field(filename)
character(len=40) :: filename

    open(unit=100,file=trim(filename),access="stream",status="old",action="read",iostat=io)
        read(100, pos=69), read_old_v
    close(100)

end subroutine



subroutine get_new_size()
character(len=40), parameter :: filename = 'dns.in'
    open(100, file=trim(filename))
        read(100, *) nx, ny, nz
        read(100, *) alfa0, beta0
        read(100, *) ni
        read(100, *) a, ymin, ymax
        read(100, *) CPI, CPI_type, gamma
        read(100, *) meanpx, meanpz
        read(100, *) meanflowx, meanflowz
        read(100, *) u0,uN
        read(100, *) deltat, cflmax, time
        read(100, *) dt_field, dt_save, t_max, time_from_restart
        read(100, *) nstep
        read(100, *) npy
    close(100)

    allocate(v(-1:(ny+1), (-nz):nz, 0:nx, 1:3))

end subroutine



subroutine get_weights(ii, lett, old_lett, im, ip, wm, wp)
integer, intent(in) :: ii
real(C_DOUBLE), intent(in) :: lett, old_lett
integer, intent(out) :: im, ip
real(C_DOUBLE), intent(out) :: wm, wp
real(C_DOUBLE) :: i_prj

    i_prj = ii * (lett/old_lett)
    im = floor(abs(i_prj)); ip = im + 1

    if (i_prj < 0) then
        ip = - im
        im = ip - 1
    end if

    wm = ((old_lett*ip)-(lett*ii))/old_lett
    wp = 1 - wm

end subroutine



end program