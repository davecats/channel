program print_time
use iso_c_binding, only: C_INT, C_DOUBLE
implicit none

integer :: i
character(len=32) :: arg
logical :: verbose = .FALSE.

    do i = 1, command_argument_count()

        call get_command_argument(i, arg)

        select case (arg)
        case('-v', '--verbose')
            verbose = .TRUE.
        case ('-h', '--help')
            call print_help()
            stop
        case default
            call get_header(trim(arg))
        end select

    end do

    if (.NOT. verbose) then
        print *
    end if

contains

subroutine print_help()
    print *, "print_header file"
    print *, "Prints header of a Dati.cart.xx.out file."
end subroutine

subroutine get_header(filename)
character(*), intent(in) :: filename
integer(C_INT) :: nx, ny, nz
real(C_DOUBLE) :: alfa0, beta0, ni, a, ymin, ymax, time
integer :: position = 3*C_INT + 6*C_DOUBLE + 1

    open(unit=777, file=filename, access="stream", status="old", action="read")
        
        read(777) nx, ny, nz
        read(777) alfa0, beta0, ni, a, ymin, ymax, time
    
    close(777)

    ! print integers
    if (verbose) then
        write(*,'(a,$)') "nx, ny, nz -----> "
    end if
    write(*,'(I5,$)') nx, ny, nz
    write(*,'(a,$)') " "
    if (verbose) then
        print *
    end if

    ! print frequency resolution
    if (verbose) then
        write(*,'(a,$)') "alfa0, beta0 -----> "
    end if
    write(*,'(F11.6,$)') alfa0, beta0
    write(*,'(a,$)') " "
    if (verbose) then
        print *
    end if

    ! print reynolds number
    if (verbose) then
        write(*,'(a,$)') "ni -----> "
    end if
    write(*,'(F11.6,$)') ni
    write(*,'(a,$)') " "
    if (verbose) then
        print *
    end if
    
    ! print y discretisation
    if (verbose) then
        write(*,'(a,$)') "a, ymin, ymax -----> "
    end if
    write(*,'(F11.6,$)') a, ymin, ymax
    write(*,'(a,$)') " "
    if (verbose) then
        print *
    end if

    ! print time
    if (verbose) then
        write(*,'(a,$)') "time -----> "
    end if
    write(*,'(F11.6,$)') time
    write(*,'(a,$)') " "
    if (verbose) then
        print *
    end if

end subroutine 

end program