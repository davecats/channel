program print_time
use iso_c_binding, only: C_INT, C_DOUBLE
implicit none

integer :: i
character(len=32) :: arg
logical :: verbose = .FALSE.

    do i = 1, command_argument_count()

        call get_command_argument(i, arg)

        select case (arg)
        case ('-v', '--verbose')
            verbose = .TRUE.
        case ('-h', '--help')
            call print_help()
            stop
        case default
            call get_time(trim(arg))
        end select

    end do

    if (.NOT. verbose) then
        print *
    end if

contains

subroutine print_help()
    print *, "print_time [-v] file_1 [file_2 ... file_n]"
    print *, "Flag '-v' is used to obtain verbose output."
end subroutine

subroutine get_time(filename)
character(*), intent(in) :: filename
real(C_DOUBLE) :: time
integer :: position = 3*C_INT + 6*C_DOUBLE + 1

    open(unit=777, file=filename, access="stream", status="old", action="read")
        read(777,pos=position) time
    close(777)

    if (verbose) then
        write(*,'(a,$)') filename, " -----> "
    end if
    write(*,'(F11.6,$)') time
    write(*,'(a,$)') " "

    if (verbose) then
        print *
    end if

end subroutine 

end program