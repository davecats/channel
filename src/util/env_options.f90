#include "build_options.h"

! Uniform access to the CHANNEL_* environment overrides.
!
! Every option used to parse itself at its point of use, which let the accepted
! spellings drift apart -- CHANNEL_OVERLAPPING understood "off" while
! CHANNEL_Y_PIPELINE_TIMING did not, and CHANNEL_YS_FORCE_CUSTOM_GPSV wanted an
! integer where its neighbours wanted a word.  Routing every lookup through this
! module keeps a single answer to "what counts as true" and to "what is a valid
! integer", and leaves the call sites saying only what they want.
module env_options

  use, intrinsic :: iso_c_binding, only: C_INT, C_INT64_T
  implicit none
  private

  integer, parameter :: ENV_VALUE_LEN = 256

  public :: env_text, env_flag, env_int, env_int64, env_int_list, lowercase

contains

  ! Raw value of <name>, trimmed.  .false. when unset or empty, in which case
  ! value is untouched.
  logical function env_text(name, value)
    character(len=*), intent(in) :: name
    character(len=*), intent(out) :: value
    character(len=ENV_VALUE_LEN) :: raw
    integer :: length, status

    env_text = .false.
    value = ""
    call get_environment_variable(name, raw, length, status)
    if (status /= 0 .or. length <= 0) return

    value = adjustl(trim(raw(:min(length, ENV_VALUE_LEN))))
    env_text = (len_trim(value) > 0)
  end function env_text

  ! Updates value when <name> holds a recognised boolean.  An unset variable
  ! leaves value alone; an unrecognised one warns and also leaves it alone, so
  ! callers can simply assign their default first and then call this.
  subroutine env_flag(name, value)
    character(len=*), intent(in) :: name
    logical, intent(inout) :: value
    character(len=ENV_VALUE_LEN) :: text
    integer :: parsed, io

    if (.not. env_text(name, text)) return

    select case (lowercase(trim(text)))
    case ("1", "true", ".true.", "yes", "on")
      value = .true.
      return
    case ("0", "false", ".false.", "no", "off")
      value = .false.
      return
    end select

    ! Historically some flags were written as plain integers, where any nonzero
    ! value meant "on".
    read (text, *, iostat=io) parsed
    if (io == 0) then
      value = (parsed /= 0)
      return
    end if

    print *, "Warning: ignoring invalid value for "//trim(name)//": ", trim(text)
  end subroutine env_flag

  ! .true. when <name> parses as an integer, which is then stored in value.
  logical function env_int(name, value)
    character(len=*), intent(in) :: name
    integer(C_INT), intent(out) :: value
    character(len=ENV_VALUE_LEN) :: text
    integer :: io

    env_int = .false.
    value = 0_C_INT
    if (.not. env_text(name, text)) return

    read (text, *, iostat=io) value
    env_int = (io == 0)
    if (.not. env_int) value = 0_C_INT
  end function env_int

  logical function env_int64(name, value)
    character(len=*), intent(in) :: name
    integer(C_INT64_T), intent(out) :: value
    character(len=ENV_VALUE_LEN) :: text
    integer :: io

    env_int64 = .false.
    value = 0_C_INT64_T
    if (.not. env_text(name, text)) return

    read (text, *, iostat=io) value
    env_int64 = (io == 0)
    if (.not. env_int64) value = 0_C_INT64_T
  end function env_int64

  ! .true. when <name> is set; values(1:count) then holds the leading integers of
  ! a list separated by any of " ,;:xX".  Parsing stops at the first entry that
  ! does not convert, so a trailing separator is harmless.
  logical function env_int_list(name, values, count)
    character(len=*), intent(in) :: name
    integer(C_INT), intent(out) :: values(:)
    integer, intent(out) :: count
    character(len=ENV_VALUE_LEN) :: text
    integer :: i, io

    values = 0_C_INT
    count = 0
    env_int_list = env_text(name, text)
    if (.not. env_int_list) return

    do i = 1, len_trim(text)
      if (index(",;:xX", text(i:i)) > 0) text(i:i) = " "
    end do

    read (text, *, iostat=io) values
    do i = 1, size(values)
      if (values(i) == 0_C_INT) exit
      count = count + 1
    end do
  end function env_int_list

  ! ASCII-only case folding, used here and by the ini parser to make keys and
  ! enum-valued settings insensitive to how they were typed.
  function lowercase(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, c

    do i = 1, len(text)
      c = iachar(text(i:i))
      if (c >= iachar("A") .and. c <= iachar("Z")) then
        out(i:i) = achar(c + iachar("a") - iachar("A"))
      else
        out(i:i) = text(i:i)
      end if
    end do
  end function lowercase

end module env_options
