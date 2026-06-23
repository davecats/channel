#include "header.h"

module roctx

  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  implicit none

  private
  integer :: n = 0

  public :: roctxpush, roctxpop

#if defined(HAVE_CUDA)
  interface
    integer(c_int) function nvtxrangepush(message) bind(c, name="nvtxRangePushA")
      use iso_c_binding, only: c_char, c_int
      implicit none
      character(c_char) :: message(*)
    end function nvtxrangepush

    integer(c_int) function nvtxrangepop() bind(c, name="nvtxRangePop")
      use iso_c_binding, only: c_int
      implicit none
    end function nvtxrangepop
  end interface
#elif defined(HAVE_HIP)
  interface
    subroutine roctxrangepush(message) bind(c, name="roctxRangePushA")
      use iso_c_binding, only: c_char
      implicit none
      character(c_char) :: message(*)
    end subroutine roctxrangepush

    subroutine roctxrangepop() bind(c, name="roctxRangePop")
      implicit none
    end subroutine roctxrangepop
  end interface
#endif

contains

  subroutine roctxPush(name)
    character(len=*), intent(in) :: name
    character(kind=c_char, len=len_trim(name) + 1) :: cname
#if defined(HAVE_CUDA)
    integer(c_int) :: ignored
#endif

    cname = trim(name)//c_null_char
    n = n + 1
#if defined(HAVE_CUDA)
    ignored = nvtxRangePush(cname)
#elif defined(HAVE_HIP)
    call roctxRangePush(cname)
#endif
  end subroutine roctxPush

  subroutine roctxPop(name)
    character(len=*), intent(in) :: name
#if defined(HAVE_CUDA)
    integer(c_int) :: ignored
#endif

    n = n - 1
    if (n < 0) then
      print *, "invalid pop for: ", name
      return
    end if
#if defined(HAVE_CUDA)
    ignored = nvtxRangePop()
#elif defined(HAVE_HIP)
    call roctxRangePop()
#endif
  end subroutine roctxPop

end module roctx
