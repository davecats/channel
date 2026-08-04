#include "build_options.h"

module roctx

  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  implicit none

  private
  integer :: n = 0

  public :: roctxpush, roctxpop

#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
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
#if (defined(HAVE_CUDA) && defined(HAVE_NVTX)) || defined(HAVE_HIP)
    ! Only built on profiled builds: the null-terminated copy is an automatic
    ! character temporary, and these calls sit in the y-solve inner loop.
    character(kind=c_char, len=len_trim(name) + 1) :: cname
#endif
#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
    integer(c_int) :: ignored
#endif

    n = n + 1
#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
    cname = trim(name)//c_null_char
    ignored = nvtxRangePush(cname)
#elif defined(HAVE_HIP)
    cname = trim(name)//c_null_char
    call roctxRangePush(cname)
#endif
  end subroutine roctxPush

  subroutine roctxPop(name)
    character(len=*), intent(in) :: name
#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
    integer(c_int) :: ignored
#endif

    n = n - 1
    if (n < 0) then
      print *, "invalid pop for: ", name
      return
    end if
#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
    ignored = nvtxRangePop()
#elif defined(HAVE_HIP)
    call roctxRangePop()
#endif
  end subroutine roctxPop

end module roctx
