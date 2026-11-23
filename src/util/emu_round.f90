#include "build_options.h"

! Single-precision emulation for numerical experiments (CHANNEL_EMU_ROUND).
!
! Rounds chosen arrays to the nearest float32 value, in place and on the
! device, at the points where single-precision storage would lose the
! information.  Everything is still computed in double; rounding the
! coefficients of a linear system perturbs it by eps*|A|, the same leading
! order as the backward error of solving it in single precision.
!
! CHANNEL_EMU_ROUND is a comma-separated list of
!   state   V after every substep (the solution stored in single)
!   vel     V velocity components only     phi  V scalar components only
!   oldrhs  the stored explicit terms (RK memory), after buildrhs
!   matv    the v-equation (biharmonic) matrix bands, after assembly
!   math    the eta and scalar (Helmholtz) matrix bands, after assembly
!   stash   rhs_stash, the transformed products buildrhs differentiates in y
module emu_round
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: emu_init, emu_round_array
  logical, public, save :: emu_vel = .false., emu_phi = .false., emu_oldrhs = .false.
  logical, public, save :: emu_matv = .false., emu_math = .false., emu_stash = .false.
  logical, public, save :: emu_any = .false.

  interface emu_round_array
    module procedure round_c4, round_r4
  end interface emu_round_array

contains

  subroutine emu_init(verbose)
    logical, intent(in) :: verbose
    character(len=256) :: s
    integer :: n, st
    call get_environment_variable("CHANNEL_EMU_ROUND", s, n, st)
    if (st /= 0 .or. n == 0) return
    s = ',' // trim(adjustl(s)) // ','
    emu_vel = has(s, 'state') .or. has(s, 'vel')
    emu_phi = has(s, 'state') .or. has(s, 'phi')
    emu_oldrhs = has(s, 'oldrhs')
    emu_matv = has(s, 'matv'); emu_math = has(s, 'math'); emu_stash = has(s, 'stash')
    emu_any = emu_vel .or. emu_phi .or. emu_oldrhs .or. emu_matv .or. emu_math .or. emu_stash
    if (verbose) print '(A,6L2)', " CHANNEL_EMU_ROUND vel phi oldrhs matv math stash:", &
      emu_vel, emu_phi, emu_oldrhs, emu_matv, emu_math, emu_stash
  contains
    logical function has(list, word)
      character(len=*), intent(in) :: list, word
      has = index(list, ',' // word // ',') > 0
    end function has
  end subroutine emu_init

  subroutine round_c4(x)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:, :, :, :)
    integer :: i, j, k, l
    !$omp target teams distribute parallel do collapse(4) default(none) shared(x) private(i, j, k, l)
    do l = 1, size(x, 4)
      do k = 1, size(x, 3)
        do j = 1, size(x, 2)
          do i = 1, size(x, 1)
            x(i, j, k, l) = cmplx(real(real(x(i, j, k, l)), C_FLOAT), real(aimag(x(i, j, k, l)), C_FLOAT), C_DOUBLE)
          end do
        end do
      end do
    end do
  end subroutine round_c4

  subroutine round_r4(x)
    real(C_DOUBLE), intent(inout) :: x(:, :, :, :)
    integer :: i, j, k, l
    !$omp target teams distribute parallel do collapse(4) default(none) shared(x) private(i, j, k, l)
    do l = 1, size(x, 4)
      do k = 1, size(x, 3)
        do j = 1, size(x, 2)
          do i = 1, size(x, 1)
            x(i, j, k, l) = real(real(x(i, j, k, l), C_FLOAT), C_DOUBLE)
          end do
        end do
      end do
    end do
  end subroutine round_r4
end module emu_round
