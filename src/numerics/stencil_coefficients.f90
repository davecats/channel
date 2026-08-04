#include "header.h"

! Compact (Pade-type) wall-normal derivative coefficients -- the numbers.
!
! The macros that *apply* these to a field are the separate file
! stencil_macros.fypph.  Coefficients here, code there.
!
! der(iy, n, j) holds the j = -2..2 stencil weights at node iy for
!   n = 0 : interpolation, 1 : D^1, 2 : D^2, 3 : D^4,
! obtained by requiring the stencil to be exact on polynomials.  The d* arrays
! are the corresponding one-sided rows at and just outside each wall, used to
! close the system there: d0/d1/d2 are the value/first/second derivative rows,
! and the m1 and np1 suffixes mark the ghost rows outside the lower and upper
! wall.
!
! Declarations only, and no `use` of any other module -- see channel_state for
! why that is load-bearing rather than stylistic.  setup_derivatives, which
! computes these and needs banded_lu, lives with its caller instead.
module stencil_coefficients

  use, intrinsic :: iso_c_binding

  implicit none

  real(C_DOUBLE), allocatable :: der(:, :, :)

  real(C_DOUBLE), dimension(-2:2) :: d040, d140, d240, d14m1, d24m1
  real(C_DOUBLE), dimension(-2:2) :: d04n, d14n, d24n, d14np1, d24np1
  !$omp declare target(d14np1, d14n, d14m1, d140)

end module stencil_coefficients
