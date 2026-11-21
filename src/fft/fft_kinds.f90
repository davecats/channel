#include "build_options.h"

! Working precision of the xz transforms: VVdz (the z transforms, FFT/IFT),
! VVdx, rVVdx and products (the x transforms, RFT/HFT, the physical-space
! fields and their products), and the xz transpose buffers.
!
! C_DOUBLE unless CMake sets CHANNEL_SINGLE.  The kernels that fill VVdz from
! V and read the transformed products back into double arrays convert.
module fft_kinds
  use, intrinsic :: iso_c_binding, only: C_FLOAT, C_DOUBLE
  implicit none
#ifdef CHANNEL_SINGLE
  integer, parameter :: fk = C_FLOAT
#else
  integer, parameter :: fk = C_DOUBLE
#endif
end module fft_kinds
