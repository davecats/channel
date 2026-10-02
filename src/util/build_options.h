! ==============================================
! Compile-time switches
! ==============================================
!
! The switches below are numerical and diagnostic choices that are the same for
! every case, so they are set here.
!
! The *physics* switch is not here any more -- it is a CMake option, because
! changing the problem should not mean editing a source file:
!
!   -DCHANNEL_HALF_CHANNEL=ON   solve a half channel (defines halfchannel)
!
! It is compile-time rather than a deck setting on purpose: it changes which
! stencil rows close the wall and how the mesh is stretched, and making it
! runtime would put a branch in kernels that run every mode, every substep.
! Do not add a `#define halfchannel` below -- it would collide with the CMake
! definition.
!
! The scalar wall condition was a second such switch, CHANNEL_PHI_NEUMANN.  It
! is `[scalars] bc` in the deck now, per scalar, since the rows are assembled
! on the host once per run and no kernel ever branches on them.

! Force (nxd,nzd) to be at most the product of a
! power of 2 and a single factor 3
#define useFFTfit

! Measure per timestep execution time
#define chron

! Verbose echo of parallel parameters
#define mpiverbose
