! ==============================================
! Compile-time switches
! ==============================================
!
! The switches below are numerical and diagnostic choices that are the same for
! every case, so they are set here.
!
! The two *physics* switches are not here any more -- they are CMake options,
! because changing the problem should not mean editing a source file:
!
!   -DCHANNEL_HALF_CHANNEL=ON   solve a half channel (defines halfchannel)
!   -DCHANNEL_PHI_NEUMANN=ON    Neumann wall condition on the passive
!                               scalars instead of Dirichlet (defines
!                               phiNeumann)
!
! Both are compile-time rather than deck settings on purpose: they change which
! stencil rows close the wall and how the mesh is stretched, and making them
! runtime would put a branch in kernels that run every mode, every substep.
! Do not add a `#define halfchannel` or `#define phiNeumann` below -- it would
! collide with the CMake definition.

! Force (nxd,nzd) to be at most the product of a
! power of 2 and a single factor 3
#define useFFTfit

! Measure per timestep execution time
#define chron

! Verbose echo of parallel parameters
#define mpiverbose
