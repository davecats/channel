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

! Storage type of the explicit terms kept between Runge-Kutta substeps (oldrhs)
! and of the transformed nonlinear products buildrhs differentiates in y
! (rhs_stash).  CMake option CHANNEL_SINGLE.  Both are computed in double
! and only stored in single: validation showed the rounding there 2-3 orders
! of magnitude below that of a single-precision state.  V and the solver
! stay double (rounding those measurably perturbs the flow).
#ifdef CHANNEL_SINGLE
#define RHS_COMPLEX complex(C_FLOAT_COMPLEX)
#else
#define RHS_COMPLEX complex(C_DOUBLE_COMPLEX)
#endif

! Coefficient type of the wall-normal band matrices (the y-solver workspace).
! Every one of them is real -- stencil weights, k^2, lambda, the wall closure
! rows -- and only the right-hand sides and solutions are complex, so real
! bands halve the solver's coefficient traffic.  cuSPARSE's gpsv takes complex
! bands, so CUDA builds keep complex ones.
#if defined(HAVE_CUDA)
#define YS_COEF complex(C_DOUBLE_COMPLEX)
#define YS_CZERO (0.0d0, 0.0d0)
#else
#define YS_COEF real(C_DOUBLE)
#define YS_CZERO 0.0d0
#endif
