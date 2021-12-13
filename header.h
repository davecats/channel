! ==============================================
! C preprocessor defines
! ==============================================
! Stop program when warnings or interactive prompts
! are met; suggested for execution on cluster.
#define warnings_are_fatal !      <-----------------------------------------------------------------------------------------
! Use nonblocking communication in XZ directions
#define nonblockingXZ
! Use nonblocking communication in Y direction
#define nonblockingY
! Force (nxd,nzd) to be at most the product of a
! power of 2 and a single factor 3
#define useFFTfit
! half or full channel
!#define halfchannel
! Add a bodyforce 
!#define bodyforce
! define a bodyforce in space (ibm)
!#define ibm
! Measure per timestep execution time
#define chron
! Distrubute processes among processor so that
! the y-dierction is contiguous 
#define ycontiguous
! Compute convection velocity
!#define convvel
! Verbose echo of parallel parameters
#define mpiverbose
! Disable code optimisation where possible (FFTW),
! in order to get code whose behaviour can be replicated
! exactly with different parallelisations. Useful for 
! testing.
!#define no_optimising_code


! ==============================================
! Define dependencies
! ==============================================
#if defined (ibm) && ! defined (bodyforce)
#define bodyforce
#endif
