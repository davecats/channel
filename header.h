! ==============================================
! C preprocessor defines
! ==============================================

! Force (nxd,nzd) to be at most the product of a
! power of 2,3,5,7
#define useFFTfit

! half or full channel
!#define halfchannel

! Add a bodyforce
!#define bodyforce

! Measure per timestep execution time
#define chron

! Verbose echo of parallel parameters
#define mpiverbose
