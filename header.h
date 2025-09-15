! ==============================================
! C preprocessor defines
! ==============================================

! Force (nxd,nzd) to be at most the product of a
! power of 2 and a single factor 3
#define useFFTfit

! half or full channel
!#define halfchannel

! Add a bodyforce 
!#define bodyforce

! Measure per timestep execution time
#define chron

! Verbose echo of parallel parameters
#define mpiverbose

#define nPhi 0
