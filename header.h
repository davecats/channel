! ==============================================
! C preprocessor defines
! ==============================================
! Which mpi transpose do you want?
#define packunpack

! Stop program when warnings or interactive prompts
! are met; suggested for execution on cluster.
!#define warnings_are_fatal !      <-----------------------------------------------------------------------------------------

! Use nonblocking communication in XZ directions
!#define nonblockingXZ

! Use nonblocking communication in Y direction
! only if the program is NOT uiuj_*
! (uiuj does not support nonblocking comm in Y)
!#ifndef forceblockingY
!#define nonblockingY
!#endif

	! Force (nxd,nzd) to be at most the product of a
! power of 2 and a single factor 3
#define useFFTfit

! half or full channel
!#define halfchannel

! Add a bodyforce 
!#define bodyforce
!#define BODYFORCE_HEADER "body_forces/am_f1/am_pardec.inc"
!#define BODYFORCE_MODULES "body_forces/am_f1/am_f1.inc"

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

! This avoids to save Dati.cart.out every dt_field during
! the simulation. The file is saved only at the end of
! runtime
!#define runtime_avoid_savefld


! ==============================================
! Define dependencies
! ==============================================
#if defined (ibm) && ! defined (bodyforce)
#define bodyforce
#endif
