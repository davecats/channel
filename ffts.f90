!============================================!
!                                            !
!        Fast Fourier Transforms (cuFFT)     !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!         (OpenMP GPU-accelerated)           !
!                                            !
!============================================!
!
! Author: Dr. Davide Gatti
! GPU Adaptation: Gemini
!
! This module provides an interface to the
! NVIDIA cuFFT library for GPU-based FFTs.
!

#include "header.h"

MODULE ffts

  USE, intrinsic :: iso_c_binding
  USE cudafor

  IMPLICIT NONE

  ! cuFFT plan handles
  INTEGER(C_INT) :: plan_pFFT, plan_pIFT, plan_pRFT, plan_pHFT


CONTAINS

   !$OMP DECLARE TARGET
   SUBROUTINE init_fft(nxd_in,nxB_in,nzd_in,nzB_in,odd_n_real,s)
     integer(C_INT), intent(in) :: nxd_in,nxB_in,nzd_in,nzB_in
     logical, optional, intent(in) :: odd_n_real
     integer, dimension(2), optional :: s
     integer(C_INT) :: n(1), inembed(1), onembed(1), idist, odist, istride, ostride, batch, rank
     INTEGER(C_INT) :: istat

     ! Z-direction transforms (C2C)
     rank = 1
     n(1) = nzd_in
     inembed = n
     onembed = n
     istride = 1
     ostride = 1
     idist = nzd_in
     odist = nzd_in
     batch = nxB_in
     istat = cufftPlanMany(plan_pFFT, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch)
     istat = cufftPlanMany(plan_pIFT, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2C, batch)

     ! X-direction transforms (C2R and R2C)
     rank = 1
     n(1) = 2*nxd_in
     inembed = [nxd_in+1]
     onembed = [2*(nxd_in+1)]
     istride = 1
     ostride = 1
     idist = nxd_in+1
     odist = 2*(nxd_in+1)
     batch = nzB_in
     istat = cufftPlanMany(plan_pRFT, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2R, batch)
     istat = cufftPlanMany(plan_pHFT, rank, n, onembed, ostride, odist, inembed, istride, idist, CUFFT_R2C, batch)

   END SUBROUTINE init_fft

   LOGICAL FUNCTION fftFIT(i) result(isFIT)
     integer(C_INT), intent(in) :: i
     ! cuFFT is efficient for sizes that are powers of 2, 3, 5, 7.
     ! This original function is still a reasonable heuristic.
     integer(C_INT) :: j
     j=i
     DO WHILE ( MOD(j,2)==0 ); j=j/2; END DO
     DO WHILE ( MOD(j,3)==0 ); j=j/3; END DO
     DO WHILE ( MOD(j,5)==0 ); j=j/5; END DO
     DO WHILE ( MOD(j,7)==0 ); j=j/7; END DO
     isFIT = (j == 1)
   END FUNCTION fftFIT

   !$OMP DECLARE TARGET
   SUBROUTINE FFT(x)
     complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:,:)
     INTEGER(C_INT) :: istat
     istat = cufftExecC2C(plan_pFFT, loc(x), loc(x), CUFFT_FORWARD)
   END SUBROUTINE FFT

   !$OMP DECLARE TARGET
   SUBROUTINE IFT(x)
     complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:,:)
     INTEGER(C_INT) :: istat
     istat = cufftExecC2C(plan_pIFT, loc(x), loc(x), CUFFT_INVERSE)
   END SUBROUTINE IFT

   !$OMP DECLARE TARGET
   SUBROUTINE RFT(x,rx)
     complex(C_DOUBLE_COMPLEX), intent(in) :: x(:,:)
     INTEGER(C_INT) :: istat
     real(C_DOUBLE), intent(out) :: rx(:,:)
     istat = cufftExecC2R(plan_pRFT, loc(x), loc(rx))
   END SUBROUTINE RFT

   !$OMP DECLARE TARGET
   SUBROUTINE HFT(rx,x)
     complex(C_DOUBLE_COMPLEX), intent(out) :: x(:,:)
     INTEGER(C_INT) :: istat
     real(C_DOUBLE), intent(in) :: rx(:,:)
     istat = cufftExecR2C(plan_pHFT, loc(rx), loc(x))
   END SUBROUTINE HFT

   !$OMP DECLARE TARGET
   SUBROUTINE free_fft()
     INTEGER(C_INT) :: istat
     istat = cufftDestroy(plan_pFFT)
     istat = cufftDestroy(plan_pIFT)
     istat = cufftDestroy(plan_pRFT)
     istat = cufftDestroy(plan_pHFT)
   END SUBROUTINE free_fft

END MODULE ffts


