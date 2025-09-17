!============================================!
!                                            !
!           Fast Fourier Transforms          !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: Dr. Davide Gatti
! Date  : 28/Jul/2015
!

#include "header.h"

MODULE ffts
#ifdef HAVE_CUDA
  use cudafor
  use cufft
#endif
  USE, intrinsic :: iso_c_binding
  IMPLICIT NONE
  INCLUDE 'fftw3.f03'
  integer, save        :: plan_type = FFTW_PATIENT
  TYPE(C_PTR), save    :: pFFT, pIFT, pRFT, pHFT, ptrVVdx, ptrVVdz, ptrFdx, ptrFdz
#ifdef HAVE_CUDA
  complex(C_DOUBLE_COMPLEX), dimension(:, :, :, :), allocatable :: VVdz, VVdx
  real(C_DOUBLE), dimension(:, :, :, :), allocatable :: rVVdx
  integer :: cu_pFFT, cu_pIFT, cu_pRFT, cu_pHFT
#endif

CONTAINS

#ifndef HAVE_CUDA
  SUBROUTINE init_fft(VVdz, VVdx, rVVdx, nxd, nxB, ny, nzd, nzB, nPhi, odd_n_real, s)
    integer(C_INT), intent(in) :: nxd, nxB, nzd, nzB, ny, nPhi
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:, :, :, :), intent(out) :: VVdx, VVdz
    real(C_DOUBLE), pointer, dimension(:, :, :, :), intent(out) :: rVVdx
    logical, optional, intent(in) :: odd_n_real
    integer, dimension(2), optional :: s
    integer, dimension(2) :: sn = 6

    integer(C_INT), dimension(1) :: n_z, n_x, rn_x
    n_z = [nzd]; n_x = [nxd]; rn_x = [2*nxd]; 
    if (present(odd_n_real)) then
      ! notice: odd_n_real is basically .FALSE. by default
      ! meaning that by default the logical size of the real transform is even
      if (odd_n_real .eqv. .TRUE.) rn_x = rn_x - 1
    end if
    if (present(s)) sn = s

    sn(2) = ny + 3
    sn(1) = 6 + 3*nPhi
    !Allocate aligned memory
    ptrVVdz = fftw_alloc_complex(int(nxB*nzd*sn(1)*sn(2), C_SIZE_T))
    ptrVVdx = fftw_alloc_complex(int((nxd + 1)*nzB*sn(1)*sn(2), C_SIZE_T))

    !Convert C to F pointer
    CALL c_f_pointer(ptrVVdz, VVdz, [nzd, nxB, sn(2), sn(1)]); 
    CALL c_f_pointer(ptrVVdx, VVdx, [nxd + 1, nzB, sn(2), sn(1)])
    CALL c_f_pointer(ptrVVdx, rVVdx, [2*(nxd + 1), nzB, sn(2), sn(1)])

    !$omp target enter data map(to: VVdz)
    !FFTs plans
    pFFT = fftw_plan_many_dft(1, n_z, nxB, VVdz(:, :, 1, 1), n_z, 1, nzd, VVdz(:, :, 1, 1), n_z, 1, nzd, FFTW_FORWARD, plan_type)
    pIFT = fftw_plan_many_dft(1, n_z, nxB, VVdz(:, :, 1, 1), n_z, 1, nzd, VVdz(:, :, 1, 1), n_z, 1, nzd, FFTW_BACKWARD, plan_type)
    pRFT = fftw_plan_many_dft_c2r(1, rn_x, nzB, VVdx(:, :, 1, 1), n_x + 1, 1, (nxd + 1), &
                                  rVVdx(:, :, 1, 1), 2*(n_x + 1), 1, 2*(nxd + 1), plan_type)
    pHFT = fftw_plan_many_dft_r2c(1, rn_x, nzB, rVVdx(:, :, 1, 1), 2*(n_x + 1), 1, 2*(nxd + 1), &
                                  VVdx(:, :, 1, 1), n_x + 1, 1, (nxd + 1), plan_type)
  END SUBROUTINE init_fft
#endif
#ifdef HAVE_CUDA
  SUBROUTINE init_cufft(nxd, nxB, ny, nzd, nzB, nPhi)
    use cufft
    IMPLICIT NONE
    integer(C_INT), intent(in) :: nxd, nxB, nzd, nzB, ny, nPhi
    integer :: istat
    integer, dimension(1) :: n, inembed, onembed
    integer :: batch, idist, odist, istride, ostride

    allocate (VVdz(nzd, nxB, ny + 3, 6 + 3*nPhi))
    allocate (VVdx(nxd + 1, nzB, ny + 3, 6 + 3*nPhi))
    allocate (rVVdx(2*(nxd + 1), nzB, ny + 3, 6 + 3*nPhi))
    !$omp target enter data map(to: VVdz, VVdx, rVVdx)

    !!$omp target enter data map(alloc: VVdz, VVdx, rVVdx)
    !FFTs plans
    istat = cufftCreate(cu_pIFT)
    istat = cufftSetAutoAllocation(cu_pIFT, 0)
    istat = cufftPlan1d(cu_pIFT, nzd, CUFFT_Z2Z, (3 + nPhi)*(ny + 3)*nxB)

    istat = cufftCreate(cu_pFFT)
    istat = cufftSetAutoAllocation(cu_pFFT, 0)
    istat = cufftPlan1d(cu_pFFT, nzd, CUFFT_Z2Z, (6 + 3*nPhi)*(ny + 3)*nxB)

    n(1) = 2*nxd            ! length
    batch = nzB*(3 + nPhi)*(ny + 3)
    istride = 1                  ! contiguous along x
    ostride = 1
    idist = nxd + 1            ! distance between consecutive complex transforms
    odist = 2*(nxd + 1)        ! distance between consecutive real outputs
    inembed(1) = nxd + 1           ! padded leading dim of complex array
    onembed(1) = 2*(nxd + 1)       ! padded leading dim of real array

    istat = cufftCreate(cu_pRFT)
    istat = cufftSetAutoAllocation(cu_pRFT, 0)
    !istat = cufftPlan1d(cu_pRFT, 2 * nxd, CUFFT_Z2D, 1)
    istat = cufftPlanMany(cu_pRFT, 1, n, inembed, istride, idist, &
                          onembed, ostride, odist, CUFFT_Z2D, batch)

    istat = cufftCreate(cu_pHFT)
    istat = cufftSetAutoAllocation(cu_pHFT, 0)
    istat = cufftPlanMany(cu_pHFT, 1, n, onembed, ostride, odist, &
                          inembed, istride, idist, CUFFT_D2Z, nzB*(6 + 3*nPhi)*(ny + 3))

  END SUBROUTINE init_cufft
#endif

  LOGICAL FUNCTION fftFIT(i) result(isFIT)
    integer(C_INT), intent(in) :: i
    integer(C_INT) :: j
    j = i
    DO WHILE (MOD(j, 2) == 0)
      j = ishft(j, -1)
    END DO
    isFIT = ((j == 1) .OR. (j == 3))
  END FUNCTION fftFIT

  SUBROUTINE FFT(x, ny, nPhi)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:, :, :, :)
    integer(C_INT), intent(in) :: ny, nPhi
    integer :: i, iV, istat
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(x)
    istat = cudaDeviceSynchronize()
    istat = cufftExecZ2Z(cu_pFFT, x(1, 1, 1, 1), x(1, 1, 1, 1), CUFFT_FORWARD)
    istat = cudaDeviceSynchronize()
    !$omp end target data
#else
    DO i = 1, ny + 3
      DO iV = 1, 6 + 3*nPhi
        CALL fftw_execute_dft(pFFT, x(:, :, i, iV), x(:, :, i, iV)); 
      END DO
    END do
#endif
  END SUBROUTINE FFT

  SUBROUTINE IFT(x, ny, nPhi)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:, :, :, :)
    integer(C_INT), intent(in) :: ny, nPhi
    integer :: i, iV, istat
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(x)
    istat = cudaDeviceSynchronize()
    istat = cufftExecZ2Z(cu_pIFT, x(1, 1, 1, 1), x(1, 1, 1, 1), CUFFT_INVERSE)
    istat = cudaDeviceSynchronize()
    !$omp end target data
#else
    DO i = 1, ny + 3
      DO iV = 1, 3 + nPhi
        CALL fftw_execute_dft(pIFT, x(:, :, i, iV), x(:, :, i, iV))
      END DO
    END DO
#endif
  END SUBROUTINE IFT

  SUBROUTINE RFT(x, rx, ny, nPhi)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX) :: x(:, :, :, :)
    integer(C_INT), intent(in) :: ny, nPhi
    real(C_DOUBLE) :: rx(:, :, :, :)
    integer :: i, iV, istat
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(x, rx)
    istat = cudaDeviceSynchronize()
    istat = cufftExecZ2D(cu_pRFT, x(1, 1, 1, 1), rx(1, 1, 1, 1))
    istat = cudaDeviceSynchronize()
    !$omp end target data
#else
    DO i = 1, ny + 3
      DO iV = 1, 3 + nPhi
        CALL fftw_execute_dft_c2r(pRFT, x(:, :, i, iV), rx(:, :, i, iV))
      END DO
    END DO
#endif
  END SUBROUTINE RFT

  SUBROUTINE HFT(rx, x, ny, nPhi)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX) :: x(:, :, :, :)
    integer(C_INT), intent(in) :: ny, nPhi
    real(C_DOUBLE) :: rx(:, :, :, :)
    integer :: i, iV, istat
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(rx, x)
    istat = cudaDeviceSynchronize()
    istat = cufftExecD2Z(cu_pHFT, rx(1, 1, 1, 1), x(1, 1, 1, 1))
    istat = cudaDeviceSynchronize()
    !$omp end target data
#else
    DO i = 1, ny + 3
      DO iV = 1, 6 + 3*nPhi
        CALL fftw_execute_dft_r2c(pHFT, rx(:, :, i, iV), x(:, :, i, iV)); 
      END DO
    END DO
#endif
  END SUBROUTINE HFT

  SUBROUTINE free_fft(VVdz, VVdx, rVVdx)
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:, :, :, :), intent(out) :: VVdx, VVdz
    real(C_DOUBLE), pointer, dimension(:, :, :, :), intent(out) :: rVVdx

    !$omp target exit data map(from: VVdz)
    CALL fftw_free(ptrVVdx); CALL fftw_free(ptrVVdz); 
  END SUBROUTINE free_fft

END MODULE ffts

