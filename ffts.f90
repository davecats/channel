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
#ifdef HAVE_HIP
  use hipfort
  use hipfort_hipfft
#endif
  USE, intrinsic :: iso_c_binding
  IMPLICIT NONE

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  complex(C_DOUBLE_COMPLEX), dimension(:, :, :, :), allocatable :: VVdz, VVdx
  real(C_DOUBLE), dimension(:, :, :, :), allocatable :: rVVdx
#elif defined(HAVE_FFTW)
  INCLUDE 'fftw3.f03'
  integer, save        :: plan_type = FFTW_PATIENT
  TYPE(C_PTR), save    :: pFFT, pIFT, pRFT, pHFT, ptrVVdx, ptrVVdz, ptrFdx, ptrFdz
#endif
#ifdef HAVE_CUDA
  integer :: cu_pFFT, cu_pIFT, cu_pRFT, cu_pHFT
#elif HAVE_HIP
  type(c_ptr) :: hip_pFFT, hip_pIFT, hip_pRFT, hip_pHFT
#endif


CONTAINS

#ifdef HAVE_FFTW
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
#elif defined HAVE_CUDA
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
    istat = cufftPlanMany(cu_pRFT, 1, n, inembed, istride, idist, &
                          onembed, ostride, odist, CUFFT_Z2D, batch)

    istat = cufftCreate(cu_pHFT)
    istat = cufftSetAutoAllocation(cu_pHFT, 0)
    istat = cufftPlanMany(cu_pHFT, 1, n, onembed, ostride, odist, &
                          inembed, istride, idist, CUFFT_D2Z, nzB*(6 + 3*nPhi)*(ny + 3))

  END SUBROUTINE init_cufft
#elif defined(HAVE_HIP)
  SUBROUTINE init_hipfft(nxd, nxB, ny, nzd, nzB, nPhi)
    use hipfort
    use hipfort_hipfft
    IMPLICIT NONE
    integer(C_INT), intent(in) :: nxd, nxB, nzd, nzB, ny, nPhi
    integer :: istat
    integer, dimension(1) :: n, inembed, onembed
    integer(C_INT) :: batch, idist, odist, istride, ostride

    allocate (VVdz(nzd, nxB, ny + 3, 6 + 3*nPhi))
    allocate (VVdx(nxd + 1, nzB, ny + 3, 6 + 3*nPhi))
    allocate (rVVdx(2*(nxd + 1), nzB, ny + 3, 6 + 3*nPhi))
    !$omp target enter data map(to: VVdz, VVdx, rVVdx)

    !FFTs plans
    istat = hipfftCreate(hip_pIFT)
    istat = hipfftSetAutoAllocation(hip_pIFT, 0)
    istat = hipfftPlan1d(hip_pIFT, nzd, HIPFFT_Z2Z, (3 + nPhi)*(ny + 3)*nxB)

    istat = hipfftCreate(hip_pFFT)
    istat = hipfftSetAutoAllocation(hip_pFFT, 0)
    istat = hipfftPlan1d(hip_pFFT, nzd, HIPFFT_Z2Z, (6 + 3*nPhi)*(ny + 3)*nxB)

    n(1) = 2*nxd            ! length
    batch = nzB*(3 + nPhi)*(ny + 3)
    istride = 1                  ! contiguous along x
    ostride = 1
    idist = nxd + 1            ! distance between consecutive complex transforms
    odist = 2*(nxd + 1)        ! distance between consecutive real outputs
    inembed(1) = nxd + 1           ! padded leading dim of complex array
    onembed(1) = 2*(nxd + 1)       ! padded leading dim of real array

    istat = hipfftCreate(hip_pRFT)
    istat = hipfftSetAutoAllocation(hip_pRFT, 0)
    istat = hipfftPlanMany(hip_pRFT, int(1, c_int), c_loc(n), c_loc(inembed), istride, idist, &
                           c_loc(onembed), ostride, odist, HIPFFT_Z2D, batch)

    istat = hipfftCreate(hip_pHFT)
    istat = hipfftSetAutoAllocation(hip_pHFT, 0)
    istat = hipfftPlanMany(hip_pHFT, int(1, c_int), c_loc(n), c_loc(onembed), ostride, odist, &
                     c_loc(inembed), istride, idist, HIPFFT_D2Z, int(nzB*(6 + 3*nPhi)*(ny + 3), c_int))

  END SUBROUTINE init_hipfft
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
#elif defined(HAVE_HIP)
    !$omp target data use_device_addr(x)
    istat = hipDeviceSynchronize()
    istat = hipfftExecZ2Z(hip_pFFT, c_loc(x(1, 1, 1, 1)), c_loc(x(1, 1, 1, 1)), HIPFFT_FORWARD)
    istat = hipDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_FFTW)
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
#elif defined(HAVE_HIP)
    !$omp target data use_device_addr(x)
    istat = hipDeviceSynchronize()
    istat = hipfftExecZ2Z(hip_pIFT, c_loc(x(1, 1, 1, 1)), c_loc(x(1, 1, 1, 1)), HIPFFT_INVERSE)
    istat = hipDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_FFTW)
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
#elif defined(HAVE_HIP)
    !$omp target data use_device_addr(x, rx)
    istat = hipDeviceSynchronize()
    istat = hipfftExecZ2D(hip_pRFT, c_loc(x(1, 1, 1, 1)), c_loc(rx(1, 1, 1, 1)))
    istat = hipDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_FFTW)
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
#elif defined(HAVE_HIP)
    !$omp target data use_device_addr(rx, x)
    istat = hipDeviceSynchronize()
    istat = hipfftExecD2Z(hip_pHFT, c_loc(rx(1, 1, 1, 1)), c_loc(x(1, 1, 1, 1)))
    istat = hipDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_FFTW)
    DO i = 1, ny + 3
      DO iV = 1, 6 + 3*nPhi
        CALL fftw_execute_dft_r2c(pHFT, rx(:, :, i, iV), x(:, :, i, iV)); 
      END DO
    END DO
#endif
  END SUBROUTINE HFT

#if defined(HAVE_FFTW)
  SUBROUTINE free_fft(VVdz, VVdx, rVVdx)
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:, :, :, :), intent(out) :: VVdx, VVdz
    real(C_DOUBLE), pointer, dimension(:, :, :, :), intent(out) :: rVVdx

    !$omp target exit data map(from: VVdz)
    CALL fftw_free(ptrVVdx); CALL fftw_free(ptrVVdz); 
  END SUBROUTINE free_fft
#endif

END MODULE ffts

