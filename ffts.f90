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
  use mpi_transpose, only: ny0, nyN
  IMPLICIT NONE

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  complex(C_DOUBLE_COMPLEX), dimension(:, :, :, :), allocatable :: VVdz, VVdx
  real(C_DOUBLE), dimension(:, :, :, :), allocatable :: rVVdx, products
#elif defined(HAVE_FFTW)
  INCLUDE 'fftw3.f03'
  integer, save        :: plan_type = FFTW_PATIENT
  real(C_DOUBLE), dimension(:, :, :, :), pointer :: products
  TYPE(C_PTR), save    :: pFFT, pIFT, pRFT, pHFT
  complex(C_DOUBLE_COMPLEX), target, allocatable, save :: fftw_VVdz(:, :, :, :), fftw_VVdx(:, :, :, :)
  real(C_DOUBLE), target, allocatable, save :: fftw_rVVdx(:, :, :, :)
#endif
#ifdef HAVE_CUDA
  integer :: cu_pFFT, cu_pIFT, cu_pRFT, cu_pHFT
#elif HAVE_HIP
  type(c_ptr) :: hip_pFFT, hip_pIFT, hip_pRFT, hip_pHFT
#endif
  integer(C_INT), save :: fft_y0, fft_yN, fft_ny
  logical, save :: debug_fft_compare = .false.
  logical, save :: debug_fft_compare_initialized = .false.
  integer(C_INT), save :: debug_rft_prints = 0
  integer(C_INT), save :: debug_hft_prints = 0

CONTAINS

  subroutine init_debug_fft_compare_flag()
    implicit none
    character(len=16) :: env_value
    integer :: status, length

    if (debug_fft_compare_initialized) return
    debug_fft_compare_initialized = .true.
    debug_fft_compare = .false.
    call get_environment_variable("CHANNEL_DEBUG_FFT_COMPARE", env_value, length, status)
    if (status /= 0) return

    select case (adjustl(trim(env_value(:length))))
    case ("1", "true", "TRUE", "yes", "YES", "on", "ON")
      debug_fft_compare = .true.
    case ("0", "false", "FALSE", "no", "NO", "off", "OFF")
      debug_fft_compare = .false.
    case default
      print *, "Warning: invalid value for CHANNEL_DEBUG_FFT_COMPARE:", trim(env_value(:length))
    end select
  end subroutine init_debug_fft_compare_flag

  subroutine debug_compare_rft(label, x, rx)
    implicit none
    character(len=*), intent(in) :: label
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:, :, :)
    real(C_DOUBLE), intent(inout) :: rx(:, :, :)
    integer(C_INT) :: nreal, ix, iz, iy, k
    real(C_DOUBLE) :: theta, max_err, ref_val
    complex(C_DOUBLE_COMPLEX) :: phase

    call init_debug_fft_compare_flag()
    if (.not. debug_fft_compare) return
    if (debug_rft_prints >= 6) return
    if (size(x, 1) == 0 .or. size(x, 2) == 0 .or. size(x, 3) == 0) return

    nreal = 2*(size(x, 1) - 1)
    if (nreal <= 0) return

    !$omp target update from(x, rx)

    iz = 1
    iy = 1
    max_err = 0.0d0
    do ix = 1, nreal
      ref_val = real(x(1, iz, iy), C_DOUBLE) + (-1.0d0)**(ix - 1)*real(x(size(x, 1), iz, iy), C_DOUBLE)
      do k = 2, size(x, 1) - 1
        theta = 2.0d0*acos(-1.0d0)*real((ix - 1)*(k - 1), C_DOUBLE)/real(nreal, C_DOUBLE)
        phase = cmplx(cos(theta), sin(theta), C_DOUBLE_COMPLEX)
        ref_val = ref_val + 2.0d0*real(x(k, iz, iy)*phase, C_DOUBLE)
      end do
      max_err = max(max_err, abs(rx(ix, iz, iy) - ref_val))
    end do

    debug_rft_prints = debug_rft_prints + 1
    print *, "FFT_DEBUG RFT ", trim(label), " n=", nreal, " max_err=", max_err, &
      " sample_out=", rx(1, iz, iy)
  end subroutine debug_compare_rft

  subroutine debug_compare_hft(label, rx, x)
    implicit none
    character(len=*), intent(in) :: label
    real(C_DOUBLE), intent(inout) :: rx(:, :, :)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:, :, :)
    integer(C_INT) :: nreal, ix, iz, iy, k
    real(C_DOUBLE) :: theta, max_err
    complex(C_DOUBLE_COMPLEX) :: ref_val, phase

    call init_debug_fft_compare_flag()
    if (.not. debug_fft_compare) return
    if (debug_hft_prints >= 6) return
    if (size(x, 1) == 0 .or. size(x, 2) == 0 .or. size(x, 3) == 0) return

    nreal = 2*(size(x, 1) - 1)
    if (nreal <= 0) return

    !$omp target update from(rx, x)

    iz = 1
    iy = 1
    max_err = 0.0d0
    do k = 1, size(x, 1)
      ref_val = cmplx(0.0d0, 0.0d0, C_DOUBLE_COMPLEX)
      do ix = 1, nreal
        theta = -2.0d0*acos(-1.0d0)*real((ix - 1)*(k - 1), C_DOUBLE)/real(nreal, C_DOUBLE)
        phase = cmplx(cos(theta), sin(theta), C_DOUBLE_COMPLEX)
        ref_val = ref_val + rx(ix, iz, iy)*phase
      end do
      max_err = max(max_err, abs(x(k, iz, iy) - ref_val))
    end do

    debug_hft_prints = debug_hft_prints + 1
    print *, "FFT_DEBUG HFT ", trim(label), " n=", nreal, " max_err=", max_err, &
      " sample_out=", x(1, iz, iy)
  end subroutine debug_compare_hft

  subroutine get_fft_memory_estimate(nxd, nxB, nzd, nzB, nPhi, overlapping, n_floats)
    implicit none
    integer(C_INT), intent(in) :: nxd, nxB, nzd, nzB, nPhi
    logical, intent(in) :: overlapping
    integer(C_INT64_T), intent(out) :: n_floats
    integer(C_INT64_T) :: nflds
    integer(C_INT64_T) :: local_y

    nflds = int(merge(2, 1, overlapping), C_INT64_T)
    local_y = int(nyN - ny0 + 5, C_INT64_T)

    n_floats = 0_C_INT64_T
    n_floats = n_floats + 2_C_INT64_T*int(nzd, C_INT64_T)*int(nxB, C_INT64_T)*local_y*nflds
    n_floats = n_floats + 2_C_INT64_T*int(nxd + 1, C_INT64_T)*int(nzB, C_INT64_T)*local_y*nflds
    n_floats = n_floats + int(2*(nxd + 1), C_INT64_T)*int(nzB, C_INT64_T)*local_y*int(3 + nPhi, C_INT64_T)
    n_floats = n_floats + int(2*(nxd + 1), C_INT64_T)*int(nzB, C_INT64_T)*local_y*nflds
  end subroutine get_fft_memory_estimate

#ifdef HAVE_FFTW
  SUBROUTINE init_fft(VVdz, VVdx, rVVdx, nxd, nxB, nzd, nzB, nPhi, overlapping, odd_n_real, s)
    integer(C_INT), intent(in) :: nxd, nxB, nzd, nzB, nPhi
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:, :, :, :), intent(out) :: VVdx, VVdz
    real(C_DOUBLE), pointer, dimension(:, :, :, :), intent(out) :: rVVdx
    logical, intent(in) :: overlapping
    logical, optional, intent(in) :: odd_n_real
    integer, dimension(2), optional :: s
    integer, dimension(2) :: sn = 6
    integer(C_INT), dimension(1) :: n_z, n_x, rn_x
    integer :: nflds
    n_z = [nzd]; n_x = [nxd]; rn_x = [2*nxd]; 
    if (present(odd_n_real)) then
      ! notice: odd_n_real is basically .FALSE. by default
      ! meaning that by default the logical size of the real transform is even
      if (odd_n_real .eqv. .TRUE.) rn_x = rn_x - 1
    end if
    if (present(s)) sn = s

    nflds = merge(2, 1, overlapping)
    fft_y0 = ny0 - 2
    fft_yN = nyN + 2
    fft_ny = fft_yN - fft_y0 + 1

    sn(2) = fft_ny
    sn(1) = 6 + 3*nPhi
    allocate (fftw_VVdz(nzd, nxB, fft_y0:fft_yN, nflds))
    allocate (fftw_VVdx(nxd + 1, nzB, fft_y0:fft_yN, nflds))
    allocate (fftw_rVVdx(2*(nxd + 1), nzB, fft_y0:fft_yN, 3 + nPhi))
    VVdz => fftw_VVdz
    VVdx => fftw_VVdx
    rVVdx => fftw_rVVdx
    allocate (products(2*(nxd + 1), nzB, fft_y0:fft_yN, nflds))

    !$omp target enter data map(to: VVdz)
    !FFTs plans
    pFFT = fftw_plan_many_dft(1, n_z, nxB, VVdz(:, :, fft_y0, 1), n_z, 1, nzd, VVdz(:, :, fft_y0, 1), n_z, 1, nzd, FFTW_FORWARD, plan_type)
    pIFT = fftw_plan_many_dft(1, n_z, nxB, VVdz(:, :, fft_y0, 1), n_z, 1, nzd, VVdz(:, :, fft_y0, 1), n_z, 1, nzd, FFTW_BACKWARD, plan_type)
    pRFT = fftw_plan_many_dft_c2r(1, rn_x, nzB, VVdx(:, :, fft_y0, 1), n_x + 1, 1, (nxd + 1), &
                                  rVVdx(:, :, fft_y0, 1), 2*(n_x + 1), 1, 2*(nxd + 1), plan_type)
    pHFT = fftw_plan_many_dft_r2c(1, rn_x, nzB, rVVdx(:, :, fft_y0, 1), 2*(n_x + 1), 1, 2*(nxd + 1), &
                                  VVdx(:, :, fft_y0, 1), n_x + 1, 1, (nxd + 1), plan_type)
  END SUBROUTINE init_fft
#elif defined HAVE_CUDA
  SUBROUTINE init_cufft(nxd, nxB, nzd, nzB, nPhi, overlapping)
    use cufft
    IMPLICIT NONE
    integer(C_INT), intent(in) :: nxd, nxB, nzd, nzB, nPhi
    logical, intent(in) :: overlapping
    integer :: istat
    integer, dimension(1) :: n, inembed, onembed
    integer :: batch, idist, odist, istride, ostride
    integer :: nflds

    nflds = merge(2, 1, overlapping)
    fft_y0 = ny0 - 2
    fft_yN = nyN + 2
    fft_ny = fft_yN - fft_y0 + 1

    allocate (VVdz(nzd, nxB, fft_y0:fft_yN, nflds))
    allocate (VVdx(nxd + 1, nzB, fft_y0:fft_yN, nflds))
    allocate (rVVdx(2*(nxd + 1), nzB, fft_y0:fft_yN, 3 + nPhi))
    allocate (products(2*(nxd + 1), nzB, fft_y0:fft_yN, nflds))
    !$omp target enter data map(to: VVdz, VVdx, rVVdx, products)

    !FFTs plans
    istat = cufftCreate(cu_pIFT)
    istat = cufftSetAutoAllocation(cu_pIFT, 0)
    istat = cufftPlan1d(cu_pIFT, nzd, CUFFT_Z2Z, fft_ny*nxB)

    istat = cufftCreate(cu_pFFT)
    istat = cufftSetAutoAllocation(cu_pFFT, 0)
    istat = cufftPlan1d(cu_pFFT, nzd, CUFFT_Z2Z, fft_ny*nxB)

    n(1) = 2*nxd            ! length
    batch = nzB*fft_ny
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
                          inembed, istride, idist, CUFFT_D2Z, nzB*fft_ny)

  END SUBROUTINE init_cufft
#elif defined(HAVE_HIP)
  SUBROUTINE init_hipfft(nxd, nxB, nzd, nzB, nPhi, overlapping)
    use hipfort
    use hipfort_hipfft
    IMPLICIT NONE
    integer(C_INT), intent(in) :: nxd, nxB, nzd, nzB, nPhi
    logical, intent(in) :: overlapping
    integer :: istat
    integer, dimension(1), target :: n, inembed, onembed
    integer(C_INT) :: batch, idist, odist, istride, ostride
    integer :: nflds

    nflds = merge(2, 1, overlapping)
    fft_y0 = ny0 - 2
    fft_yN = nyN + 2
    fft_ny = fft_yN - fft_y0 + 1

    allocate (VVdz(nzd, nxB, fft_y0:fft_yN, nflds))
    allocate (VVdx(nxd + 1, nzB, fft_y0:fft_yN, nflds))
    allocate (rVVdx(2*(nxd + 1), nzB, fft_y0:fft_yN, 3 + nPhi))
    allocate (products(2*(nxd + 1), nzB, fft_y0:fft_yN, nflds))
    !$omp target enter data map(to: VVdz, VVdx, rVVdx, products)

    !FFTs plans
    istat = hipfftCreate(hip_pIFT)
    istat = hipfftSetAutoAllocation(hip_pIFT, 0)
    istat = hipfftPlan1d(hip_pIFT, nzd, HIPFFT_Z2Z, fft_ny*nxB)

    istat = hipfftCreate(hip_pFFT)
    istat = hipfftSetAutoAllocation(hip_pFFT, 0)
    istat = hipfftPlan1d(hip_pFFT, nzd, HIPFFT_Z2Z, fft_ny*nxB)

    n(1) = 2*nxd            ! length
    batch = nzB*fft_ny
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
                           c_loc(inembed), istride, idist, HIPFFT_D2Z, int(nzB*fft_ny, c_int))

    call init_debug_fft_compare_flag()
    if (debug_fft_compare) then
      print *, "FFT_DEBUG init_hipfft nzd=", nzd, " nxd=", nxd, " nxB=", nxB, " nzB=", nzB, &
        " fft_ny=", fft_ny, " RFT batch=", batch, " idist=", idist, " odist=", odist, &
        " inembed=", inembed(1), " onembed=", onembed(1)
    end if

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

  SUBROUTINE FFT(x)
#if defined(HAVE_HIP)
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: x(:, :, ny0 - 2:)
#else
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:, :, ny0 - 2:)
#endif
    integer :: y0
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    integer :: istat
#elif defined(HAVE_FFTW)
    integer :: i
#endif
    y0 = lbound(x, 3)
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(x)
    istat = cudaDeviceSynchronize()
    istat = cufftExecZ2Z(cu_pFFT, x(1, 1, y0), x(1, 1, y0), CUFFT_FORWARD)
    istat = cudaDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_HIP)
    !$omp target data use_device_addr(x)
    istat = hipDeviceSynchronize()
    istat = hipfftExecZ2Z(hip_pFFT, c_loc(x(1, 1, y0)), c_loc(x(1, 1, y0)), HIPFFT_FORWARD)
    istat = hipDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_FFTW)
    DO i = fft_y0, fft_yN
      CALL fftw_execute_dft(pFFT, x(:, :, i), x(:, :, i)); 
    END DO
#endif
  END SUBROUTINE FFT

  SUBROUTINE IFT(x)
#if defined(HAVE_HIP)
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: x(:, :, ny0 - 2:)
#else
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:, :, ny0 - 2:)
#endif
    integer :: y0
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    integer :: istat
#elif defined(HAVE_FFTW)
    integer :: i
#endif
    y0 = lbound(x, 3)
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(x)
    istat = cudaDeviceSynchronize()
    istat = cufftExecZ2Z(cu_pIFT, x(1, 1, y0), x(1, 1, y0), CUFFT_INVERSE)
    istat = cudaDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_HIP)
    !$omp target data use_device_addr(x)
    istat = hipDeviceSynchronize()
    istat = hipfftExecZ2Z(hip_pIFT, c_loc(x(1, 1, y0)), c_loc(x(1, 1, y0)), HIPFFT_INVERSE)
    istat = hipDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_FFTW)
    DO i = fft_y0, fft_yN
      CALL fftw_execute_dft(pIFT, x(:, :, i), x(:, :, i))
    END DO
#endif
  END SUBROUTINE IFT

  SUBROUTINE RFT(x, rx)
    IMPLICIT NONE
#if defined(HAVE_HIP)
    complex(C_DOUBLE_COMPLEX), target :: x(:, :, ny0 - 2:)
    real(C_DOUBLE), target :: rx(:, :, ny0 - 2:)
#else
    complex(C_DOUBLE_COMPLEX) :: x(:, :, ny0 - 2:)
    real(C_DOUBLE) :: rx(:, :, ny0 - 2:)
#endif
    integer :: x_y0, rx_y0, nreal
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    integer :: istat
#elif defined(HAVE_FFTW)
    integer :: i
#endif
    x_y0 = lbound(x, 3)
    rx_y0 = lbound(rx, 3)
    nreal = 2*(size(x, 1) - 1)
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(x, rx)
    istat = cudaDeviceSynchronize()
    istat = cufftExecZ2D(cu_pRFT, x(1, 1, x_y0), rx(1, 1, rx_y0))
    istat = cudaDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_HIP)
    !$omp target data use_device_addr(x, rx)
    istat = hipDeviceSynchronize()
    istat = hipfftExecZ2D(hip_pRFT, c_loc(x(1, 1, x_y0)), c_loc(rx(1, 1, rx_y0)))
    istat = hipDeviceSynchronize()
    if (nreal < size(rx, 1)) then
      rx(nreal + 1:size(rx, 1), :, :) = 0.0d0
    end if
    !$omp end target data
    call debug_compare_rft("hipfftExecZ2D", x, rx)
#elif defined(HAVE_FFTW)
    DO i = fft_y0, fft_yN
      CALL fftw_execute_dft_c2r(pRFT, x(:, :, i), rx(:, :, i))
    END DO
#endif
  END SUBROUTINE RFT

  SUBROUTINE HFT(rx, x)
    IMPLICIT NONE
#if defined(HAVE_HIP)
    complex(C_DOUBLE_COMPLEX), target :: x(:, :, ny0 - 2:)
    real(C_DOUBLE), target :: rx(:, :, ny0 - 2:)
#else
    complex(C_DOUBLE_COMPLEX) :: x(:, :, ny0 - 2:)
    real(C_DOUBLE) :: rx(:, :, ny0 - 2:)
#endif
    integer :: x_y0, rx_y0
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    integer :: istat
#elif defined(HAVE_FFTW)
    integer :: i
#endif
    x_y0 = lbound(x, 3)
    rx_y0 = lbound(rx, 3)
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(rx, x)
    istat = cudaDeviceSynchronize()
    istat = cufftExecD2Z(cu_pHFT, rx(1, 1, rx_y0), x(1, 1, x_y0))
    istat = cudaDeviceSynchronize()
    !$omp end target data
#elif defined(HAVE_HIP)
    !$omp target data use_device_addr(rx, x)
    istat = hipDeviceSynchronize()
    istat = hipfftExecD2Z(hip_pHFT, c_loc(rx(1, 1, rx_y0)), c_loc(x(1, 1, x_y0)))
    istat = hipDeviceSynchronize()
    !$omp end target data
    call debug_compare_hft("hipfftExecD2Z", rx, x)
#elif defined(HAVE_FFTW)
    DO i = fft_y0, fft_yN
      CALL fftw_execute_dft_r2c(pHFT, rx(:, :, i), x(:, :, i)); 
    END DO
#endif
  END SUBROUTINE HFT

#if defined(HAVE_FFTW)
  SUBROUTINE free_fft(VVdz, VVdx, rVVdx)
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:, :, :, :), intent(out) :: VVdx, VVdz
    real(C_DOUBLE), pointer, dimension(:, :, :, :), intent(out) :: rVVdx

    !$omp target exit data map(from: VVdz)
    if (associated(products)) deallocate (products)
    if (allocated(fftw_rVVdx)) deallocate (fftw_rVVdx)
    if (allocated(fftw_VVdx)) deallocate (fftw_VVdx)
    if (allocated(fftw_VVdz)) deallocate (fftw_VVdz)
    nullify (VVdz, VVdx, rVVdx)
  END SUBROUTINE free_fft
#endif

END MODULE ffts
