#include "header.h"

MODULE pressure_output

  USE, intrinsic :: iso_c_binding
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  USE dnsdata, ONLY: V, der, k2, ialfa, ibeta, d140, d240, d24n, ni, alfa0, beta0, factor, &
                     ny, nz, nxd, nPhi, izd, d040, zero_bc, overlapping, solve_compact_component_current_layout
  USE ffts, ONLY: FFT, IFT, RFT, HFT, VVdz, VVdx, get_fft_workspace_bytes, get_fft_workspace_bytes_for_dims, &
                  bind_fft_workspace, unbind_fft_workspace
#else
  USE dnsdata, ONLY: V, der, k2, ialfa, ibeta, d140, d240, d24n, ni, alfa0, beta0, factor, &
                     ny, nz, nxd, nPhi, izd, VVdz, VVdx, d040, zero_bc, overlapping, solve_compact_component_current_layout
  USE ffts, ONLY: FFT, IFT, RFT, HFT
#endif
  USE byte_workspace, ONLY: workspace_request, workspace_release, workspace_slice, workspace_align_offset
  USE mpi_transpose, ONLY: ny0, nyN, nx0, nxN, nxB, nzB, nzd, nx, npy_grid, ipy, iproc, ierr, &
                           sendbuf, recvbuf, pack_zTOx, unpack_zTOx, pack_xTOz, unpack_xTOz, alltoall, &
                           fft_transpose_is_local, repack_zTOx_local, repack_xTOz_local, &
                           roctxPush, roctxPop, &
                           MPI_Request, MPI_Status, MPI_Wait
  USE y_line_solvers, ONLY: ys_gpsv_owner_matrix, ys_gpsv_owner_rhs, ys_lower_ghost_owner, ys_lower_boundary_owner, &
                            ys_upper_boundary_owner, ys_upper_ghost_owner, ys_eqm1_owner, ys_eq0_owner, &
                            ys_eqn_owner, ys_eqnp1_owner

  IMPLICIT NONE

  private

  integer, parameter :: iu = 1, iv = 2, iw = 3
  integer, parameter :: iux = 4, ivx = 5, ivz = 6, iwx = 7, iuz = 8, iwz = 9
  integer, parameter :: iuxx = 10, iuxz = 11, iwxz = 12, iwzz = 13

  logical, save :: pressure_initialized = .false.

  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: pressure_src0(:, :, :)
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: pressure_src1(:, :, :)
  real(C_DOUBLE), pointer, contiguous, save :: pressure_real0(:, :, :)
  real(C_DOUBLE), pointer, contiguous, save :: pressure_real1(:, :, :)
  real(C_DOUBLE), pointer, contiguous, save :: pressure_h0(:, :, :)
  real(C_DOUBLE), pointer, contiguous, save :: pressure_h1(:, :, :)

  public :: init_pressure_output, free_pressure_output
  public :: compute_pressure_output, compute_poisson, compute_dpdy
  public :: get_pressure_memory_estimate, get_pressure_workspace_estimate

CONTAINS

  subroutine get_pressure_memory_estimate(n_floats)
    implicit none
    integer(C_INT64_T), intent(out) :: n_floats

    n_floats = 0_C_INT64_T
  end subroutine get_pressure_memory_estimate

  subroutine get_pressure_workspace_estimate(nbytes)
    implicit none
    integer(C_SIZE_T), intent(out) :: nbytes
    integer(C_SIZE_T) :: pressure_bytes
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    integer(C_SIZE_T) :: fft_offset, fft_bytes
#endif

    call get_pressure_local_workspace_bytes(pressure_bytes)
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    call get_fft_workspace_bytes_for_dims(nxd, nxB, nzd, nzB, nPhi, overlapping, fft_bytes)
    fft_offset = workspace_align_offset(pressure_bytes)
    nbytes = workspace_align_offset(fft_offset + fft_bytes)
#else
    nbytes = pressure_bytes
#endif
  end subroutine get_pressure_workspace_estimate

  SUBROUTINE init_pressure_output()
    IMPLICIT NONE

    if (pressure_initialized) return

    nullify (pressure_src0, pressure_src1, pressure_real0, pressure_real1, pressure_h0, pressure_h1)

    pressure_initialized = .true.
  END SUBROUTINE init_pressure_output

  SUBROUTINE free_pressure_output()
    IMPLICIT NONE

    if (.not. pressure_initialized) return

    if (associated(pressure_src0)) call release_pressure_workspace("pressure_output")
    pressure_initialized = .false.
  END SUBROUTINE free_pressure_output

  SUBROUTINE compute_poisson(p)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(out) :: p(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)

    call compute_pressure_output(p_out=p)
  END SUBROUTINE compute_poisson

  SUBROUTINE compute_dpdy(dpdy)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(out) :: dpdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)

    call compute_pressure_output(dpdy_out=dpdy)
  END SUBROUTINE compute_dpdy

  SUBROUTINE compute_pressure_output(p_out, dpdy_out)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), optional, intent(out) :: p_out(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), optional, intent(out) :: dpdy_out(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)

    if (.not. pressure_initialized) call init_pressure_output()

    ! Pressure is only needed on output steps, so we rebuild its source terms from the current spectral state.
    call acquire_pressure_workspace("pressure_output")
    call assemble_pressure_sources()

    if (present(p_out)) then
      call assemble_pressure_rhs(pressure_src0, pressure_src1, p_out)
    end if
    if (present(dpdy_out)) then
      call assemble_dpdy_rhs(pressure_src0, pressure_src1, dpdy_out)
    end if
    call release_pressure_workspace("pressure_output")

    if (present(p_out)) then
      call solve_compact_component_current_layout(p_out, assemble_pressure_operator, assemble_pressure_boundaries, &
                                                  0.0d0, 0.0d0, p_out, &
                                                  solve_label="pressure current-layout gpsv", &
                                                  symmetric_operator=.true., transpose_derivative=.true.)
    end if
    if (present(dpdy_out)) then
      call solve_compact_component_current_layout(dpdy_out, assemble_dpdy_operator, assemble_dpdy_boundaries, &
                                                  0.0d0, 0.0d0, dpdy_out, &
                                                  solve_label="pressure current-layout gpsv", symmetric_operator=.true.)
    end if
  END SUBROUTINE compute_pressure_output

  subroutine get_pressure_workspace_bytes(pressure_bytes, fft_offset, total_bytes)
    implicit none
    integer(C_SIZE_T), intent(out) :: pressure_bytes, fft_offset, total_bytes
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    integer(C_SIZE_T) :: fft_bytes
#endif

    call get_pressure_local_workspace_bytes(pressure_bytes)

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    call get_fft_workspace_bytes(fft_bytes)
    fft_offset = workspace_align_offset(pressure_bytes)
    total_bytes = workspace_align_offset(fft_offset + fft_bytes)
#else
    fft_offset = 0_C_SIZE_T
    total_bytes = pressure_bytes
#endif
  end subroutine get_pressure_workspace_bytes

  subroutine get_pressure_local_workspace_bytes(nbytes)
    implicit none
    integer(C_SIZE_T), intent(out) :: nbytes
    integer(C_SIZE_T) :: offset, n_src, n_real

    n_src = int(nyN - ny0 + 5, C_SIZE_T)*int(2*nz + 1, C_SIZE_T)*int(nxN - nx0 + 1, C_SIZE_T)
    n_real = int(2*(nxd + 1), C_SIZE_T)*int(nzB, C_SIZE_T)*int(nyN - ny0 + 5, C_SIZE_T)

    offset = 0_C_SIZE_T
    offset = workspace_align_offset(offset + n_src*int(C_SIZEOF((0.0_C_DOUBLE, 0.0_C_DOUBLE)), C_SIZE_T))
    offset = workspace_align_offset(offset + n_src*int(C_SIZEOF((0.0_C_DOUBLE, 0.0_C_DOUBLE)), C_SIZE_T))
    offset = workspace_align_offset(offset + n_real*int(C_SIZEOF(0.0_C_DOUBLE), C_SIZE_T))
    offset = workspace_align_offset(offset + n_real*int(C_SIZEOF(0.0_C_DOUBLE), C_SIZE_T))
    offset = workspace_align_offset(offset + n_real*int(C_SIZEOF(0.0_C_DOUBLE), C_SIZE_T))
    nbytes = workspace_align_offset(offset + n_real*int(C_SIZEOF(0.0_C_DOUBLE), C_SIZE_T))
  end subroutine get_pressure_local_workspace_bytes

  subroutine acquire_pressure_workspace(owner)
    implicit none
    character(len=*), intent(in) :: owner
    type(C_PTR) :: base
    integer(C_SIZE_T) :: pressure_bytes, fft_offset, total_bytes

    call get_pressure_workspace_bytes(pressure_bytes, fft_offset, total_bytes)
    call workspace_request(total_bytes, owner, base)
    call bind_pressure_workspace()
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    call bind_fft_workspace(fft_offset)
#endif
  end subroutine acquire_pressure_workspace

  subroutine release_pressure_workspace(owner)
    implicit none
    character(len=*), intent(in) :: owner

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    call unbind_fft_workspace()
#endif
    call unbind_pressure_workspace()
    call workspace_release(owner)
  end subroutine release_pressure_workspace

  subroutine bind_pressure_workspace()
    implicit none
    type(C_PTR) :: ptr
    complex(C_DOUBLE_COMPLEX), pointer :: cbuf(:)
    real(C_DOUBLE), pointer :: rbuf(:)
    integer(C_SIZE_T) :: offset, n_src, n_real

    n_src = int(nyN - ny0 + 5, C_SIZE_T)*int(2*nz + 1, C_SIZE_T)*int(nxN - nx0 + 1, C_SIZE_T)
    n_real = int(2*(nxd + 1), C_SIZE_T)*int(nzB, C_SIZE_T)*int(nyN - ny0 + 5, C_SIZE_T)

    offset = 0_C_SIZE_T
    call workspace_slice(offset, ptr)
    call c_f_pointer(ptr, cbuf, [int(n_src)])
    pressure_src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN) => cbuf
    offset = workspace_align_offset(offset + n_src*int(C_SIZEOF((0.0_C_DOUBLE, 0.0_C_DOUBLE)), C_SIZE_T))

    call workspace_slice(offset, ptr)
    call c_f_pointer(ptr, cbuf, [int(n_src)])
    pressure_src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN) => cbuf
    offset = workspace_align_offset(offset + n_src*int(C_SIZEOF((0.0_C_DOUBLE, 0.0_C_DOUBLE)), C_SIZE_T))

    call workspace_slice(offset, ptr)
    call c_f_pointer(ptr, rbuf, [int(n_real)])
    pressure_real0(1:2*(nxd + 1), 1:nzB, ny0 - 2:nyN + 2) => rbuf
    offset = workspace_align_offset(offset + n_real*int(C_SIZEOF(0.0_C_DOUBLE), C_SIZE_T))

    call workspace_slice(offset, ptr)
    call c_f_pointer(ptr, rbuf, [int(n_real)])
    pressure_real1(1:2*(nxd + 1), 1:nzB, ny0 - 2:nyN + 2) => rbuf
    offset = workspace_align_offset(offset + n_real*int(C_SIZEOF(0.0_C_DOUBLE), C_SIZE_T))

    call workspace_slice(offset, ptr)
    call c_f_pointer(ptr, rbuf, [int(n_real)])
    pressure_h0(1:2*(nxd + 1), 1:nzB, ny0 - 2:nyN + 2) => rbuf
    offset = workspace_align_offset(offset + n_real*int(C_SIZEOF(0.0_C_DOUBLE), C_SIZE_T))

    call workspace_slice(offset, ptr)
    call c_f_pointer(ptr, rbuf, [int(n_real)])
    pressure_h1(1:2*(nxd + 1), 1:nzB, ny0 - 2:nyN + 2) => rbuf

    !$omp target enter data map(to: pressure_src0, pressure_src1, pressure_real0, pressure_real1, pressure_h0, pressure_h1)
    !$omp target
    pressure_src0(ny0 - 2, -nz, nx0) = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
    pressure_src1(ny0 - 2, -nz, nx0) = (0.0_C_DOUBLE, 0.0_C_DOUBLE)
    pressure_real0(1, 1, ny0 - 2) = 0.0_C_DOUBLE
    pressure_real1(1, 1, ny0 - 2) = 0.0_C_DOUBLE
    pressure_h0(1, 1, ny0 - 2) = 0.0_C_DOUBLE
    pressure_h1(1, 1, ny0 - 2) = 0.0_C_DOUBLE
    !$omp end target
  end subroutine bind_pressure_workspace

  subroutine unbind_pressure_workspace()
    implicit none

    !$omp target exit data map(release: pressure_src0, pressure_src1, pressure_real0, pressure_real1, pressure_h0, pressure_h1)
    nullify (pressure_src0, pressure_src1, pressure_real0, pressure_real1, pressure_h0, pressure_h1)
  end subroutine unbind_pressure_workspace

  SUBROUTINE assemble_pressure_sources()
    IMPLICIT NONE
    integer(C_INT) :: iy, iz, ix, y_first, y_last
    y_first = ny0
    y_last = nyN

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(pressure_h0, pressure_h1, ny, y_first, y_last, nzB, nxd) private(iy, iz, ix)
    do iy = y_first - 2, y_last + 2
      do iz = 1, nzB
        do ix = 1, 2*(nxd + 1)
          pressure_h0(ix, iz, iy) = 0.0d0
          pressure_h1(ix, iz, iy) = 0.0d0
        end do
      end do
    end do
    call build_pressure_real_x(iux, pressure_real0)
    call build_pressure_real_x(iwz, pressure_real1)
    call accumulate_h0_diagonal(pressure_h0, pressure_real0, pressure_real1)

    call build_pressure_real_x(iwx, pressure_real0)
    call build_pressure_real_x(iuz, pressure_real1)
    call accumulate_scaled_product(pressure_h0, pressure_real0, pressure_real1, -2.0d0*factor)

    call build_pressure_real_x(iu, pressure_real0)
    call build_pressure_real_x(iuxx, pressure_real1)
    call accumulate_scaled_product(pressure_h0, pressure_real0, pressure_real1, -2.0d0*factor)
    call build_pressure_real_x(iwxz, pressure_real1)
    call accumulate_scaled_product(pressure_h0, pressure_real0, pressure_real1, -2.0d0*factor)

    call build_pressure_real_x(iw, pressure_real0)
    call build_pressure_real_x(iuxz, pressure_real1)
    call accumulate_scaled_product(pressure_h0, pressure_real0, pressure_real1, -2.0d0*factor)
    call build_pressure_real_x(iwzz, pressure_real1)
    call accumulate_scaled_product(pressure_h0, pressure_real0, pressure_real1, -2.0d0*factor)

    call build_pressure_real_x(iu, pressure_real0)
    call build_pressure_real_x(ivx, pressure_real1)
    call accumulate_scaled_product(pressure_h1, pressure_real0, pressure_real1, -2.0d0*factor)

    call build_pressure_real_x(iw, pressure_real0)
    call build_pressure_real_x(ivz, pressure_real1)
    call accumulate_scaled_product(pressure_h1, pressure_real0, pressure_real1, -2.0d0*factor)

    call real_x_to_spectral_field(pressure_h0, pressure_src0)
    call real_x_to_spectral_field(pressure_h1, pressure_src1)
  END SUBROUTINE assemble_pressure_sources

  SUBROUTINE accumulate_h0_diagonal(dst, dudx, dwdz)
    IMPLICIT NONE
    real(C_DOUBLE), intent(inout) :: dst(:, :, ny0 - 2:)
    real(C_DOUBLE), intent(in) :: dudx(:, :, ny0 - 2:)
    real(C_DOUBLE), intent(in) :: dwdz(:, :, ny0 - 2:)
    integer(C_INT) :: iy, iz, ix, y_first, y_last
    y_first = ny0
    y_last = nyN

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(dst, dudx, dwdz, ny, y_first, y_last, nzB, nxd, factor) private(iy, iz, ix)
    do iy = y_first - 2, y_last + 2
      do iz = 1, nzB
        do ix = 1, 2*nxd
          dst(ix, iz, iy) = dst(ix, iz, iy) - 2.0d0*( &
                            dudx(ix, iz, iy)*dudx(ix, iz, iy) + &
                            dudx(ix, iz, iy)*dwdz(ix, iz, iy) + &
                            dwdz(ix, iz, iy)*dwdz(ix, iz, iy))*factor
        end do
      end do
    end do
  END SUBROUTINE accumulate_h0_diagonal

  SUBROUTINE accumulate_scaled_product(dst, lhs, rhs, scale)
    IMPLICIT NONE
    real(C_DOUBLE), intent(inout) :: dst(:, :, ny0 - 2:)
    real(C_DOUBLE), intent(in) :: lhs(:, :, ny0 - 2:)
    real(C_DOUBLE), intent(in) :: rhs(:, :, ny0 - 2:)
    real(C_DOUBLE), intent(in) :: scale
    integer(C_INT) :: iy, iz, ix, y_first, y_last
    y_first = ny0
    y_last = nyN

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(dst, lhs, rhs, scale, ny, y_first, y_last, nzB, nxd) private(iy, iz, ix)
    do iy = y_first - 2, y_last + 2
      do iz = 1, nzB
        do ix = 1, 2*nxd
          dst(ix, iz, iy) = dst(ix, iz, iy) + scale*lhs(ix, iz, iy)*rhs(ix, iz, iy)
        end do
      end do
    end do
  END SUBROUTINE accumulate_scaled_product

  SUBROUTINE build_pressure_real_x(kind, rx)
    IMPLICIT NONE
    integer, intent(in) :: kind
    real(C_DOUBLE), intent(out) :: rx(:, :, ny0 - 2:)

    call load_pressure_field_to_zbuf(kind)
    call spectral_field_to_real_x(rx)
  END SUBROUTINE build_pressure_real_x

  SUBROUTINE load_pressure_field_to_zbuf(kind)
    IMPLICIT NONE
    integer, intent(in) :: kind
    integer(C_INT) :: ix, iz, iy, jx, izd_idx, y_first, y_last
    y_first = ny0
    y_last = nyN

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz, nzd, nxB, y_first, y_last) private(iy, jx, izd_idx)
    do iy = y_first - 2, y_last + 2
      do jx = 1, nxB
        do izd_idx = 1, nzd
          VVdz(izd_idx, jx, iy, 1) = (0.0d0, 0.0d0)
        end do
      end do
    end do

    select case (kind)
    case (iu)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (iv)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = V(iy, iz, ix, 2)
          end do
        end do
      end do
    case (iw)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = V(iy, iz, ix, 3)
          end do
        end do
      end do
    case (iux)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ialfa(ix)*V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (ivx)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ialfa(ix)*V(iy, iz, ix, 2)
          end do
        end do
      end do
    case (ivz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ibeta, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ibeta(iz)*V(iy, iz, ix, 2)
          end do
        end do
      end do
    case (iwx)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ialfa(ix)*V(iy, iz, ix, 3)
          end do
        end do
      end do
    case (iuz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ibeta, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ibeta(iz)*V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (iwz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ibeta, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ibeta(iz)*V(iy, iz, ix, 3)
          end do
        end do
      end do
    case (iuxx)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ialfa(ix)*ialfa(ix)*V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (iuxz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, ibeta, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ialfa(ix)*ibeta(iz)*V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (iwxz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, ibeta, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ialfa(ix)*ibeta(iz)*V(iy, iz, ix, 3)
          end do
        end do
      end do
    case (iwzz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ibeta, izd, ny, y_first, y_last, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = y_first - 2, y_last + 2
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy, 1) = ibeta(iz)*ibeta(iz)*V(iy, iz, ix, 3)
          end do
        end do
      end do
    end select
  END SUBROUTINE load_pressure_field_to_zbuf

  SUBROUTINE spectral_field_to_real_x(rx)
    IMPLICIT NONE
    real(C_DOUBLE), intent(out) :: rx(:, :, ny0 - 2:)
#ifdef HAVE_MPI
    type(MPI_Request) :: request
    type(MPI_Status) :: status
#endif
    integer(C_INT) :: ix, iz, iy, y_first, y_last
    y_first = ny0
    y_last = nyN

    call IFT(VVdz(:, :, :, 1))
    if (fft_transpose_is_local) then
      call repack_zTOx_local(VVdz(:, :, :, 1), VVdx(:, :, :, 1))
    else
      call pack_zTOx(VVdz(:, :, :, 1), sendbuf(:, 1))
      call alltoall(sendbuf(:, 1), recvbuf(:, 1), request, "zTOx pressure_spectral_to_real")
    end if
#ifdef HAVE_MPI
    if (.not. fft_transpose_is_local) then
      call roctxPush("MPI_Wait zTOx pressure_spectral_to_real")
      call MPI_Wait(request, status, ierr)
      call roctxPop("MPI_Wait zTOx pressure_spectral_to_real")
    end if
#endif
    if (.not. fft_transpose_is_local) call unpack_zTOx(recvbuf(:, 1), VVdx(:, :, :, 1))
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdx, nx, nxd, nzB, y_first, y_last) private(ix, iz, iy)
    do iy = y_first - 2, y_last + 2
      do iz = 1, nzB
        do ix = nx + 2, nxd + 1
          VVdx(ix, iz, iy, 1) = (0.0d0, 0.0d0)
        end do
      end do
    end do
    call RFT(VVdx(:, :, :, 1), rx)
  END SUBROUTINE spectral_field_to_real_x

  SUBROUTINE real_x_to_spectral_field(rx, field)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: rx(:, :, ny0 - 2:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
#ifdef HAVE_MPI
    type(MPI_Request) :: request
    type(MPI_Status) :: status
#endif
    integer(C_INT) :: ix, iz, iy, y_first, y_last
    y_first = ny0
    y_last = nyN

    call HFT(rx, VVdx(:, :, :, 1))
    if (fft_transpose_is_local) then
      call repack_xTOz_local(VVdx(:, :, :, 1), VVdz(:, :, :, 1))
    else
      call pack_xTOz(VVdx(:, :, :, 1), sendbuf(:, 1))
      call alltoall(sendbuf(:, 1), recvbuf(:, 1), request, "xTOz pressure_real_to_spectral")
    end if
#ifdef HAVE_MPI
    if (.not. fft_transpose_is_local) then
      call roctxPush("MPI_Wait xTOz pressure_real_to_spectral")
      call MPI_Wait(request, status, ierr)
      call roctxPop("MPI_Wait xTOz pressure_real_to_spectral")
    end if
#endif
    if (.not. fft_transpose_is_local) call unpack_xTOz(recvbuf(:, 1), VVdz(:, :, :, 1))
    call FFT(VVdz(:, :, :, 1))
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz, nx0, nxN, ny, y_first, y_last, nz, field) private(ix, iz, iy)
    do ix = nx0, nxN
      do iy = y_first - 2, y_last + 2
        do iz = 0, nz
          field(iy, iz, ix) = VVdz(iz + 1, ix - nx0 + 1, iy, 1)
        end do
      end do
    end do
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz,nx0, nxN, ny, y_first, y_last, nz, field, izd) private(ix, iz, iy)
    do ix = nx0, nxN
      do iy = y_first - 2, y_last + 2
        do iz = -nz, -1
          field(iy, iz, ix) = VVdz(izd(iz) + 1, ix - nx0 + 1, iy, 1)
        end do
      end do
    end do
  END SUBROUTINE real_x_to_spectral_field

  subroutine pressure_velocity_view(owner_src, u_owner, v_owner, w_owner, ix0_owner, ixN_owner)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    complex(C_DOUBLE_COMPLEX), pointer :: u_owner(:, :, :), v_owner(:, :, :), w_owner(:, :, :)
    integer(C_INT), intent(out) :: ix0_owner, ixN_owner

    ix0_owner = lbound(owner_src, 3)
    ixN_owner = ubound(owner_src, 3)

    u_owner(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN) => V(:, :, :, 1)
    v_owner(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN) => V(:, :, :, 2)
    w_owner(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN) => V(:, :, :, 3)
  end subroutine pressure_velocity_view

  subroutine assemble_pressure_rhs(src0, src1, rhs)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: rhs(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer(C_INT) :: ix, iz, iy, y_first, y_last

    y_first = ny0
    y_last = nyN

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(src0, src1, rhs, der, y_first, y_last, nz, nx0, nxN) &
    !$omp private(ix, iz, iy)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = y_first, y_last
          rhs(iy, iz, ix) = sum(der(iy, 0, -2:2)*src0(iy - 2:iy + 2, iz, ix)) + &
                            sum(der(iy, 1, -2:2)*src1(iy - 2:iy + 2, iz, ix))
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    if (ny0 == 1 .and. nx0 <= 0 .and. nxN >= 0) then
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(rhs, src0, src1) private(ix)
      do ix = 0, 0
        rhs(0, 0, 0) = src0(0, 0, 0) + src1(0, 0, 0)
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine assemble_pressure_rhs

  subroutine assemble_dpdy_rhs(src0, src1, rhs)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: rhs(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer(C_INT) :: ix, iz, iy, y_first, y_last

    y_first = ny0
    y_last = nyN

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(src0, src1, rhs, der, y_first, y_last, nz, nx0, nxN) &
    !$omp private(ix, iz, iy)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = y_first, y_last
          rhs(iy, iz, ix) = sum(der(iy, 1, -2:2)*src0(iy - 2:iy + 2, iz, ix)) + &
                            sum(der(iy, 2, -2:2)*src1(iy - 2:iy + 2, iz, ix))
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine assemble_dpdy_rhs

  subroutine assemble_pressure_operator(owner_src, lambda_coeff, diffusion_coeff, row_start, row_end, derivative_order)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff
    integer(C_INT), intent(in) :: row_start, row_end
    integer(C_INT), optional, intent(in) :: derivative_order
    integer(C_INT) :: ix, iz, iy, ix0_owner, ixN_owner
    associate (unused_lambda => lambda_coeff, unused_diffusion => diffusion_coeff)
    end associate

    ix0_owner = lbound(owner_src, 3)
    ixN_owner = ubound(owner_src, 3)

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp& shared(owner_src, der, k2, row_start, row_end, nz, ix0_owner, ixN_owner) &
    !$omp& shared(ys_gpsv_owner_matrix, ys_gpsv_owner_rhs) &
    !$omp private(ix, iz, iy)
    do ix = ix0_owner, ixN_owner
      do iz = -nz, nz
        do iy = row_start, row_end
          ys_gpsv_owner_rhs(iz, ix, iy) = owner_src(iy, iz, ix)
          if (ix == 0 .and. iz == 0) then
            ys_gpsv_owner_matrix(iz, ix, iy, -2:2) = cmplx(der(iy, 2, -2:2), 0.0d0, kind=C_DOUBLE)
          else
            ys_gpsv_owner_matrix(iz, ix, iy, -2:2) = cmplx(der(iy, 2, -2:2) - k2(iz, ix)*der(iy, 0, -2:2), &
                                                           0.0d0, kind=C_DOUBLE)
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine assemble_pressure_operator

  subroutine assemble_dpdy_operator(owner_src, lambda_coeff, diffusion_coeff, row_start, row_end, derivative_order)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff
    integer(C_INT), intent(in) :: row_start, row_end
    integer(C_INT), optional, intent(in) :: derivative_order
    integer(C_INT) :: ix, iz, iy, ix0_owner, ixN_owner
    associate (unused_lambda => lambda_coeff, unused_diffusion => diffusion_coeff)
    end associate

    ix0_owner = lbound(owner_src, 3)
    ixN_owner = ubound(owner_src, 3)

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp& shared(owner_src, der, k2, row_start, row_end, nz, ix0_owner, ixN_owner) &
    !$omp& shared(ys_gpsv_owner_matrix, ys_gpsv_owner_rhs) &
    !$omp private(ix, iz, iy)
    do ix = ix0_owner, ixN_owner
      do iz = -nz, nz
        do iy = row_start, row_end
          ys_gpsv_owner_rhs(iz, ix, iy) = owner_src(iy, iz, ix)
          ys_gpsv_owner_matrix(iz, ix, iy, -2:2) = cmplx(der(iy, 2, -2:2) - k2(iz, ix)*der(iy, 0, -2:2), &
                                                         0.0d0, kind=C_DOUBLE)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine assemble_dpdy_operator

  subroutine assemble_pressure_boundaries(owner_src, row_start, row_end, has_lower_boundary, has_upper_boundary, derivative_order)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    integer(C_INT), intent(in) :: row_start, row_end
    logical, intent(in) :: has_lower_boundary, has_upper_boundary
    integer(C_INT), optional, intent(in) :: derivative_order
    integer(C_INT) :: ix, iz, ix0_owner, ixN_owner
    complex(C_DOUBLE_COMPLEX), pointer :: u_owner(:, :, :), v_owner(:, :, :), w_owner(:, :, :)
    associate (unused_row_start => row_start, unused_row_end => row_end)
    end associate

    call pressure_velocity_view(owner_src, u_owner, v_owner, w_owner, ix0_owner, ixN_owner)

    if (has_lower_boundary) then
      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp& shared(owner_src, u_owner, v_owner, w_owner, d140, d240, der) &
      !$omp& shared(ialfa, ibeta, k2, ni, nz, ix0_owner, ixN_owner) &
      !$omp& shared(ys_eqm1_owner, ys_eq0_owner, ys_lower_ghost_owner, ys_lower_boundary_owner) &
      !$omp private(ix, iz)
      do ix = ix0_owner, ixN_owner
        do iz = -nz, nz
          if (ix == 0 .and. iz == 0) then
            ys_eqm1_owner(:, iz, ix) = d140
            ys_eq0_owner(:, iz, ix) = d240
            ys_lower_ghost_owner(iz, ix) = ni*sum(d240(-2:2)*v_owner(-1:3, iz, ix))
            ys_lower_boundary_owner(iz, ix) = owner_src(0, iz, ix)
          else
            ys_eqm1_owner(:, iz, ix) = der(1, 3, :)
            ys_eq0_owner(:, iz, ix) = 0.0d0
            ys_eq0_owner(-1, iz, ix) = 1.0d0
            ys_lower_ghost_owner(iz, ix) = (0.0d0, 0.0d0)
            ys_lower_boundary_owner(iz, ix) = -ni*( &
                                              ialfa(ix)*sum(d240(-2:2)*u_owner(-1:3, iz, ix)) + &
                                              ibeta(iz)*sum(d240(-2:2)*w_owner(-1:3, iz, ix)))/k2(iz, ix)
          end if
        end do
      end do
      !$omp end target teams distribute parallel do
    end if

    if (has_upper_boundary) then
      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp& shared(u_owner, w_owner, d24n, der, ialfa, ibeta, k2, ni, ny, nz, ix0_owner, ixN_owner) &
      !$omp& shared(ys_eqn_owner, ys_eqnp1_owner, ys_upper_boundary_owner, ys_upper_ghost_owner) &
      !$omp private(ix, iz)
      do ix = ix0_owner, ixN_owner
        do iz = -nz, nz
          ys_eqn_owner(:, iz, ix) = 0.0d0
          ys_eqn_owner(1, iz, ix) = 1.0d0
          ys_eqnp1_owner(:, iz, ix) = der(ny - 1, 3, :)
          ys_upper_ghost_owner(iz, ix) = (0.0d0, 0.0d0)
          if (ix == 0 .and. iz == 0) then
            ys_upper_boundary_owner(iz, ix) = (0.0d0, 0.0d0)
          else
            ys_upper_boundary_owner(iz, ix) = -ni*( &
                                              ialfa(ix)*sum(d24n(-2:2)*u_owner(ny - 3:ny + 1, iz, ix)) + &
                                              ibeta(iz)*sum(d24n(-2:2)*w_owner(ny - 3:ny + 1, iz, ix)))/k2(iz, ix)
          end if
        end do
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine assemble_pressure_boundaries

  subroutine assemble_dpdy_boundaries(owner_src, row_start, row_end, has_lower_boundary, has_upper_boundary, derivative_order)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    integer(C_INT), intent(in) :: row_start, row_end
    logical, intent(in) :: has_lower_boundary, has_upper_boundary
    integer(C_INT), optional, intent(in) :: derivative_order
    integer(C_INT) :: ix, iz, ix0_owner, ixN_owner
    complex(C_DOUBLE_COMPLEX), pointer :: u_owner(:, :, :), v_owner(:, :, :), w_owner(:, :, :)
    associate (unused_row_start => row_start, unused_row_end => row_end)
    end associate

    call pressure_velocity_view(owner_src, u_owner, v_owner, w_owner, ix0_owner, ixN_owner)

    if (has_lower_boundary) then
      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp& shared(v_owner, d240, der, ni, nz, ix0_owner, ixN_owner) &
      !$omp& shared(ys_eqm1_owner, ys_eq0_owner, ys_lower_ghost_owner, ys_lower_boundary_owner) &
      !$omp private(ix, iz)
      do ix = ix0_owner, ixN_owner
        do iz = -nz, nz
          ys_eq0_owner(:, iz, ix) = 0.0d0
          ys_eq0_owner(-1, iz, ix) = 1.0d0
          ys_eqm1_owner(:, iz, ix) = der(1, 3, :)
          ys_lower_ghost_owner(iz, ix) = (0.0d0, 0.0d0)
          ys_lower_boundary_owner(iz, ix) = ni*sum(d240(-2:2)*v_owner(-1:3, iz, ix))
        end do
      end do
      !$omp end target teams distribute parallel do
    end if

    if (has_upper_boundary) then
      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp& shared(v_owner, d24n, der, ni, ny, nz, ix0_owner, ixN_owner) &
      !$omp& shared(ys_eqn_owner, ys_eqnp1_owner, ys_upper_boundary_owner, ys_upper_ghost_owner) &
      !$omp private(ix, iz)
      do ix = ix0_owner, ixN_owner
        do iz = -nz, nz
          ys_eqn_owner(:, iz, ix) = 0.0d0
          ys_eqn_owner(1, iz, ix) = 1.0d0
          ys_eqnp1_owner(:, iz, ix) = der(ny - 1, 3, :)
          ys_upper_ghost_owner(iz, ix) = (0.0d0, 0.0d0)
          ys_upper_boundary_owner(iz, ix) = ni*sum(d24n(-2:2)*v_owner(ny - 3:ny + 1, iz, ix))
        end do
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine assemble_dpdy_boundaries

END MODULE pressure_output
