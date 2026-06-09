#include "header.h"

MODULE pressure_output

  USE, intrinsic :: iso_c_binding
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  USE dnsdata, ONLY: V, der, k2, ialfa, ibeta, d140, d240, d24n, ni, alfa0, beta0, factor, &
                     ny, nz, nxd, izd, eliminate_assembled_boundaries, &
                     pack_assembled_pentadiagonal, unpack_assembled_pentadiagonal, &
                     reconstruct_assembled_boundaries
  USE ffts, ONLY: FFT, IFT, RFT, HFT, VVdz, VVdx
#else
  USE dnsdata, ONLY: V, der, k2, ialfa, ibeta, d140, d240, d24n, ni, alfa0, beta0, factor, &
                     ny, nz, nxd, izd, VVdz, VVdx, eliminate_assembled_boundaries, &
                     pack_assembled_pentadiagonal, unpack_assembled_pentadiagonal, &
                     reconstruct_assembled_boundaries
  USE ffts, ONLY: FFT, IFT, RFT, HFT
#endif
  USE mpi_transpose, ONLY: ny0, nyN, nx0, nxN, nxB, nzB, nzd, nx, npy_grid, ierr, &
                           sendbuf, recvbuf, pack_zTOx, unpack_zTOx, pack_xTOz, unpack_xTOz, alltoall, &
                           fft_transpose_is_local, repack_zTOx_local, repack_xTOz_local
  USE roctx, ONLY: roctxPush, roctxPop
  USE y_line_solvers, ONLY: ys_prepare_assembled_workspace, ys_solve_endpoint_schur, ys_solve_packed_pentadiagonal, ys_local_rhs, ys_local_operator, &
                            ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
                            ys_eqm1, ys_eq0, ys_eqn, ys_eqnp1
#ifdef HAVE_MPI
  USE mpi_f08
#endif

  IMPLICIT NONE

  private

  integer, parameter :: iu = 1, iv = 2, iw = 3
  integer, parameter :: iux = 4, ivx = 5, ivz = 6, iwx = 7, iuz = 8, iwz = 9
  integer, parameter :: iuxx = 10, iuxz = 11, iwxz = 12, iwzz = 13

  logical, save :: pressure_initialized = .false.

  ! Pressure output is recomputed on demand, so these buffers are kept alive between output calls.
  complex(C_DOUBLE_COMPLEX), allocatable, save :: pressure_src0(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: pressure_src1(:, :, :)
  real(C_DOUBLE), allocatable, save :: pressure_real0(:, :, :)
  real(C_DOUBLE), allocatable, save :: pressure_real1(:, :, :)
  real(C_DOUBLE), allocatable, save :: pressure_h0(:, :, :)
  real(C_DOUBLE), allocatable, save :: pressure_h1(:, :, :)

  public :: init_pressure_output, free_pressure_output
  public :: compute_pressure_output, compute_poisson, compute_dpdy
  public :: get_pressure_memory_estimate

CONTAINS

  subroutine get_pressure_memory_estimate(n_floats)
    implicit none
    integer(C_INT64_T), intent(out) :: n_floats
    integer(C_INT64_T) :: local_y, spectral_planes, real_planes

    local_y = int(nyN - ny0 + 5, C_INT64_T)
    spectral_planes = local_y*int(2*nz + 1, C_INT64_T)*int(nxN - nx0 + 1, C_INT64_T)
    real_planes = int(2*(nxd + 1), C_INT64_T)*int(nzB, C_INT64_T)*int(nyN - ny0 + 5, C_INT64_T)

    n_floats = 4_C_INT64_T*spectral_planes + 4_C_INT64_T*real_planes
  end subroutine get_pressure_memory_estimate

  SUBROUTINE init_pressure_output()
    IMPLICIT NONE

    if (pressure_initialized) return

    allocate (pressure_src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    allocate (pressure_src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    allocate (pressure_real0(2*(nxd + 1), nzB, ny0 - 2:nyN + 2))
    allocate (pressure_real1(2*(nxd + 1), nzB, ny0 - 2:nyN + 2))
    allocate (pressure_h0(2*(nxd + 1), nzB, ny0 - 2:nyN + 2))
    allocate (pressure_h1(2*(nxd + 1), nzB, ny0 - 2:nyN + 2))

    pressure_src0 = (0.0d0, 0.0d0)
    pressure_src1 = (0.0d0, 0.0d0)
    pressure_real0 = 0.0d0
    pressure_real1 = 0.0d0
    pressure_h0 = 0.0d0
    pressure_h1 = 0.0d0

    !$omp target enter data map(alloc: pressure_src0, pressure_src1, pressure_real0, pressure_real1, pressure_h0, pressure_h1)

    pressure_initialized = .true.
  END SUBROUTINE init_pressure_output

  SUBROUTINE free_pressure_output()
    IMPLICIT NONE

    if (.not. pressure_initialized) return

    !$omp target exit data map(delete: pressure_src0, pressure_src1, pressure_real0, pressure_real1, pressure_h0, pressure_h1)
    deallocate (pressure_src0, pressure_src1, pressure_real0, pressure_real1, pressure_h0, pressure_h1)

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
    call assemble_pressure_sources()

    if (present(p_out)) then
      call solve_pressure_field(pressure_src0, pressure_src1, p_out)
      !$omp target update to(p_out)
    end if
    if (present(dpdy_out)) then
      call solve_dpdy_field(pressure_src0, pressure_src1, dpdy_out)
      !$omp target update to(dpdy_out)
    end if
  END SUBROUTINE compute_pressure_output

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
    integer(C_INT) :: ix, iz, iy, jx, y_first, y_last
    y_first = ny0
    y_last = nyN

    call IFT(VVdz(:, :, :, 1))
    if (fft_transpose_is_local) then
      call repack_zTOx_local(VVdz(:, :, :, 1), VVdx(:, :, :, 1), ny)
    else
      call pack_zTOx(VVdz(:, :, :, 1), sendbuf(:, 1), ny)
      call alltoall(sendbuf(:, 1), recvbuf(:, 1), request, "zTOx pressure_spectral_to_real")
    end if
#ifdef HAVE_MPI
    if (.not. fft_transpose_is_local) then
      call roctxPush("MPI_Wait zTOx pressure_spectral_to_real")
      call MPI_Wait(request, status, ierr)
      call roctxPop("MPI_Wait zTOx pressure_spectral_to_real")
    end if
#endif
    if (.not. fft_transpose_is_local) call unpack_zTOx(recvbuf(:, 1), VVdx(:, :, :, 1), ny)
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
      call repack_xTOz_local(VVdx(:, :, :, 1), VVdz(:, :, :, 1), ny)
    else
      call pack_xTOz(VVdx(:, :, :, 1), sendbuf(:, 1), ny)
      call alltoall(sendbuf(:, 1), recvbuf(:, 1), request, "xTOz pressure_real_to_spectral")
    end if
#ifdef HAVE_MPI
    if (.not. fft_transpose_is_local) then
      call roctxPush("MPI_Wait xTOz pressure_real_to_spectral")
      call MPI_Wait(request, status, ierr)
      call roctxPop("MPI_Wait xTOz pressure_real_to_spectral")
    end if
#endif
    if (.not. fft_transpose_is_local) call unpack_xTOz(recvbuf(:, 1), VVdz(:, :, :, 1), ny)
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

  SUBROUTINE solve_pressure_field(src0, src1, p)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: p(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    call solve_pressure_like_field(src0, src1, p, .false.)
  END SUBROUTINE solve_pressure_field

  SUBROUTINE solve_dpdy_field(src0, src1, dpdy)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dpdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    call solve_pressure_like_field(src0, src1, dpdy, .true.)
  END SUBROUTINE solve_dpdy_field

  subroutine solve_pressure_like_field(src0, src1, dst, solve_dpdy)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dst(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    logical, intent(in) :: solve_dpdy
    integer(C_INT) :: ix, iz, iy, iline, nlines_z, row_start, row_end, line_count
    logical :: has_lower_boundary, has_upper_boundary

    nlines_z = 2*nz + 1
    row_start = ny0
    row_end = nyN
    line_count = (nxN - nx0 + 1)*(2*nz + 1)
    call ys_prepare_assembled_workspace(ny, nz, ny0, nyN, line_count, .true.)

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(src0, src1, V, der, k2, ialfa, ibeta, ni, d140, d240, d24n, ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, &
    !$omp& ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, ys_eqm1, ys_eq0, ys_eqn, ys_eqnp1, &
    !$omp& row_start, row_end, solve_dpdy, nlines_z, nx0, nxN, nz, ny) &
    !$omp private(ix, iz, iy, iline)
    do ix = nx0, nxN
      do iz = -nz, nz
        iline = (ix - nx0)*nlines_z + iz + nz + 1

        do iy = row_start, row_end
          if (solve_dpdy) then
            ys_local_operator(iy, -2:2, iline) = der(iy, 2, -2:2) - k2(iz, ix)*der(iy, 0, -2:2)
            ys_local_rhs(iy, iline) = sum(der(iy, 1, -2:2)*src0(iy - 2:iy + 2, iz, ix)) + &
                                      sum(der(iy, 2, -2:2)*src1(iy - 2:iy + 2, iz, ix))
          else
            ys_local_rhs(iy, iline) = sum(der(iy, 0, -2:2)*src0(iy - 2:iy + 2, iz, ix)) + &
                                      sum(der(iy, 1, -2:2)*src1(iy - 2:iy + 2, iz, ix))
            if (ix == 0 .and. iz == 0) then
              ys_local_operator(iy, -2:2, iline) = der(iy, 2, -2:2)
            else
              ys_local_operator(iy, -2:2, iline) = der(iy, 2, -2:2) - k2(iz, ix)*der(iy, 0, -2:2)
            end if
          end if
        end do

        if (solve_dpdy) then
          if (row_start == 1) then
            ys_eq0(:, iline) = 0.0d0
            ys_eq0(-1, iline) = 1.0d0
            ys_eqm1(:, iline) = der(1, 3, :)
            ys_lower_ghost_rhs(iline) = (0.0d0, 0.0d0)
            ys_lower_boundary_rhs(iline) = ni*sum(d240(-2:2)*V(-1:3, iz, ix, 2))
          end if
          if (row_end == ny - 1) then
            ys_eqn(:, iline) = 0.0d0
            ys_eqn(1, iline) = 1.0d0
            ys_eqnp1(:, iline) = der(ny - 1, 3, :)
            ys_upper_ghost_rhs(iline) = (0.0d0, 0.0d0)
            ys_upper_boundary_rhs(iline) = ni*sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, 2))
          end if
        else if (ix == 0 .and. iz == 0) then
          if (row_start == 1) then
            ys_eqm1(:, iline) = d140
            ys_eq0(:, iline) = d240
            ys_lower_ghost_rhs(iline) = ni*sum(d240(-2:2)*V(-1:3, iz, ix, 2))
            ys_lower_boundary_rhs(iline) = src0(0, iz, ix) + src1(0, iz, ix)
          end if
          if (row_end == ny - 1) then
            ys_eqn(:, iline) = 0.0d0
            ys_eqn(1, iline) = 1.0d0
            ys_eqnp1(:, iline) = der(ny - 1, 3, :)
            ys_upper_boundary_rhs(iline) = (0.0d0, 0.0d0)
            ys_upper_ghost_rhs(iline) = (0.0d0, 0.0d0)
          end if
        else
          if (row_start == 1) then
            ys_eqm1(:, iline) = der(1, 3, :)
            ys_eq0(:, iline) = 0.0d0
            ys_eq0(-1, iline) = 1.0d0
            ys_lower_ghost_rhs(iline) = (0.0d0, 0.0d0)
            ys_lower_boundary_rhs(iline) = -ni*(ialfa(ix)*sum(d240(-2:2)*V(-1:3, iz, ix, 1)) + &
                                                ibeta(iz)*sum(d240(-2:2)*V(-1:3, iz, ix, 3)))/k2(iz, ix)
          end if
          if (row_end == ny - 1) then
            ys_eqn(:, iline) = 0.0d0
            ys_eqn(1, iline) = 1.0d0
            ys_eqnp1(:, iline) = der(ny - 1, 3, :)
            ys_upper_ghost_rhs(iline) = (0.0d0, 0.0d0)
            ys_upper_boundary_rhs(iline) = -ni*(ialfa(ix)*sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, 1)) + &
                                                ibeta(iz)*sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, 3)))/k2(iz, ix)
          end if
        end if
      end do
    end do
    !$omp end target teams distribute parallel do

    has_lower_boundary = (ny0 == 1)
    has_upper_boundary = (nyN == ny - 1)
    call eliminate_assembled_boundaries(ny0, nyN, line_count, has_lower_boundary, has_upper_boundary)
    if (npy_grid == 1) then
      call pack_assembled_pentadiagonal(ny0, nyN, line_count)
      call ys_solve_packed_pentadiagonal(nyN - ny0 + 1, line_count, "pressure current-layout gpsv")
      call unpack_assembled_pentadiagonal(dst(ny0 - 2, -nz, nx0), ny0, nyN, line_count)
    else
      call ys_solve_endpoint_schur(dst, .true.)
    end if
    call reconstruct_assembled_boundaries(dst(ny0 - 2, -nz, nx0), ny0, nyN, line_count, has_lower_boundary, has_upper_boundary)
  end subroutine solve_pressure_like_field

END MODULE pressure_output
