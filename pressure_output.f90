#include "header.h"

MODULE pressure_output

  USE, intrinsic :: iso_c_binding
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  USE dnsdata, ONLY: V, der, k2, ialfa, ibeta, d140, d240, d24n, ni, alfa0, beta0, factor, &
                     ny, nz, nxd, izd
  USE ffts, ONLY: FFT, IFT, RFT, HFT, VVdz, VVdx
#else
  USE dnsdata, ONLY: V, der, k2, ialfa, ibeta, d140, d240, d24n, ni, alfa0, beta0, factor, &
                     ny, nz, nxd, izd, VVdz, VVdx
  USE ffts, ONLY: FFT, IFT, RFT, HFT
#endif
  USE mpi_transpose, ONLY: ny0, nyN, nx0, nxN, nxB, nzB, nzd, nx, ierr, ylB, zpy0, zpyB, &
                           sendbuf, recvbuf, pack_zTOx, unpack_zTOx, pack_xTOz, unpack_xTOz, alltoall, &
                           fft_transpose_is_local, repack_zTOx_local, repack_xTOz_local, &
                           transpose_xz_to_y_pencil, transpose_y_pencil_to_xz
  USE y_line_solvers, ONLY: ys_solve_ghost_system
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
    real_planes = int(2*(nxd + 1), C_INT64_T)*int(nzB, C_INT64_T)*int(ny + 3, C_INT64_T)

    n_floats = 4_C_INT64_T*spectral_planes + 4_C_INT64_T*real_planes
  end subroutine get_pressure_memory_estimate

  SUBROUTINE init_pressure_output()
    IMPLICIT NONE

    if (pressure_initialized) return

    allocate (pressure_src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    allocate (pressure_src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    allocate (pressure_real0(2*(nxd + 1), nzB, ny + 3))
    allocate (pressure_real1(2*(nxd + 1), nzB, ny + 3))
    allocate (pressure_h0(2*(nxd + 1), nzB, ny + 3))
    allocate (pressure_h1(2*(nxd + 1), nzB, ny + 3))

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
    !$omp target update from(V)
    call assemble_pressure_sources()

    if (present(p_out)) call solve_pressure_field(pressure_src0, pressure_src1, p_out)
    if (present(dpdy_out)) call solve_dpdy_field(pressure_src0, pressure_src1, dpdy_out)
  END SUBROUTINE compute_pressure_output

  SUBROUTINE assemble_pressure_sources()
    IMPLICIT NONE
    integer(C_INT) :: iy, iz, ix

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(pressure_h0, pressure_h1, ny, nzB, nxd) private(iy, iz, ix)
    do iy = 1, ny + 3
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
    real(C_DOUBLE), intent(inout) :: dst(2*(nxd + 1), nzB, ny + 3)
    real(C_DOUBLE), intent(in) :: dudx(2*(nxd + 1), nzB, ny + 3)
    real(C_DOUBLE), intent(in) :: dwdz(2*(nxd + 1), nzB, ny + 3)
    integer(C_INT) :: iy, iz, ix

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(dst, dudx, dwdz, ny, nzB, nxd, factor) private(iy, iz, ix)
    do iy = 1, ny + 3
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
    real(C_DOUBLE), intent(inout) :: dst(2*(nxd + 1), nzB, ny + 3)
    real(C_DOUBLE), intent(in) :: lhs(2*(nxd + 1), nzB, ny + 3)
    real(C_DOUBLE), intent(in) :: rhs(2*(nxd + 1), nzB, ny + 3)
    real(C_DOUBLE), intent(in) :: scale
    integer(C_INT) :: iy, iz, ix

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(dst, lhs, rhs, scale, ny, nzB, nxd) private(iy, iz, ix)
    do iy = 1, ny + 3
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
    real(C_DOUBLE), intent(out) :: rx(2*(nxd + 1), nzB, ny + 3)

    call load_pressure_field_to_zbuf(kind)
    call spectral_field_to_real_x(rx)
  END SUBROUTINE build_pressure_real_x

  SUBROUTINE load_pressure_field_to_zbuf(kind)
    IMPLICIT NONE
    integer, intent(in) :: kind
    integer(C_INT) :: ix, iz, iy, jx, izd_idx

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz, nzd, nxB, ny) private(iy, jx, izd_idx)
    do iy = 1, ny + 3
      do jx = 1, nxB
        do izd_idx = 1, nzd
          VVdz(izd_idx, jx, iy, 1) = (0.0d0, 0.0d0)
        end do
      end do
    end do

    select case (kind)
    case (iu)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (iv)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = V(iy, iz, ix, 2)
          end do
        end do
      end do
    case (iw)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = V(iy, iz, ix, 3)
          end do
        end do
      end do
    case (iux)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ialfa(ix)*V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (ivx)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ialfa(ix)*V(iy, iz, ix, 2)
          end do
        end do
      end do
    case (ivz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ibeta, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ibeta(iz)*V(iy, iz, ix, 2)
          end do
        end do
      end do
    case (iwx)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ialfa(ix)*V(iy, iz, ix, 3)
          end do
        end do
      end do
    case (iuz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ibeta, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ibeta(iz)*V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (iwz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ibeta, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ibeta(iz)*V(iy, iz, ix, 3)
          end do
        end do
      end do
    case (iuxx)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ialfa(ix)*ialfa(ix)*V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (iuxz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, ibeta, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ialfa(ix)*ibeta(iz)*V(iy, iz, ix, 1)
          end do
        end do
      end do
    case (iwxz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ialfa, ibeta, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ialfa(ix)*ibeta(iz)*V(iy, iz, ix, 3)
          end do
        end do
      end do
    case (iwzz)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(VVdz, V, ibeta, izd, ny, nz, nx0, nxN) private(ix, iz, iy, jx)
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = -1, ny + 1
            jx = ix - nx0 + 1
            VVdz(izd(iz) + 1, jx, iy + 2, 1) = ibeta(iz)*ibeta(iz)*V(iy, iz, ix, 3)
          end do
        end do
      end do
    end select
  END SUBROUTINE load_pressure_field_to_zbuf

  SUBROUTINE spectral_field_to_real_x(rx)
    IMPLICIT NONE
    real(C_DOUBLE), intent(out) :: rx(2*(nxd + 1), nzB, ny + 3)
#ifdef HAVE_MPI
    type(MPI_Request) :: request
    type(MPI_Status) :: status
#endif
    integer(C_INT) :: ix, iz, iy, jx

    call IFT(VVdz(:, :, :, 1), ny)
    if (fft_transpose_is_local) then
      call repack_zTOx_local(VVdz(:, :, :, 1), VVdx(:, :, :, 1), ny)
    else
      call pack_zTOx(VVdz(:, :, :, 1), sendbuf(:, 1), ny)
      call alltoall(sendbuf(:, 1), recvbuf(:, 1), request)
    end if
#ifdef HAVE_MPI
    if (.not. fft_transpose_is_local) call MPI_Wait(request, status, ierr)
#endif
    if (.not. fft_transpose_is_local) call unpack_zTOx(recvbuf(:, 1), VVdx(:, :, :, 1), ny)
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdx, nx, nxd, nzB, ny) private(ix, iz, iy)
    do iy = 1, ny + 3
      do iz = 1, nzB
        do ix = nx + 2, nxd + 1
          VVdx(ix, iz, iy, 1) = (0.0d0, 0.0d0)
        end do
      end do
    end do
    call RFT(VVdx(:, :, :, 1), rx, ny)
  END SUBROUTINE spectral_field_to_real_x

  SUBROUTINE real_x_to_spectral_field(rx, field)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: rx(2*(nxd + 1), nzB, ny + 3)
    complex(C_DOUBLE_COMPLEX), intent(out) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
#ifdef HAVE_MPI
    type(MPI_Request) :: request
    type(MPI_Status) :: status
#endif
    integer(C_INT) :: ix, iz, iy

    call HFT(rx, VVdx(:, :, :, 1), ny)
    if (fft_transpose_is_local) then
      call repack_xTOz_local(VVdx(:, :, :, 1), VVdz(:, :, :, 1), ny)
    else
      call pack_xTOz(VVdx(:, :, :, 1), sendbuf(:, 1), ny)
      call alltoall(sendbuf(:, 1), recvbuf(:, 1), request)
    end if
#ifdef HAVE_MPI
    if (.not. fft_transpose_is_local) call MPI_Wait(request, status, ierr)
#endif
    if (.not. fft_transpose_is_local) call unpack_xTOz(recvbuf(:, 1), VVdz(:, :, :, 1), ny)
    call FFT(VVdz(:, :, :, 1), ny)

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz, nx0, nxN, ny, nz, field) private(ix, iz, iy)
    do ix = nx0, nxN
      do iy = -1, ny + 1
        do iz = 0, nz
          field(iy, iz, ix) = VVdz(iz + 1, ix - nx0 + 1, iy + 2, 1)
        end do
      end do
    end do
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz,nx0, nxN, ny, nz, field, izd) private(ix, iz, iy)
    do ix = nx0, nxN
      do iy = -1, ny + 1
        do iz = -nz, -1
          field(iy, iz, ix) = VVdz(izd(iz) + 1, ix - nx0 + 1, iy + 2, 1)
        end do
      end do
    end do
  END SUBROUTINE real_x_to_spectral_field

  SUBROUTINE solve_pressure_field(src0, src1, p)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: p(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), allocatable :: src0_xz(:, :, :), src1_xz(:, :, :), p_xz(:, :, :)
    complex(C_DOUBLE_COMPLEX), allocatable :: src0_y(:, :, :), src1_y(:, :, :), p_y(:, :, :)
    real(C_DOUBLE) :: pmat(1:ny + 1, -2:2), eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    complex(C_DOUBLE_COMPLEX) :: tmp, tmp2
    integer(C_INT) :: ix, iz, iy, ix_local, iz_local, ix_global, iz_global

    allocate (src0_xz(ylB, 2*nz + 1, nxB), src1_xz(ylB, 2*nz + 1, nxB), p_xz(ylB, 2*nz + 1, nxB))
    allocate (src0_y(ny + 3, zpyB, nxB), src1_y(ny + 3, zpyB, nxB), p_y(ny + 3, zpyB, nxB))

    do ix = nx0, nxN
      do iz = -nz, nz
        src0_xz(:, iz + nz + 1, ix - nx0 + 1) = src0(:, iz, ix)
        src1_xz(:, iz + nz + 1, ix - nx0 + 1) = src1(:, iz, ix)
      end do
    end do

    call transpose_xz_to_y_pencil(src0_xz, src0_y)
    call transpose_xz_to_y_pencil(src1_xz, src1_y)

    do ix_local = 1, nxB
      ix_global = nx0 + ix_local - 1
      do iz_local = 1, zpyB
        iz_global = zpy0 + iz_local - 1 - (nz + 1)
        pmat = 0.0d0

        do iy = 1, ny - 1
          p_y(iy + 2, iz_local, ix_local) = sum(der(iy, 0, -2:2)*src0_y(iy:iy + 4, iz_local, ix_local)) + &
                                            sum(der(iy, 1, -2:2)*src1_y(iy:iy + 4, iz_local, ix_local))
        end do

        if (ix_global == 0 .and. iz_global == 0) then
          do iy = 1, ny - 1
            pmat(iy, -2:2) = der(iy, 2, -2:2)
          end do
          eqm1(:) = d140
          p_y(1, iz_local, ix_local) = ni*sum(d240(-2:2)*V(-1:3, iz_global, ix_global, 2))
          eq0(:) = d240
          p_y(2, iz_local, ix_local) = src0_y(2, iz_local, ix_local) + src1_y(2, iz_local, ix_local)
          p_y(ny + 2, iz_local, ix_local) = 0.0d0
          eqn(:) = 0.0d0
          eqn(1) = 1.0d0
          eqnp1(:) = der(ny - 1, 3, :)
          p_y(ny + 3, iz_local, ix_local) = 0.0d0
        else
          do iy = 1, ny - 1
            pmat(iy, -2:2) = der(iy, 2, -2:2) - k2(iz_global, ix_global)*der(iy, 0, -2:2)
          end do
          eqm1(:) = der(1, 3, :)
          p_y(1, iz_local, ix_local) = 0.0d0
          tmp = sum(d240(-2:2)*V(-1:3, iz_global, ix_global, 1))
          tmp2 = sum(d240(-2:2)*V(-1:3, iz_global, ix_global, 3))
          p_y(2, iz_local, ix_local) = -ni*(ialfa(ix_global)*tmp + ibeta(iz_global)*tmp2)/k2(iz_global, ix_global)
          eq0(:) = 0.0d0
          eq0(-1) = 1.0d0
          tmp = sum(d24n(-2:2)*V(ny - 3:ny + 1, iz_global, ix_global, 1))
          tmp2 = sum(d24n(-2:2)*V(ny - 3:ny + 1, iz_global, ix_global, 3))
          p_y(ny + 2, iz_local, ix_local) = -ni*(ialfa(ix_global)*tmp + ibeta(iz_global)*tmp2)/k2(iz_global, ix_global)
          eqn(:) = 0.0d0
          eqn(1) = 1.0d0
          p_y(ny + 3, iz_local, ix_local) = 0.0d0
          eqnp1(:) = der(ny - 1, 3, :)
        end if

        call ys_solve_ghost_system(p_y(:, iz_local, ix_local), pmat, eqm1, eq0, eqn, eqnp1, ny)
      end do
    end do

    call transpose_y_pencil_to_xz(p_y, p_xz)

    do ix = nx0, nxN
      do iz = -nz, nz
        p(:, iz, ix) = p_xz(:, iz + nz + 1, ix - nx0 + 1)
      end do
    end do

    deallocate (src0_xz, src1_xz, p_xz, src0_y, src1_y, p_y)
  END SUBROUTINE solve_pressure_field

  SUBROUTINE solve_dpdy_field(src0, src1, dpdy)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dpdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), allocatable :: src0_xz(:, :, :), src1_xz(:, :, :), dpdy_xz(:, :, :)
    complex(C_DOUBLE_COMPLEX), allocatable :: src0_y(:, :, :), src1_y(:, :, :), dpdy_y(:, :, :)
    real(C_DOUBLE) :: pmat(1:ny + 1, -2:2), eqm1(-2:2), eq0(-2:2), eqnp1(-2:2), eqn(-2:2)
    integer(C_INT) :: ix, iz, iy, ix_local, iz_local, ix_global, iz_global

    allocate (src0_xz(ylB, 2*nz + 1, nxB), src1_xz(ylB, 2*nz + 1, nxB), dpdy_xz(ylB, 2*nz + 1, nxB))
    allocate (src0_y(ny + 3, zpyB, nxB), src1_y(ny + 3, zpyB, nxB), dpdy_y(ny + 3, zpyB, nxB))

    do ix = nx0, nxN
      do iz = -nz, nz
        src0_xz(:, iz + nz + 1, ix - nx0 + 1) = src0(:, iz, ix)
        src1_xz(:, iz + nz + 1, ix - nx0 + 1) = src1(:, iz, ix)
      end do
    end do

    call transpose_xz_to_y_pencil(src0_xz, src0_y)
    call transpose_xz_to_y_pencil(src1_xz, src1_y)

    do ix_local = 1, nxB
      ix_global = nx0 + ix_local - 1
      do iz_local = 1, zpyB
        iz_global = zpy0 + iz_local - 1 - (nz + 1)
        pmat = 0.0d0

        do iy = 1, ny - 1
          pmat(iy, -2:2) = der(iy, 2, -2:2) - k2(iz_global, ix_global)*der(iy, 0, -2:2)
          dpdy_y(iy + 2, iz_local, ix_local) = sum(der(iy, 1, -2:2)*src0_y(iy:iy + 4, iz_local, ix_local)) + &
                                               sum(der(iy, 2, -2:2)*src1_y(iy:iy + 4, iz_local, ix_local))
        end do

        dpdy_y(2, iz_local, ix_local) = ni*sum(d240(-2:2)*V(-1:3, iz_global, ix_global, 2))
        eq0 = 0.0d0
        eq0(-1) = 1.0d0
        dpdy_y(1, iz_local, ix_local) = 0.0d0
        eqm1 = der(1, 3, :)
        dpdy_y(ny + 2, iz_local, ix_local) = ni*sum(d24n(-2:2)*V(ny - 3:ny + 1, iz_global, ix_global, 2))
        eqn = 0.0d0
        eqn(1) = 1.0d0
        dpdy_y(ny + 3, iz_local, ix_local) = 0.0d0
        eqnp1 = der(ny - 1, 3, :)

        call ys_solve_ghost_system(dpdy_y(:, iz_local, ix_local), pmat, eqm1, eq0, eqn, eqnp1, ny)
      end do
    end do

    call transpose_y_pencil_to_xz(dpdy_y, dpdy_xz)

    do ix = nx0, nxN
      do iz = -nz, nz
        dpdy(:, iz, ix) = dpdy_xz(:, iz + nz + 1, ix - nx0 + 1)
      end do
    end do

    deallocate (src0_xz, src1_xz, dpdy_xz, src0_y, src1_y, dpdy_y)
  END SUBROUTINE solve_dpdy_field

END MODULE pressure_output
