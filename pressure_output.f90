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
  USE mpi_transpose, ONLY: ny0, nyN, nx0, nxN, nxB, nzB, nzd, nx, ierr, &
                           sendbuf, recvbuf, pack_zTOx, unpack_zTOx, pack_xTOz, unpack_xTOz, alltoall, &
                           fft_transpose_is_local, repack_zTOx_local, repack_xTOz_local
  USE y_line_solvers, ONLY: ys_prepare_ghost_field_workspace, ys_solve_ghost_field, ys_fill_ghost_padded_field, ys_rhs_store, ys_matrix_store, &
                            ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
          do iy = ny0 - 2, nyN + 2
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
      do iy = ny0 - 2, nyN + 2
        do iz = 0, nz
          field(iy, iz, ix) = VVdz(iz + 1, ix - nx0 + 1, iy + 2, 1)
        end do
      end do
    end do
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz,nx0, nxN, ny, nz, field, izd) private(ix, iz, iy)
    do ix = nx0, nxN
      do iy = ny0 - 2, nyN + 2
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
    integer(C_INT) :: ix, iz, iy, iline, nlines_z, ix_first, ix_last, iz_first, iz_last
    integer(C_INT) :: src_y_first, src_y_last, rhs_y_first, rhs_y_last, mat_y_first, mat_y_last

    call ys_prepare_ghost_field_workspace(ny, nz, nxB)
    nlines_z = 2*nz + 1
    ix_first = nx0
    ix_last = nxN
    iz_first = -nz
    iz_last = nz
    src_y_first = lbound(src0, 1)
    src_y_last = ubound(src0, 1)
    rhs_y_first = lbound(ys_rhs_store, 1)
    rhs_y_last = ubound(ys_rhs_store, 1)
    mat_y_first = lbound(ys_matrix_store, 1)
    mat_y_last = ubound(ys_matrix_store, 1)

    !$omp target enter data map(alloc: ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store)
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(src0, src1, V, der, d140, d240, d24n, ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store, &
    !$omp& k2, ialfa, ibeta, ni, ny, ix_first, ix_last, iz_first, iz_last, nlines_z, src_y_first, src_y_last, rhs_y_first, rhs_y_last, &
    !$omp& mat_y_first, mat_y_last) private(ix, iz, iy, iline)
    do ix = ix_first, ix_last
      do iz = iz_first, iz_last
        iline = (ix - ix_first)*nlines_z + (iz - iz_first + 1)

        do iy = rhs_y_first, rhs_y_last
          ys_rhs_store(iy, iline) = (0.0d0, 0.0d0)
        end do
        do iy = mat_y_first, mat_y_last
          ys_matrix_store(iy, -2:2, iline) = 0.0d0
        end do
        ys_eqm1_store(-2:2, iline) = 0.0d0
        ys_eq0_store(-2:2, iline) = 0.0d0
        ys_eqn_store(-2:2, iline) = 0.0d0
        ys_eqnp1_store(-2:2, iline) = 0.0d0

        do iy = max(1_C_INT, src_y_first + 2), min(ny - 1, src_y_last - 2)
          ys_rhs_store(iy, iline) = sum(der(iy, 0, -2:2)*src0(iy - 2:iy + 2, iz, ix)) + &
                                    sum(der(iy, 1, -2:2)*src1(iy - 2:iy + 2, iz, ix))
        end do

        if (ix == 0 .and. iz == 0) then
          do iy = max(1_C_INT, src_y_first + 2), min(ny - 1, src_y_last - 2)
            ys_matrix_store(iy, -2:2, iline) = der(iy, 2, -2:2)
          end do
          if (src_y_first <= -1 .and. src_y_last >= 3) then
            ys_eqm1_store(:, iline) = d140
            ys_rhs_store(-1, iline) = ni*sum(d240(-2:2)*V(-1:3, iz, ix, 2))
            ys_eq0_store(:, iline) = d240
            ys_rhs_store(0, iline) = src0(0, iz, ix) + src1(0, iz, ix)
          end if
          ys_eqn_store(1, iline) = 1.0d0
          if (src_y_first <= ny - 3 .and. src_y_last >= ny + 1) then
            ys_eqnp1_store(:, iline) = der(ny - 1, 3, :)
          end if
        else
          do iy = max(1_C_INT, src_y_first + 2), min(ny - 1, src_y_last - 2)
            ys_matrix_store(iy, -2:2, iline) = der(iy, 2, -2:2) - k2(iz, ix)*der(iy, 0, -2:2)
          end do
          if (src_y_first <= -1 .and. src_y_last >= 3) then
            ys_eqm1_store(:, iline) = der(1, 3, :)
            ys_rhs_store(0, iline) = -ni*(ialfa(ix)*sum(d240(-2:2)*V(-1:3, iz, ix, 1)) + &
                                          ibeta(iz)*sum(d240(-2:2)*V(-1:3, iz, ix, 3)))/k2(iz, ix)
            ys_eq0_store(-1, iline) = 1.0d0
          end if
          if (src_y_first <= ny - 3 .and. src_y_last >= ny + 1) then
            ys_rhs_store(ny, iline) = -ni*(ialfa(ix)*sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, 1)) + &
                                           ibeta(iz)*sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, 3)))/k2(iz, ix)
            ys_eqn_store(1, iline) = 1.0d0
            ys_eqnp1_store(:, iline) = der(ny - 1, 3, :)
          end if
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    !$omp target update from(ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store)
    !$omp target exit data map(delete: ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store)

    call ys_solve_ghost_field(p(ny0:nyN, :, :), ny, nz)
    call ys_fill_ghost_padded_field(p, ny, nz)
  END SUBROUTINE solve_pressure_field

  SUBROUTINE solve_dpdy_field(src0, src1, dpdy)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dpdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer(C_INT) :: ix, iz, iy, iline, nlines_z, ix_first, ix_last, iz_first, iz_last
    integer(C_INT) :: src_y_first, src_y_last, rhs_y_first, rhs_y_last, mat_y_first, mat_y_last

    call ys_prepare_ghost_field_workspace(ny, nz, nxB)
    nlines_z = 2*nz + 1
    ix_first = nx0
    ix_last = nxN
    iz_first = -nz
    iz_last = nz
    src_y_first = lbound(src0, 1)
    src_y_last = ubound(src0, 1)
    rhs_y_first = lbound(ys_rhs_store, 1)
    rhs_y_last = ubound(ys_rhs_store, 1)
    mat_y_first = lbound(ys_matrix_store, 1)
    mat_y_last = ubound(ys_matrix_store, 1)

    !$omp target enter data map(alloc: ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store)
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(src0, src1, V, der, d240, d24n, ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store, &
    !$omp& k2, ni, ny, ix_first, ix_last, iz_first, iz_last, nlines_z, src_y_first, src_y_last, rhs_y_first, rhs_y_last, mat_y_first, mat_y_last) &
    !$omp& private(ix, iz, iy, iline)
    do ix = ix_first, ix_last
      do iz = iz_first, iz_last
        iline = (ix - ix_first)*nlines_z + (iz - iz_first + 1)

        do iy = rhs_y_first, rhs_y_last
          ys_rhs_store(iy, iline) = (0.0d0, 0.0d0)
        end do
        do iy = mat_y_first, mat_y_last
          ys_matrix_store(iy, -2:2, iline) = 0.0d0
        end do
        ys_eqm1_store(-2:2, iline) = 0.0d0
        ys_eq0_store(-2:2, iline) = 0.0d0
        ys_eqn_store(-2:2, iline) = 0.0d0
        ys_eqnp1_store(-2:2, iline) = 0.0d0

        do iy = max(1_C_INT, src_y_first + 2), min(ny - 1, src_y_last - 2)
          ys_matrix_store(iy, -2:2, iline) = der(iy, 2, -2:2) - k2(iz, ix)*der(iy, 0, -2:2)
          ys_rhs_store(iy, iline) = sum(der(iy, 1, -2:2)*src0(iy - 2:iy + 2, iz, ix)) + &
                                    sum(der(iy, 2, -2:2)*src1(iy - 2:iy + 2, iz, ix))
        end do

        if (src_y_first <= -1 .and. src_y_last >= 3) then
          ys_rhs_store(0, iline) = ni*sum(d240(-2:2)*V(-1:3, iz, ix, 2))
          ys_eq0_store(-1, iline) = 1.0d0
          ys_eqm1_store(:, iline) = der(1, 3, :)
        end if
        if (src_y_first <= ny - 3 .and. src_y_last >= ny + 1) then
          ys_rhs_store(ny, iline) = ni*sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, 2))
          ys_eqn_store(1, iline) = 1.0d0
          ys_eqnp1_store(:, iline) = der(ny - 1, 3, :)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    !$omp target update from(ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store)
    !$omp target exit data map(delete: ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store)

    call ys_solve_ghost_field(dpdy(ny0:nyN, :, :), ny, nz)
    call ys_fill_ghost_padded_field(dpdy, ny, nz)
  END SUBROUTINE solve_dpdy_field

END MODULE pressure_output
