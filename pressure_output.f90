#include "header.h"

MODULE pressure_output

  USE, intrinsic :: iso_c_binding
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  USE dnsdata, ONLY: V, der, k2, ialfa, ibeta, d140, d240, d24n, ni, alfa0, beta0, factor, &
                     ny, nz, nxd, izd, d040, zero_bc, use_yslab_linsolve, solve_compact_component_current_layout
  USE ffts, ONLY: FFT, IFT, RFT, HFT, VVdz, VVdx
#else
  USE dnsdata, ONLY: V, der, k2, ialfa, ibeta, d140, d240, d24n, ni, alfa0, beta0, factor, &
                     ny, nz, nxd, izd, VVdz, VVdx, d040, zero_bc, use_yslab_linsolve, solve_compact_component_current_layout
  USE ffts, ONLY: FFT, IFT, RFT, HFT
#endif
  USE mpi_transpose, ONLY: ny0, nyN, nx0, nxN, nxB, nzB, nzd, nx, npy_grid, ipy, ierr, &
                           sendbuf, recvbuf, pack_zTOx, unpack_zTOx, pack_xTOz, unpack_xTOz, alltoall, &
                           fft_transpose_is_local, repack_zTOx_local, repack_xTOz_local, &
                           yslab_workspace, prepare_yslab_scratch, yslab_copy_to_full, &
                           yslab_owned_first_line, yslab_owned_line_count
  USE roctx, ONLY: roctxPush, roctxPop
  USE y_line_solvers, ONLY: ys_gpsv_owner_matrix, ys_gpsv_owner_rhs, ys_lower_ghost_owner, ys_lower_boundary_owner, &
                            ys_upper_boundary_owner, ys_upper_ghost_owner, ys_eqm1_owner, ys_eq0_owner, &
                            ys_eqn_owner, ys_eqnp1_owner
#ifdef HAVE_MPI
  USE mpi_transpose, ONLY: yslab_transpose_to_full
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

  ! Views into spare column blocks of mpi_transpose:yslab_workspace.
  ! They are valid only between transpose_pressure_velocity_to_yslab() and the end of the current pressure solve.
  complex(C_DOUBLE_COMPLEX), pointer, save :: pressure_vslab_u(:, :, :) => null()
  complex(C_DOUBLE_COMPLEX), pointer, save :: pressure_vslab_v(:, :, :) => null()
  complex(C_DOUBLE_COMPLEX), pointer, save :: pressure_vslab_w(:, :, :) => null()

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
    nullify (pressure_vslab_u, pressure_vslab_v, pressure_vslab_w)

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

    ! For the y-slab solve, the RHS is transposed inside solve_compact_component_current_layout().
    ! The pressure BC callbacks also need wall values from V, so transpose V once here into spare yslab_workspace blocks.
    if (use_yslab_linsolve) call transpose_pressure_velocity_to_yslab()

    if (present(p_out)) then
      call solve_pressure_field(pressure_src0, pressure_src1, p_out)
    end if
    if (present(dpdy_out)) then
      call solve_dpdy_field(pressure_src0, pressure_src1, dpdy_out)
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
    integer(C_INT) :: ix, iz, iy, y_first, y_last
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

  subroutine transpose_pressure_velocity_to_yslab()
    implicit none
    integer(C_INT) :: nlines_z, total_line_count
    integer(C_INT) :: ix0_owner, ixN_owner, block_count

    nlines_z = 2*nz + 1
    total_line_count = (nxN - nx0 + 1)*nlines_z
    if (total_line_count <= 0) return

    ! Column block 1 is reserved for solve_compact_component_current_layout().
    ! Pressure stores transposed u/v/w in blocks 2:4 and lets the common solve reuse block 1.
    block_count = max(1_C_INT, yslab_owned_line_count)
    call prepare_yslab_scratch(ny + 3, 4_C_INT*block_count)

    nullify (pressure_vslab_u, pressure_vslab_v, pressure_vslab_w)

    if (yslab_owned_line_count <= 0) then
#ifdef HAVE_MPI
      ! Still participate in the collective y-slab transposes even when this rank owns no reduced-x columns.
      call yslab_transpose_to_full(V(:, :, :, 1), yslab_workspace(:, 1:1), ny, nz, total_line_count, nlines_z, .true.)
      call yslab_transpose_to_full(V(:, :, :, 2), yslab_workspace(:, 1:1), ny, nz, total_line_count, nlines_z, .true.)
      call yslab_transpose_to_full(V(:, :, :, 3), yslab_workspace(:, 1:1), ny, nz, total_line_count, nlines_z, .true.)
#endif
      return
    end if

    call roctxPush("pressure transpose V_u to yslab")
#ifdef HAVE_MPI
    call yslab_transpose_to_full(V(:, :, :, 1), yslab_workspace(:, block_count + 1:2*block_count), &
                                 ny, nz, total_line_count, nlines_z, .true.)
#else
    call yslab_copy_to_full(V(:, :, :, 1), yslab_workspace(:, block_count + 1:2*block_count), &
                            ny, nz, yslab_owned_first_line, yslab_owned_line_count, nlines_z)
#endif
    call roctxPop("pressure transpose V_u to yslab")

    call roctxPush("pressure transpose V_v to yslab")
#ifdef HAVE_MPI
    call yslab_transpose_to_full(V(:, :, :, 2), yslab_workspace(:, 2*block_count + 1:3*block_count), &
                                 ny, nz, total_line_count, nlines_z, .true.)
#else
    call yslab_copy_to_full(V(:, :, :, 2), yslab_workspace(:, 2*block_count + 1:3*block_count), &
                            ny, nz, yslab_owned_first_line, yslab_owned_line_count, nlines_z)
#endif
    call roctxPop("pressure transpose V_v to yslab")

    call roctxPush("pressure transpose V_w to yslab")
#ifdef HAVE_MPI
    call yslab_transpose_to_full(V(:, :, :, 3), yslab_workspace(:, 3*block_count + 1:4*block_count), &
                                 ny, nz, total_line_count, nlines_z, .true.)
#else
    call yslab_copy_to_full(V(:, :, :, 3), yslab_workspace(:, 3*block_count + 1:4*block_count), &
                            ny, nz, yslab_owned_first_line, yslab_owned_line_count, nlines_z)
#endif
    call roctxPop("pressure transpose V_w to yslab")

    ix0_owner = nx0 + (yslab_owned_first_line - 1)/nlines_z
    ixN_owner = ix0_owner + yslab_owned_line_count/nlines_z - 1
    pressure_vslab_u(-1:ny + 1, -nz:nz, ix0_owner:ixN_owner) => yslab_workspace(:, block_count + 1:2*block_count)
    pressure_vslab_v(-1:ny + 1, -nz:nz, ix0_owner:ixN_owner) => yslab_workspace(:, 2*block_count + 1:3*block_count)
    pressure_vslab_w(-1:ny + 1, -nz:nz, ix0_owner:ixN_owner) => yslab_workspace(:, 3*block_count + 1:4*block_count)
  end subroutine transpose_pressure_velocity_to_yslab

  subroutine pressure_velocity_view(owner_src, u_owner, v_owner, w_owner, ix0_owner, ixN_owner)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    complex(C_DOUBLE_COMPLEX), pointer :: u_owner(:, :, :), v_owner(:, :, :), w_owner(:, :, :)
    integer(C_INT), intent(out) :: ix0_owner, ixN_owner

    ix0_owner = lbound(owner_src, 3)
    ixN_owner = ubound(owner_src, 3)

    if (use_yslab_linsolve) then
      u_owner => pressure_vslab_u
      v_owner => pressure_vslab_v
      w_owner => pressure_vslab_w
    else
      u_owner(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN) => V(:, :, :, 1)
      v_owner(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN) => V(:, :, :, 2)
      w_owner(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN) => V(:, :, :, 3)
    end if
  end subroutine pressure_velocity_view

  SUBROUTINE solve_pressure_field(src0, src1, p)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: p(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    call assemble_pressure_rhs(src0, src1, p)
    call solve_compact_component_current_layout(p, assemble_pressure_operator, assemble_pressure_boundaries, 0.0d0, 0.0d0, p, &
                                 solve_label="pressure current-layout gpsv", symmetric_operator=.true., transpose_derivative=.true.)
  END SUBROUTINE solve_pressure_field

  SUBROUTINE solve_dpdy_field(src0, src1, dpdy)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dpdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    call assemble_dpdy_rhs(src0, src1, dpdy)
    call solve_compact_component_current_layout(dpdy, assemble_dpdy_operator, assemble_dpdy_boundaries, 0.0d0, 0.0d0, dpdy, &
                                                solve_label="pressure current-layout gpsv", symmetric_operator=.true.)
  END SUBROUTINE solve_dpdy_field

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

  subroutine assemble_pressure_operator(owner_src, lambda_coeff, diffusion_coeff, row_start, row_end)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff
    integer(C_INT), intent(in) :: row_start, row_end
    integer(C_INT) :: ix, iz, iy, ix0_owner, ixN_owner

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

  subroutine assemble_dpdy_operator(owner_src, lambda_coeff, diffusion_coeff, row_start, row_end)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff
    integer(C_INT), intent(in) :: row_start, row_end
    integer(C_INT) :: ix, iz, iy, ix0_owner, ixN_owner

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

  subroutine assemble_pressure_boundaries(owner_src, row_start, row_end, has_lower_boundary, has_upper_boundary)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    integer(C_INT), intent(in) :: row_start, row_end
    logical, intent(in) :: has_lower_boundary, has_upper_boundary
    integer(C_INT) :: ix, iz, ix0_owner, ixN_owner
    complex(C_DOUBLE_COMPLEX), pointer :: u_owner(:, :, :), v_owner(:, :, :), w_owner(:, :, :)

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

  subroutine assemble_dpdy_boundaries(owner_src, row_start, row_end, has_lower_boundary, has_upper_boundary)
    implicit none
    complex(C_DOUBLE_COMPLEX), pointer, intent(in) :: owner_src(:, :, :)
    integer(C_INT), intent(in) :: row_start, row_end
    logical, intent(in) :: has_lower_boundary, has_upper_boundary
    integer(C_INT) :: ix, iz, ix0_owner, ixN_owner
    complex(C_DOUBLE_COMPLEX), pointer :: u_owner(:, :, :), v_owner(:, :, :), w_owner(:, :, :)

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
