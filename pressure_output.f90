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
    real(C_DOUBLE) :: pmat(ny0:nyN + 2, -2:2), eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    complex(C_DOUBLE_COMPLEX) :: sol_solve(-1:ny + 1), tmp, tmp2
    integer(C_INT) :: ix, iz, iy

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(src0, src1, p, der, ny, ny0, nyN, nx0, nxN, nz) &
    !$omp private(ix, iz, iy)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0, nyN
          p(iy, iz, ix) = sum(der(iy, 0, -2:2)*src0(iy - 2:iy + 2, iz, ix)) + &
                          sum(der(iy, 1, -2:2)*src1(iy - 2:iy + 2, iz, ix))
        end do
      end do
    end do

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(src0, src1, p, der, k2, d140, d240, d24n, V, ni, ialfa, ibeta, ny, ny0, nyN, nx0, nxN, nz) &
    !$omp private(ix, iz, iy, pmat, eqm1, eq0, eqn, eqnp1, sol_solve, tmp, tmp2)
    do ix = nx0, nxN
      do iz = -nz, nz

        pmat = 0

        if (ix == 0 .and. iz == 0) then

          ! LHS = d2p/dy2
          do iy = ny0, nyN
            pmat(iy, -2:2) = der(iy, 2, -2:2)
          end do

          ! Bottom Ghost BC: momentum_y
          eqm1(:) = d140
          p(-1, iz, ix) = ni*sum(d240(-2:2)*V(-1:3, iz, ix, 2))

          ! Bottom BC: usual poisson equation in the vertical direction
          eq0(:) = d240
          p(0, iz, ix) = src0(0, iz, ix) + src1(0, iz, ix)

          ! Top  BC: Dirichlet p(ny)=0
          p(ny, iz, ix) = 0.0d0
          eqn(:) = 0.d0; eqn(1) = 1.d0

          ! Bottom Ghost BC, ghost not needed d4(p)|ny-1 = 0
          eqnp1(:) = der(ny - 1, 3, :)
          p(ny + 1, iz, ix) = 0.0d0

        else

          ! LHS = lapl(p) = D2(p)- k2*D0(p)
          do iy = ny0, nyN
            pmat(iy, -2:2) = der(iy, 2, -2:2) - k2(iz, ix)*der(iy, 0, -2:2)
          end do

          ! Bottom Ghost BC, ghost not needed d4(p)|1 = 0
          eqm1(:) = der(1, 3, :)
          p(-1, iz, ix) = 0.0d0

          ! Bottom BC, dirichlet from d/dx(momentum_x) + d/dz(momentum_z) solved for p
          tmp = sum(d240(-2:2)*V(-1:3, iz, ix, 1))
          tmp2 = sum(d240(-2:2)*V(-1:3, iz, ix, 3))
          p(0, iz, ix) = -ni*(ialfa(ix)*tmp + ibeta(iz)*tmp2)/k2(iz, ix)
          eq0(:) = 0.d0; eq0(-1) = 1.d0

          ! Top BC, dirichlet from d/dx(momentum_x) + d/dz(momentum_z) solved for p
          tmp = sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, 1))
          tmp2 = sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, 3))
          p(ny, iz, ix) = -ni*(ialfa(ix)*tmp + ibeta(iz)*tmp2)/k2(iz, ix)
          eqn(:) = 0.d0; eqn(1) = 1.d0

          ! Top Ghost BC, ghost not needed d4(p)|ny-1 = 0
          p(ny + 1, iz, ix) = 0.0d0
          eqnp1(:) = der(ny - 1, 3, :)

        end if

        call ys_solve_ghost_system(p(:, iz, ix), pmat, eqm1, eq0, eqn, eqnp1, ny)

      end do
    end do
    !$omp target update from(p)
  END SUBROUTINE solve_pressure_field

  SUBROUTINE solve_dpdy_field(src0, src1, dpdy)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: src1(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dpdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    real(C_DOUBLE) :: pmat(ny0:nyN + 2, -2:2), eqm1(-2:2), eq0(-2:2), eqnp1(-2:2), eqn(-2:2)
    integer(C_INT) :: ix, iz, iy

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(src0, src1, dpdy, der, k2, d140, d240, d24n, V, ni, ialfa, ibeta, ny, ny0, nyN, nx0, nxN, nz) &
    !$omp private(ix, iz, iy, pmat, eqm1, eqnp1, eq0, eqn)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0, nyN
          pmat(iy, -2:2) = der(iy, 2, -2:2) - k2(iz, ix)*der(iy, 0, -2:2)
          dpdy(iy, iz, ix) = sum(der(iy, 1, -2:2)*src0(iy - 2:iy + 2, iz, ix)) + &
                             sum(der(iy, 2, -2:2)*src1(iy - 2:iy + 2, iz, ix))
        end do

        ! Bottom Dirichlet B.C. : Wall value q(0) = dp/dy|bottom from wall-normal momentum
        ! q = nu*d2v/dy2 for nonzero Fourier modes, and q = 0 for the mean mode.
        dpdy(0, iz, ix) = ni*sum(d240(-2:2)*V(-1:3, iz, ix, 2))
        eq0 = 0.d0; eq0(-1) = 1.d0

        ! Bottom ghost closure. This is not a physical BC; it is the d4(q)=0-style
        ! stencil used to remove q(-1) from the first interior row.
        dpdy(-1, iz, ix) = 0.0d0
        eqm1 = der(1, 3, :)

        ! Top Dirichlet B.C.: Wall value q(ny) = dp/dy|top from wall-normal momentum
        ! q = nu*d2v/dy2 for nonzero Fourier modes, and q = 0 for the mean mode.
        dpdy(ny, iz, ix) = ni*sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, 2))
        eqn = 0.d0; eqn(1) = 1.d0

        ! Top ghost closure. This is the numerical relation used to eliminate q(ny+1)
        ! from the last interior row before solving.
        dpdy(ny + 1, iz, ix) = 0.0d0
        eqnp1 = der(ny - 1, 3, :)

        call ys_solve_ghost_system(dpdy(:, iz, ix), pmat, eqm1, eq0, eqn, eqnp1, ny)

      end do
    end do
    !$omp target update from(dpdy)
  END SUBROUTINE solve_dpdy_field

END MODULE pressure_output
