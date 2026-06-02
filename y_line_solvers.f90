#include "header.h"

module y_line_solvers

  use, intrinsic :: iso_c_binding
  use mpi_transpose, only: ny0, nx0, nxN, nxB, yl0, ylN, ylB, npy_grid, ierr, ipy, MPI_COMM_Y, &
                           allgather_y_blocks_to_xz_full
#ifdef HAVE_MPI
  use mpi_f08
#endif

  implicit none
  private

  public :: ys_lu5decomp, ys_leftlu5div
  public :: ys_solve_compact_derivative, ys_solve_compact_system, ys_solve_ghost_system
  public :: ys_solve_ghost_field, ys_solve_ghost_field_reduced

  abstract interface
    subroutine ys_build_ghost_line(ix_global, iz_global, src0_line, src1_line, line, a, eqm1, eq0, eqn, eqnp1, ny)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(C_INT), intent(in) :: ix_global, iz_global, ny
      complex(C_DOUBLE_COMPLEX), intent(in) :: src0_line(:), src1_line(:)
      complex(C_DOUBLE_COMPLEX), intent(out) :: line(-1:ny + 1)
      real(C_DOUBLE), intent(out) :: a(1:ny + 1, -2:2)
      real(C_DOUBLE), intent(out) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    end subroutine ys_build_ghost_line
  end interface

contains

  subroutine ys_solve_ghost_field(src0, src1, dst, ny, nz, build_line)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(:, :, :), src1(:, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dst(:, :, :)
    procedure(ys_build_ghost_line) :: build_line

    ! The reduced interface solve keeps two boundary-adjacent unknowns per side,
    ! so each rank must own at least four physical/ghost rows.
    if (npy_grid > 1 .and. ylB < 4) error stop "ys_solve_ghost_field requires at least four y rows per rank"
    call ys_solve_ghost_field_reduced(src0, src1, dst, ny, nz, build_line)
  end subroutine ys_solve_ghost_field

  subroutine ys_lu5decomp(a)
    real(C_DOUBLE), intent(inout) :: a(0:, -2:)
    integer(C_INT) :: hi1, hi2
    real(C_DOUBLE) :: piv
    integer :: i, k, j

    hi1 = size(a, 1) - 1
    hi2 = size(a, 2) - 3
    a(hi1 - 2, 1:2) = 0
    a(hi1 - 3, 2) = 0
    do i = hi1 - hi2, 0, -1
      do k = hi2, 1, -1
        piv = a(i, k)
        do j = -1, -2, -1
          a(i, j + k) = a(i, j + k) - piv*a(i + k, j)
        end do
      end do
      piv = 1.0d0/a(i, 0)
      a(i, 0) = piv
      a(i, -2:-1) = a(i, -2:-1)*piv
    end do
    a(0, -2:-1) = 0
    a(1, -2) = 0
  end subroutine ys_lu5decomp

  subroutine ys_leftlu5div(x, a)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(-2:)
    real(C_DOUBLE), intent(in) :: a(0:, -2:)
    integer(C_INT) :: hi1, hi2, i

    hi1 = size(a, 1) - 1
    hi2 = size(a, 2) - 3

    do i = hi1 - hi2, 0, -1
      x(i) = x(i) - (a(i, 1)*x(i + 1) + a(i, 2)*x(i + 2))
      x(i) = x(i)*a(i, 0)
    end do

    do i = 0, hi1
      x(i) = x(i) - (a(i, -2)*x(i - 2) + a(i, -1)*x(i - 1))
    end do
  end subroutine ys_leftlu5div

  subroutine ys_solve_compact_derivative(f0, f1, der, d0mat, d140, d14m1, d14n, d14np1, ny, ny0, nyN)
    integer(C_INT), intent(in) :: ny, ny0, nyN
    complex(C_DOUBLE_COMPLEX), intent(in) :: f0(-1:ny + 1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: f1(-1:ny + 1)
    real(C_DOUBLE), intent(in) :: der(ny0:nyN, 0:3, -2:2)
    real(C_DOUBLE), intent(in) :: d0mat(ny0:nyN + 2, -2:2)
    real(C_DOUBLE), intent(in) :: d140(-2:2), d14m1(-2:2), d14n(-2:2), d14np1(-2:2)
    integer(C_INT) :: iy

    f1(0) = sum(d140(-2:2)*f0(-1:3))
    f1(-1) = sum(d14m1(-2:2)*f0(-1:3))
    f1(ny) = sum(d14n(-2:2)*f0(ny - 3:ny + 1))
    f1(ny + 1) = sum(d14np1(-2:2)*f0(ny - 3:ny + 1))
    do iy = ny0, nyN
      f1(iy) = sum(der(iy, 1, -2:2)*f0(iy - 2:iy + 2))
    end do
    f1(1) = f1(1) - (der(1, 0, -1)*f1(0) + der(1, 0, -2)*f1(-1))
    f1(2) = f1(2) - der(2, 0, -2)*f1(0)
    f1(ny - 1) = f1(ny - 1) - (der(ny - 1, 0, 1)*f1(ny) + der(ny - 1, 0, 2)*f1(ny + 1))
    f1(ny - 2) = f1(ny - 2) - der(ny - 2, 0, 2)*f1(ny)
    call ys_leftlu5div(f1, d0mat)
  end subroutine ys_solve_compact_derivative

  subroutine ys_solve_compact_system(x, a, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                     rhs_lower, rhs_lower_ghost, rhs_upper, rhs_upper_ghost, ny, ny0, nyN)
    integer(C_INT), intent(in) :: ny, ny0, nyN
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(-1:ny + 1)
    real(C_DOUBLE), intent(inout) :: a(ny0:nyN + 2, -2:2)
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: rhs_lower, rhs_lower_ghost, rhs_upper, rhs_upper_ghost
    real(C_DOUBLE) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)

    eqm1 = lower_ghost_bc
    eq0 = lower_bc
    eqn = upper_bc
    eqnp1 = upper_ghost_bc
    x(-1) = rhs_lower_ghost
    x(0) = rhs_lower
    x(ny) = rhs_upper
    x(ny + 1) = rhs_upper_ghost
    call ys_solve_ghost_system(x, a, eqm1, eq0, eqn, eqnp1, ny)
  end subroutine ys_solve_compact_system

  subroutine ys_solve_ghost_system(x, a, eqm1, eq0, eqn, eqnp1, ny)
    integer(C_INT), intent(in) :: ny
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(-1:ny + 1)
    real(C_DOUBLE), intent(inout) :: a(1:ny + 1, -2:2)
    real(C_DOUBLE), intent(inout) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)

    x(0) = x(0) - x(-1)*eq0(-2)/eqm1(-2)
    eq0(-2:2) = eq0(-2:2) - eqm1(-2:2)*eq0(-2)/eqm1(-2)
    eq0(-2) = 0.0d0

    x(1) = x(1) - x(-1)*a(1, -2)/eqm1(-2)
    a(1, -2:2) = a(1, -2:2) - eqm1(-2:2)*a(1, -2)/eqm1(-2)
    a(1, -2) = 0.0d0

    x(1) = x(1) - x(0)*a(1, -1)/eq0(-1)
    a(1, -2:2) = a(1, -2:2) - eq0(-2:2)*a(1, -1)/eq0(-1)
    a(1, -1) = 0.0d0

    x(2) = x(2) - x(0)*a(2, -2)/eq0(-1)
    a(2, -2:1) = a(2, -2:1) - eq0(-1:2)*a(2, -2)/eq0(-1)
    a(2, -2) = 0.0d0

    x(ny) = x(ny) - x(ny + 1)*eqn(2)/eqnp1(2)
    eqn(-2:2) = eqn(-2:2) - eqnp1(-2:2)*eqn(2)/eqnp1(2)
    eqn(2) = 0.0d0

    x(ny - 1) = x(ny - 1) - x(ny + 1)*a(ny - 1, 2)/eqnp1(2)
    a(ny - 1, -2:2) = a(ny - 1, -2:2) - eqnp1(-2:2)*a(ny - 1, 2)/eqnp1(2)
    a(ny - 1, 2) = 0.0d0

    x(ny - 1) = x(ny - 1) - x(ny)*a(ny - 1, 1)/eqn(1)
    a(ny - 1, -2:2) = a(ny - 1, -2:2) - eqn(-2:2)*a(ny - 1, 1)/eqn(1)
    a(ny - 1, 1) = 0.0d0

    x(ny - 2) = x(ny - 2) - x(ny)*a(ny - 2, 2)/eqn(1)
    a(ny - 2, -1:2) = a(ny - 2, -1:2) - eqn(-2:1)*a(ny - 2, 2)/eqn(1)
    a(ny - 2, 2) = 0.0d0

    call ys_lu5decomp(a)
    call ys_leftlu5div(x, a)

    x(0) = (x(0) - sum(eq0(0:2)*x(1:3)))/eq0(-1)
    x(-1) = (x(-1) - sum(eqm1(-1:2)*x(0:3)))/eqm1(-2)
    x(ny) = (x(ny) - sum(eqn(-2:0)*x(ny - 3:ny - 1)))/eqn(1)
    x(ny + 1) = (x(ny + 1) - sum(eqnp1(-2:1)*x(ny - 3:ny)))/eqnp1(2)
  end subroutine ys_solve_ghost_system

  subroutine ys_solve_ghost_field_reduced(src0, src1, dst, ny, nz, build_line)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(:, :, :), src1(:, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dst(:, :, :)
    procedure(ys_build_ghost_line) :: build_line
    complex(C_DOUBLE_COMPLEX), allocatable :: dst_xz(:, :, :), dst_xz_full(:, :, :)
    complex(C_DOUBLE_COMPLEX) :: line(-1:ny + 1), local_line(-1:ny + 1)
    complex(C_DOUBLE_COMPLEX) :: local_block(ylB)
    complex(C_DOUBLE_COMPLEX) :: packed_remote(20), u_left(2), u_right(2)
    real(C_DOUBLE), allocatable :: factored_a(:, :)
    complex(C_DOUBLE_COMPLEX), allocatable :: rhs_base(:), left_couple(:, :), right_couple(:, :)
    complex(C_DOUBLE_COMPLEX), allocatable :: reduced_rhs(:, :), reduced_store(:, :, :)
    complex(C_DOUBLE_COMPLEX), allocatable :: packed_send(:, :), packed_recv(:, :, :), mat(:, :), rhs(:), sol(:)
    complex(C_DOUBLE_COMPLEX) :: lower_rhsm1, lower_rhs0, upper_rhsn, upper_rhsnp1
    real(C_DOUBLE) :: a(1:ny + 1, -2:2), eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    real(C_DOUBLE) :: local_a(1:ny + 1, -2:2), local_eqm1(-2:2), local_eq0(-2:2), local_eqn(-2:2), local_eqnp1(-2:2)
    real(C_DOUBLE) :: lower_eq0(-2:2), upper_eqn(-2:2)
    integer(C_INT) :: ix, iz, iline, nlines, nlines_z, row_start, row_end, active_n, row, iblock, row0
    logical :: has_lower_boundary, has_upper_boundary
    real(C_DOUBLE), allocatable :: eqm1_store(:, :), eqnp1_store(:, :), lower_eq0_store(:, :), upper_eqn_store(:, :)
    complex(C_DOUBLE_COMPLEX), allocatable :: lower_rhsm1_store(:), lower_rhs0_store(:), upper_rhsn_store(:), upper_rhsnp1_store(:)

    ! Solve each x/z line on the native distributed-y layout, then gather the
    ! owned y blocks back into the full x/z-shaped output field.
    allocate (dst_xz(ylB, 2*nz + 1, nxB))
    allocate (dst_xz_full(ny + 3, 2*nz + 1, nxB))

    row_start = max(1_C_INT, yl0 - 2)
    row_end = min(ny - 1, ylN - 2)
    active_n = row_end - row_start + 1
    has_lower_boundary = (row_start == 1)
    has_upper_boundary = (row_end == ny - 1)

    if (npy_grid == 1 .or. active_n < 4) then
      do ix = nx0, nxN
        do iz = -nz, nz
          call build_line(ix, iz, src0(:, iz + nz + 1, ix - nx0 + 1), src1(:, iz + nz + 1, ix - nx0 + 1), &
                          line, a, eqm1, eq0, eqn, eqnp1, ny)
          local_line = line
          local_a = a
          local_eqm1 = eqm1
          local_eq0 = eq0
          local_eqn = eqn
          local_eqnp1 = eqnp1
          call ys_solve_ghost_system(local_line, local_a, local_eqm1, local_eq0, local_eqn, local_eqnp1, ny)
          dst_xz(:, iz + nz + 1, ix - nx0 + 1) = local_line(yl0 - 2:ylN - 2)
        end do
      end do
    else
      nlines_z = size(src0, 2)
      nlines = size(src0, 3)*nlines_z
      allocate (factored_a(0:active_n - 1, -2:2), rhs_base(0:active_n - 1))
      allocate (left_couple(0:active_n - 1, 2), right_couple(0:active_n - 1, 2))
      allocate (packed_send(20, nlines), packed_recv(20, nlines, npy_grid))
      allocate (reduced_rhs(0:active_n - 1, 5), reduced_store(0:active_n - 1, 5, nlines))
      allocate (eqm1_store(-2:2, nlines), eqnp1_store(-2:2, nlines), lower_eq0_store(-2:2, nlines), upper_eqn_store(-2:2, nlines))
      allocate (lower_rhsm1_store(nlines), lower_rhs0_store(nlines), upper_rhsn_store(nlines), upper_rhsnp1_store(nlines))
      allocate (mat(4*npy_grid, 4*npy_grid), rhs(4*npy_grid), sol(4*npy_grid))

      do ix = nx0, nxN
        do iz = -nz, nz
          iline = (ix - nx0)*nlines_z + (iz + nz + 1)

          call build_line(ix, iz, src0(:, iz + nz + 1, ix - nx0 + 1), src1(:, iz + nz + 1, ix - nx0 + 1), &
                          line, a, eqm1, eq0, eqn, eqnp1, ny)
          call ys_build_reduced_local_block( &
            line, a, eqm1, eq0, eqn, eqnp1, row_start, row_end, active_n, ny, &
            has_lower_boundary, has_upper_boundary, factored_a, rhs_base, left_couple, right_couple, &
            lower_rhsm1, lower_rhs0, lower_eq0, upper_rhsn, upper_rhsnp1, upper_eqn)
          call ys_factor_penta(factored_a)
          ! One block solve produces the particular solution plus the four
          ! interface response vectors needed by the reduced system.
          reduced_rhs(:, 1) = rhs_base
          reduced_rhs(:, 2) = -left_couple(:, 1)
          reduced_rhs(:, 3) = -left_couple(:, 2)
          reduced_rhs(:, 4) = -right_couple(:, 1)
          reduced_rhs(:, 5) = -right_couple(:, 2)
          call ys_solve_factored_penta_multi(reduced_rhs, factored_a)
          reduced_store(:, :, iline) = reduced_rhs
          eqm1_store(:, iline) = eqm1
          eqnp1_store(:, iline) = eqnp1
          lower_eq0_store(:, iline) = lower_eq0
          upper_eqn_store(:, iline) = upper_eqn
          lower_rhsm1_store(iline) = lower_rhsm1
          lower_rhs0_store(iline) = lower_rhs0
          upper_rhsn_store(iline) = upper_rhsn
          upper_rhsnp1_store(iline) = upper_rhsnp1

          packed_send(1:4, iline) = [reduced_rhs(0, 1), reduced_rhs(1, 1), &
                                     reduced_rhs(active_n - 2, 1), reduced_rhs(active_n - 1, 1)]
          packed_send(5:12, iline) = [reduced_rhs(0, 2), reduced_rhs(1, 2), &
                                      reduced_rhs(active_n - 2, 2), reduced_rhs(active_n - 1, 2), &
                                      reduced_rhs(0, 3), reduced_rhs(1, 3), &
                                      reduced_rhs(active_n - 2, 3), reduced_rhs(active_n - 1, 3)]
          packed_send(13:20, iline) = [reduced_rhs(0, 4), reduced_rhs(1, 4), &
                                       reduced_rhs(active_n - 2, 4), reduced_rhs(active_n - 1, 4), &
                                       reduced_rhs(0, 5), reduced_rhs(1, 5), &
                                       reduced_rhs(active_n - 2, 5), reduced_rhs(active_n - 1, 5)]
        end do
      end do

#ifdef HAVE_MPI
      call MPI_Allgather(packed_send, 20*nlines, MPI_DOUBLE_COMPLEX, packed_recv, 20*nlines, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr)

      do ix = nx0, nxN
        do iz = -nz, nz
          iline = (ix - nx0)*nlines_z + (iz + nz + 1)
          mat = (0.0d0, 0.0d0)
          rhs = (0.0d0, 0.0d0)
          do iblock = 0, npy_grid - 1
            row0 = 4*iblock
            packed_remote = packed_recv(:, iline, iblock + 1)
            mat(row0 + 1:row0 + 4, row0 + 1:row0 + 4) = (0.0d0, 0.0d0)
            mat(row0 + 1, row0 + 1) = (1.0d0, 0.0d0)
            mat(row0 + 2, row0 + 2) = (1.0d0, 0.0d0)
            mat(row0 + 3, row0 + 3) = (1.0d0, 0.0d0)
            mat(row0 + 4, row0 + 4) = (1.0d0, 0.0d0)
            rhs(row0 + 1:row0 + 4) = packed_remote(1:4)
            if (iblock > 0) then
              mat(row0 + 1:row0 + 4, row0 - 1) = -packed_remote(5:8)
              mat(row0 + 1:row0 + 4, row0) = -packed_remote(9:12)
            end if
            if (iblock < npy_grid - 1) then
              mat(row0 + 1:row0 + 4, row0 + 5) = -packed_remote(13:16)
              mat(row0 + 1:row0 + 4, row0 + 6) = -packed_remote(17:20)
            end if
          end do

          sol = rhs
          call ys_solve_dense_complex(mat, sol)
          row0 = 4*ipy
          u_left = (0.0d0, 0.0d0)
          u_right = (0.0d0, 0.0d0)
          if (ipy > 0) u_left = sol(row0 - 1:row0)
          if (ipy < npy_grid - 1) u_right = sol(row0 + 5:row0 + 6)

          reduced_rhs(:, 1) = reduced_store(:, 1, iline) + reduced_store(:, 2, iline)*u_left(1) + &
                              reduced_store(:, 3, iline)*u_left(2) + reduced_store(:, 4, iline)*u_right(1) + &
                              reduced_store(:, 5, iline)*u_right(2)

          local_block = (0.0d0, 0.0d0)
          do row = row_start, row_end
            local_block(row - (yl0 - 2) + 1) = reduced_rhs(row - row_start, 1)
          end do

          ! Reconstruct the wall and ghost rows in the same order as
          ! ys_solve_ghost_system once the local interior block is known.
          if (has_lower_boundary) then
            local_block(2) = (lower_rhs0_store(iline) - &
                              sum(lower_eq0_store(0:2, iline)*[local_block(3), local_block(4), local_block(5)]))/ &
                             lower_eq0_store(-1, iline)
            local_block(1) = (lower_rhsm1_store(iline) - &
                              sum(eqm1_store(-1:2, iline)*[local_block(2), local_block(3), local_block(4), local_block(5)]))/ &
                             eqm1_store(-2, iline)
          end if
          if (has_upper_boundary) then
            local_block(ny - (yl0 - 2) + 1) = (upper_rhsn_store(iline) - &
                                               sum(upper_eqn_store(-2:0, iline)*[ &
                                                   local_block(ny - 3 - (yl0 - 2) + 1), &
                                                   local_block(ny - 2 - (yl0 - 2) + 1), &
                                                   local_block(ny - 1 - (yl0 - 2) + 1)]))/upper_eqn_store(1, iline)
            local_block(ny + 1 - (yl0 - 2) + 1) = (upper_rhsnp1_store(iline) - &
                                                   sum(eqnp1_store(-2:1, iline)*[ &
                                                       local_block(ny - 3 - (yl0 - 2) + 1), &
                                                       local_block(ny - 2 - (yl0 - 2) + 1), &
                                                       local_block(ny - 1 - (yl0 - 2) + 1), &
                                                       local_block(ny - (yl0 - 2) + 1)]))/eqnp1_store(2, iline)
          end if

          dst_xz(:, iz + nz + 1, ix - nx0 + 1) = local_block
        end do
      end do
#endif

      deallocate (factored_a, rhs_base, left_couple, right_couple, packed_send, packed_recv, reduced_rhs, reduced_store, &
                  mat, rhs, sol)
      deallocate (eqm1_store, eqnp1_store, lower_eq0_store, upper_eqn_store)
      deallocate (lower_rhsm1_store, lower_rhs0_store, upper_rhsn_store, upper_rhsnp1_store)
    end if

    call allgather_y_blocks_to_xz_full(dst_xz, dst_xz_full)

    do ix = nx0, nxN
      do iz = -nz, nz
        dst(:, iz + nz + 1, ix - nx0 + 1) = dst_xz_full(:, iz + nz + 1, ix - nx0 + 1)
      end do
    end do

    deallocate (dst_xz, dst_xz_full)
  end subroutine ys_solve_ghost_field_reduced

  subroutine ys_build_reduced_local_block(line, a, eqm1, eq0, eqn, eqnp1, row_start, row_end, active_n, ny, &
                                          has_lower_boundary, has_upper_boundary, factored_a, rhs_base, left_couple, right_couple, &
                                          lower_rhsm1, lower_rhs0, lower_eq0, upper_rhsn, upper_rhsnp1, upper_eqn)
    integer(C_INT), intent(in) :: row_start, row_end, active_n, ny
    logical, intent(in) :: has_lower_boundary, has_upper_boundary
    complex(C_DOUBLE_COMPLEX), intent(in) :: line(-1:ny + 1)
    real(C_DOUBLE), intent(in) :: a(1:ny + 1, -2:2), eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    real(C_DOUBLE), intent(out) :: factored_a(0:active_n - 1, -2:2), lower_eq0(-2:2), upper_eqn(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: rhs_base(0:active_n - 1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: left_couple(0:active_n - 1, 2), right_couple(0:active_n - 1, 2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: lower_rhsm1, lower_rhs0, upper_rhsn, upper_rhsnp1
    real(C_DOUBLE) :: row_coeffs(-2:2)
    integer(C_INT) :: row, idx, col, global_col

    factored_a = 0.0d0
    rhs_base = (0.0d0, 0.0d0)
    left_couple = (0.0d0, 0.0d0)
    right_couple = (0.0d0, 0.0d0)

    ! Mirror the elimination order of ys_solve_ghost_system up to the point
    ! where the local interior block becomes a pure pentadiagonal solve.
    lower_rhsm1 = line(-1)
    lower_rhs0 = line(0) - line(-1)*eq0(-2)/eqm1(-2)
    lower_eq0 = eq0 - eqm1*eq0(-2)/eqm1(-2)
    lower_eq0(-2) = 0.0d0
    upper_rhsnp1 = line(ny + 1)
    upper_rhsn = line(ny) - line(ny + 1)*eqn(2)/eqnp1(2)
    upper_eqn = eqn - eqnp1*eqn(2)/eqnp1(2)
    upper_eqn(2) = 0.0d0

    do row = row_start, row_end
      rhs_base(row - row_start) = line(row)
    end do

    if (has_lower_boundary) then
      rhs_base(0) = rhs_base(0) - line(-1)*a(1, -2)/eqm1(-2)
      row_coeffs = a(1, -2:2) - eqm1*a(1, -2)/eqm1(-2)
      row_coeffs(-2) = 0.0d0
      rhs_base(0) = rhs_base(0) - lower_rhs0*row_coeffs(-1)/lower_eq0(-1)
      row_coeffs = row_coeffs - lower_eq0*row_coeffs(-1)/lower_eq0(-1)
      row_coeffs(-1) = 0.0d0
      factored_a(0, 0:2) = row_coeffs(0:2)

      rhs_base(1) = rhs_base(1) - lower_rhs0*a(2, -2)/lower_eq0(-1)
      row_coeffs = a(2, -2:2)
      row_coeffs(-2:1) = row_coeffs(-2:1) - lower_eq0(-1:2)*a(2, -2)/lower_eq0(-1)
      row_coeffs(-2) = 0.0d0
      factored_a(1, -1:2) = row_coeffs(-1:2)
    end if

    if (has_upper_boundary) then
      rhs_base(active_n - 1) = rhs_base(active_n - 1) - line(ny + 1)*a(ny - 1, 2)/eqnp1(2)
      row_coeffs = a(ny - 1, -2:2) - eqnp1*a(ny - 1, 2)/eqnp1(2)
      row_coeffs(2) = 0.0d0
      rhs_base(active_n - 1) = rhs_base(active_n - 1) - upper_rhsn*row_coeffs(1)/upper_eqn(1)
      row_coeffs = row_coeffs - upper_eqn*row_coeffs(1)/upper_eqn(1)
      row_coeffs(1) = 0.0d0
      factored_a(active_n - 1, -2:0) = row_coeffs(-2:0)

      rhs_base(active_n - 2) = rhs_base(active_n - 2) - upper_rhsn*a(ny - 2, 2)/upper_eqn(1)
      row_coeffs = a(ny - 2, -2:2)
      row_coeffs(-1:2) = row_coeffs(-1:2) - upper_eqn(-2:1)*a(ny - 2, 2)/upper_eqn(1)
      row_coeffs(2) = 0.0d0
      factored_a(active_n - 2, -2:1) = row_coeffs(-2:1)
    end if

    do row = row_start, row_end
      idx = row - row_start
      if ((has_lower_boundary .and. row <= 2) .or. (has_upper_boundary .and. row >= ny - 2)) cycle
      row_coeffs = a(row, -2:2)
      do col = -2, 2
        global_col = row + col
        if (global_col < row_start) then
          if (global_col == row_start - 2) left_couple(idx, 1) = cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
          if (global_col == row_start - 1) left_couple(idx, 2) = cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
        else if (global_col > row_end) then
          if (global_col == row_end + 1) right_couple(idx, 1) = cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
          if (global_col == row_end + 2) right_couple(idx, 2) = cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
        else
          factored_a(idx, col) = row_coeffs(col)
        end if
      end do
    end do
  end subroutine ys_build_reduced_local_block

  subroutine ys_solve_dense_complex(mat, rhs)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: mat(:, :)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX) :: factor, pivot_row(size(mat, 2)), rhs_tmp
    real(C_DOUBLE) :: pivot_abs, cand_abs
    integer(C_INT) :: n, i, j, pivot

    n = int(size(rhs), C_INT)
    do i = 1, n
      pivot = i
      pivot_abs = abs(mat(i, i))
      do j = i + 1, n
        cand_abs = abs(mat(j, i))
        if (cand_abs > pivot_abs) then
          pivot = j
          pivot_abs = cand_abs
        end if
      end do
      if (pivot /= i) then
        pivot_row = mat(i, :)
        mat(i, :) = mat(pivot, :)
        mat(pivot, :) = pivot_row
        rhs_tmp = rhs(i)
        rhs(i) = rhs(pivot)
        rhs(pivot) = rhs_tmp
      end if
      do j = i + 1, n
        factor = mat(j, i)/mat(i, i)
        mat(j, i:n) = mat(j, i:n) - factor*mat(i, i:n)
        rhs(j) = rhs(j) - factor*rhs(i)
      end do
    end do
    do i = n, 1, -1
      if (i < n) rhs(i) = rhs(i) - sum(mat(i, i + 1:n)*rhs(i + 1:n))
      rhs(i) = rhs(i)/mat(i, i)
    end do
  end subroutine ys_solve_dense_complex

  subroutine ys_factor_penta(a)
    real(C_DOUBLE), intent(inout) :: a(0:, -2:)
    integer(C_INT) :: n, i
    real(C_DOUBLE) :: piv, factor

    ! Generic pentadiagonal LU on the local reduced block. The wall/ghost rows
    ! are already eliminated outside this factorization, so this is a pure
    ! interior solve rather than the boundary-aware variant used elsewhere.
    n = size(a, 1)
    do i = 0, n - 1
      piv = a(i, 0)
      a(i, 0) = 1.0d0/piv

      if (i + 1 < n) then
        factor = a(i + 1, -1)*a(i, 0)
        a(i + 1, -1) = factor
        a(i + 1, 0) = a(i + 1, 0) - factor*a(i, 1)
        if (i + 2 < n) a(i + 1, 1) = a(i + 1, 1) - factor*a(i, 2)
      end if

      if (i + 2 < n) then
        factor = a(i + 2, -2)*a(i, 0)
        a(i + 2, -2) = factor
        a(i + 2, -1) = a(i + 2, -1) - factor*a(i, 1)
        a(i + 2, 0) = a(i + 2, 0) - factor*a(i, 2)
      end if
    end do
  end subroutine ys_factor_penta

  subroutine ys_solve_factored_penta_multi(rhs, a)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(0:, :)
    real(C_DOUBLE), intent(in) :: a(0:, -2:)
    integer(C_INT) :: n, nrhs, i

    n = size(a, 1)
    nrhs = size(rhs, 2)

    do i = 0, n - 1
      if (i >= 2) rhs(i, 1:nrhs) = rhs(i, 1:nrhs) - a(i, -2)*rhs(i - 2, 1:nrhs)
      if (i >= 1) rhs(i, 1:nrhs) = rhs(i, 1:nrhs) - a(i, -1)*rhs(i - 1, 1:nrhs)
    end do

    do i = n - 1, 0, -1
      if (i + 1 < n) rhs(i, 1:nrhs) = rhs(i, 1:nrhs) - a(i, 1)*rhs(i + 1, 1:nrhs)
      if (i + 2 < n) rhs(i, 1:nrhs) = rhs(i, 1:nrhs) - a(i, 2)*rhs(i + 2, 1:nrhs)
      rhs(i, 1:nrhs) = rhs(i, 1:nrhs)*a(i, 0)
    end do
  end subroutine ys_solve_factored_penta_multi

end module y_line_solvers
