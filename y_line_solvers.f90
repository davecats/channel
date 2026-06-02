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
    complex(C_DOUBLE_COMPLEX) :: line(-1:ny + 1)
    complex(C_DOUBLE_COMPLEX) :: local_block(ylB)
    real(C_DOUBLE) :: a(1:ny + 1, -2:2), eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    integer(C_INT) :: ix, iz

    ! Solve each x/z line on the native distributed-y layout, then gather the
    ! owned y blocks back into the full x/z-shaped output field.
    allocate (dst_xz(ylB, 2*nz + 1, nxB))
    allocate (dst_xz_full(ny + 3, 2*nz + 1, nxB))

    do ix = nx0, nxN
      do iz = -nz, nz
        call build_line(ix, iz, src0(:, iz + nz + 1, ix - nx0 + 1), src1(:, iz + nz + 1, ix - nx0 + 1), &
                        line, a, eqm1, eq0, eqn, eqnp1, ny)
        call ys_solve_ghost_block_reduced(line, a, eqm1, eq0, eqn, eqnp1, ny, local_block)
        dst_xz(:, iz + nz + 1, ix - nx0 + 1) = local_block
      end do
    end do

    call allgather_y_blocks_to_xz_full(dst_xz, dst_xz_full)

    do ix = nx0, nxN
      do iz = -nz, nz
        dst(:, iz + nz + 1, ix - nx0 + 1) = dst_xz_full(:, iz + nz + 1, ix - nx0 + 1)
      end do
    end do

    deallocate (dst_xz, dst_xz_full)
  end subroutine ys_solve_ghost_field_reduced

  subroutine ys_solve_ghost_block_reduced(line, a, eqm1, eq0, eqn, eqnp1, ny, local_block)
    integer(C_INT), intent(in) :: ny
    complex(C_DOUBLE_COMPLEX), intent(in) :: line(-1:ny + 1)
    real(C_DOUBLE), intent(in) :: a(1:ny + 1, -2:2)
    real(C_DOUBLE), intent(in) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: local_block(ylB)
    integer(C_INT), parameter :: niface = 4
    complex(C_DOUBLE_COMPLEX) :: c_local(niface), tl_local(niface, 2), tr_local(niface, 2)
    complex(C_DOUBLE_COMPLEX) :: u_left(2), u_right(2)
    complex(C_DOUBLE_COMPLEX), allocatable :: factored_a(:, :)
    complex(C_DOUBLE_COMPLEX), allocatable :: rhs_local(:), rhs_base(:), left_couple(:, :), right_couple(:, :)
    complex(C_DOUBLE_COMPLEX) :: lower_rhsm1, lower_rhs0, upper_rhsn, upper_rhsnp1
    complex(C_DOUBLE_COMPLEX) :: row_left(2), row_right(2)
    complex(C_DOUBLE_COMPLEX), allocatable :: particular(:), left_resp(:, :), right_resp(:, :)
    real(C_DOUBLE) :: row_coeffs(-2:2), lower_eq0(-2:2), upper_eqn(-2:2)
    integer(C_INT) :: row_start, row_end, active_n, row, idx, col, global_col
    logical :: has_lower_boundary, has_upper_boundary
    complex(C_DOUBLE_COMPLEX), allocatable :: gathered(:), mat(:, :), rhs(:), sol(:)
    integer(C_INT) :: iblock, row0, recv0

    row_start = max(1_C_INT, yl0 - 2)
    row_end = min(ny - 1, ylN - 2)
    active_n = row_end - row_start + 1
    has_lower_boundary = (row_start == 1)
    has_upper_boundary = (row_end == ny - 1)
    if (npy_grid == 1 .or. active_n < 4) then
      call ys_solve_ghost_block_local(line, a, eqm1, eq0, eqn, eqnp1, ny, local_block)
      return
    end if

    allocate (factored_a(active_n, active_n), rhs_local(active_n), rhs_base(active_n), &
              left_couple(active_n, 2), right_couple(active_n, 2))
    allocate (particular(active_n), left_resp(active_n, 2), right_resp(active_n, 2))
    factored_a = (0.0d0, 0.0d0)
    rhs_local = (0.0d0, 0.0d0)
    left_couple = (0.0d0, 0.0d0)
    right_couple = (0.0d0, 0.0d0)
    lower_rhsm1 = line(-1)
    lower_rhs0 = line(0) - line(-1)*eq0(-2)/eqm1(-2)
    lower_eq0 = eq0 - eqm1*eq0(-2)/eqm1(-2)
    lower_eq0(-2) = 0.0d0
    upper_rhsnp1 = line(ny + 1)
    upper_rhsn = line(ny) - line(ny + 1)*eqn(2)/eqnp1(2)
    upper_eqn = eqn - eqnp1*eqn(2)/eqnp1(2)
    upper_eqn(2) = 0.0d0

    do row = row_start, row_end
      idx = row - row_start + 1
      rhs_local(idx) = line(row)
      row_coeffs = a(row, -2:2)
      row_left = (0.0d0, 0.0d0)
      row_right = (0.0d0, 0.0d0)

      if (has_lower_boundary) then
        if (row == 1) then
          rhs_local(idx) = rhs_local(idx) - line(-1)*a(1, -2)/eqm1(-2)
          row_coeffs = row_coeffs - eqm1*a(1, -2)/eqm1(-2)
          row_coeffs(-2) = 0.0d0
          rhs_local(idx) = rhs_local(idx) - lower_rhs0*row_coeffs(-1)/lower_eq0(-1)
          row_coeffs = row_coeffs - lower_eq0*row_coeffs(-1)/lower_eq0(-1)
          row_coeffs(-1) = 0.0d0
        else if (row == 2) then
          rhs_local(idx) = rhs_local(idx) - lower_rhs0*a(2, -2)/lower_eq0(-1)
          row_coeffs(-2:1) = row_coeffs(-2:1) - lower_eq0(-1:2)*a(2, -2)/lower_eq0(-1)
          row_coeffs(-2) = 0.0d0
        end if
      end if

      if (has_upper_boundary) then
        if (row == ny - 1) then
          rhs_local(idx) = rhs_local(idx) - line(ny + 1)*a(ny - 1, 2)/eqnp1(2)
          row_coeffs = row_coeffs - eqnp1*a(ny - 1, 2)/eqnp1(2)
          row_coeffs(2) = 0.0d0
          rhs_local(idx) = rhs_local(idx) - upper_rhsn*row_coeffs(1)/upper_eqn(1)
          row_coeffs = row_coeffs - upper_eqn*row_coeffs(1)/upper_eqn(1)
          row_coeffs(1) = 0.0d0
        else if (row == ny - 2) then
          rhs_local(idx) = rhs_local(idx) - upper_rhsn*a(ny - 2, 2)/upper_eqn(1)
          row_coeffs(-1:2) = row_coeffs(-1:2) - upper_eqn(-2:1)*a(ny - 2, 2)/upper_eqn(1)
          row_coeffs(2) = 0.0d0
        end if
      end if

      do col = -2, 2
        global_col = row + col
        if (global_col < row_start) then
          if (global_col == row_start - 2) row_left(1) = cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
          if (global_col == row_start - 1) row_left(2) = cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
        else if (global_col > row_end) then
          if (global_col == row_end + 1) row_right(1) = cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
          if (global_col == row_end + 2) row_right(2) = cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
        else
          factored_a(idx, global_col - row_start + 1) = cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
        end if
      end do
      left_couple(idx, :) = row_left
      right_couple(idx, :) = row_right
    end do

    ! The local solve now only sees rows 1:ny-1 after boundary/ghost
    ! elimination, so it is genuinely pentadiagonal again.
    call ys_factor_banded_dense(factored_a, 2_C_INT)

    rhs_base = rhs_local
    call ys_solve_factored_banded_dense(factored_a, rhs_local, 2_C_INT)
    particular = rhs_local
    rhs_local = rhs_base

    rhs_local = -left_couple(:, 1)
    call ys_solve_factored_banded_dense(factored_a, rhs_local, 2_C_INT)
    left_resp(:, 1) = rhs_local

    rhs_local = -left_couple(:, 2)
    call ys_solve_factored_banded_dense(factored_a, rhs_local, 2_C_INT)
    left_resp(:, 2) = rhs_local

    rhs_local = -right_couple(:, 1)
    call ys_solve_factored_banded_dense(factored_a, rhs_local, 2_C_INT)
    right_resp(:, 1) = rhs_local

    rhs_local = -right_couple(:, 2)
    call ys_solve_factored_banded_dense(factored_a, rhs_local, 2_C_INT)
    right_resp(:, 2) = rhs_local

    rhs_local = rhs_base

#ifdef HAVE_MPI
    c_local = [particular(1), particular(2), particular(active_n - 1), particular(active_n)]
    tl_local(:, 1) = [left_resp(1, 1), left_resp(2, 1), left_resp(active_n - 1, 1), left_resp(active_n, 1)]
    tl_local(:, 2) = [left_resp(1, 2), left_resp(2, 2), left_resp(active_n - 1, 2), left_resp(active_n, 2)]
    tr_local(:, 1) = [right_resp(1, 1), right_resp(2, 1), right_resp(active_n - 1, 1), right_resp(active_n, 1)]
    tr_local(:, 2) = [right_resp(1, 2), right_resp(2, 2), right_resp(active_n - 1, 2), right_resp(active_n, 2)]

    allocate (gathered(20*npy_grid))
    call MPI_Allgather([c_local, reshape(tl_local, [8]), reshape(tr_local, [8])], 20, MPI_DOUBLE_COMPLEX, &
                       gathered, 20, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr)

    ! The current backend replicates the tiny interface system on every rank.
    ! This keeps communication to interface data only; a future PCR step can
    ! replace this dense replicated solve without changing the local reduction.
    allocate (mat(niface*npy_grid, niface*npy_grid), rhs(niface*npy_grid), sol(niface*npy_grid))
    mat = (0.0d0, 0.0d0)
    rhs = (0.0d0, 0.0d0)
    do iblock = 0, npy_grid - 1
      row0 = iblock*niface
      recv0 = 20*iblock
      c_local = gathered(recv0 + 1:recv0 + 4)
      tl_local = reshape(gathered(recv0 + 5:recv0 + 12), [4, 2])
      tr_local = reshape(gathered(recv0 + 13:recv0 + 20), [4, 2])
      mat(row0 + 1:row0 + niface, row0 + 1:row0 + niface) = 0.0d0
      mat(row0 + 1, row0 + 1) = (1.0d0, 0.0d0)
      mat(row0 + 2, row0 + 2) = (1.0d0, 0.0d0)
      mat(row0 + 3, row0 + 3) = (1.0d0, 0.0d0)
      mat(row0 + 4, row0 + 4) = (1.0d0, 0.0d0)
      rhs(row0 + 1:row0 + niface) = c_local
      if (iblock > 0) then
        mat(row0 + 1:row0 + niface, row0 - 1) = -tl_local(:, 1)
        mat(row0 + 1:row0 + niface, row0) = -tl_local(:, 2)
      end if
      if (iblock < npy_grid - 1) then
        mat(row0 + 1:row0 + niface, row0 + 5) = -tr_local(:, 1)
        mat(row0 + 1:row0 + niface, row0 + 6) = -tr_local(:, 2)
      end if
    end do

    sol = rhs
    call ys_solve_dense_complex(mat, sol)
    row0 = 4*(ipy)
    u_left = (0.0d0, 0.0d0)
    u_right = (0.0d0, 0.0d0)
    if (ipy > 0) u_left = sol(row0 - 1:row0)
    if (ipy < npy_grid - 1) u_right = sol(row0 + 5:row0 + 6)

    rhs_local = rhs_base - left_couple(:, 1)*u_left(1) - left_couple(:, 2)*u_left(2) - &
                right_couple(:, 1)*u_right(1) - right_couple(:, 2)*u_right(2)
    call ys_solve_factored_banded_dense(factored_a, rhs_local, 2_C_INT)

    local_block = (0.0d0, 0.0d0)
    do iblock = row_start, row_end
      local_block(iblock - (yl0 - 2) + 1) = rhs_local(iblock - row_start + 1)
    end do

    if (has_lower_boundary) then
      local_block(2) = (lower_rhs0 - sum(lower_eq0(0:2)*[local_block(3), local_block(4), local_block(5)]))/lower_eq0(-1)
      local_block(1) = (lower_rhsm1 - sum(eqm1(-1:2)*[local_block(2), local_block(3), local_block(4), local_block(5)]))/eqm1(-2)
    end if
    if (has_upper_boundary) then
      local_block(ny - (yl0 - 2) + 1) = (upper_rhsn - sum(upper_eqn(-2:0)*[ &
                                                          local_block(ny - 3 - (yl0 - 2) + 1), &
                                                          local_block(ny - 2 - (yl0 - 2) + 1), &
                                                          local_block(ny - 1 - (yl0 - 2) + 1)]))/upper_eqn(1)
      local_block(ny + 1 - (yl0 - 2) + 1) = (upper_rhsnp1 - sum(eqnp1(-2:1)*[ &
                                                                local_block(ny - 3 - (yl0 - 2) + 1), &
                                                                local_block(ny - 2 - (yl0 - 2) + 1), &
                                                                local_block(ny - 1 - (yl0 - 2) + 1), &
                                                                local_block(ny - (yl0 - 2) + 1)]))/eqnp1(2)
    end if
    deallocate (gathered, mat, rhs, sol)
    deallocate (particular, left_resp, right_resp, factored_a, rhs_local, rhs_base, left_couple, right_couple)
#else
    call ys_solve_ghost_block_local(line, a, eqm1, eq0, eqn, eqnp1, ny, local_block)
#endif
  end subroutine ys_solve_ghost_block_reduced

  subroutine ys_solve_ghost_block_local(line, a, eqm1, eq0, eqn, eqnp1, ny, local_block)
    integer(C_INT), intent(in) :: ny
    complex(C_DOUBLE_COMPLEX), intent(in) :: line(-1:ny + 1)
    real(C_DOUBLE), intent(in) :: a(1:ny + 1, -2:2)
    real(C_DOUBLE), intent(in) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: local_block(ylB)
    complex(C_DOUBLE_COMPLEX) :: work(-1:ny + 1)
    real(C_DOUBLE) :: amat(1:ny + 1, -2:2), beqm1(-2:2), beq0(-2:2), beqn(-2:2), beqnp1(-2:2)

    work = line
    amat = a
    beqm1 = eqm1
    beq0 = eq0
    beqn = eqn
    beqnp1 = eqnp1
    call ys_solve_ghost_system(work, amat, beqm1, beq0, beqn, beqnp1, ny)
    local_block = work(yl0 - 2:ylN - 2)
  end subroutine ys_solve_ghost_block_local

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

  subroutine ys_factor_banded_dense(mat, half_bw)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: mat(:, :)
    integer(C_INT), intent(in) :: half_bw
    integer(C_INT) :: n, i, j, k, jmax, kmax
    complex(C_DOUBLE_COMPLEX) :: factor

    n = int(size(mat, 1), C_INT)
    do i = 1, n
      jmax = min(n, i + half_bw)
      do j = i + 1, jmax
        factor = mat(j, i)/mat(i, i)
        mat(j, i) = factor
        kmax = min(n, min(i + half_bw, j + half_bw))
        do k = i + 1, kmax
          mat(j, k) = mat(j, k) - factor*mat(i, k)
        end do
      end do
    end do
  end subroutine ys_factor_banded_dense

  subroutine ys_solve_factored_banded_dense(mat, rhs, half_bw)
    complex(C_DOUBLE_COMPLEX), intent(in) :: mat(:, :)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    integer(C_INT), intent(in) :: half_bw
    integer(C_INT) :: n, i, j

    n = int(size(rhs), C_INT)
    do i = 1, n
      do j = max(1_C_INT, i - half_bw), i - 1
        rhs(i) = rhs(i) - mat(i, j)*rhs(j)
      end do
    end do

    do i = n, 1, -1
      do j = i + 1, min(n, i + half_bw)
        rhs(i) = rhs(i) - mat(i, j)*rhs(j)
      end do
      rhs(i) = rhs(i)/mat(i, i)
    end do
  end subroutine ys_solve_factored_banded_dense

end module y_line_solvers
