#include "header.h"

module y_line_solvers

  use, intrinsic :: iso_c_binding
  use mpi_transpose, only: ny0, nx0, nxN, nxB, yl0, ylN, ylB, zpy0, zpyB, npy_grid, ierr, ipy, MPI_COMM_Y, &
                           transpose_xz_to_y_pencil, transpose_y_pencil_to_xz, allgather_y_blocks_to_xz_full
#ifdef HAVE_MPI
  use mpi_f08
#endif

  implicit none
  private

  public :: ys_lu5decomp, ys_leftlu5div
  public :: ys_solve_compact_derivative, ys_solve_compact_system, ys_solve_ghost_system
  public :: ys_solve_ghost_field, ys_solve_ghost_field_with_y_pencil, ys_solve_ghost_field_reduced

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

    if (npy_grid > 1 .and. ylB >= 4) then
      call ys_solve_ghost_field_reduced(src0, src1, dst, ny, nz, build_line)
    else
      call ys_solve_ghost_field_with_y_pencil(src0, src1, dst, ny, nz, build_line)
    end if
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

  subroutine ys_solve_ghost_field_with_y_pencil(src0, src1, dst, ny, nz, build_line)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0(:, :, :), src1(:, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dst(:, :, :)
    procedure(ys_build_ghost_line) :: build_line
    complex(C_DOUBLE_COMPLEX), allocatable :: src0_xz(:, :, :), src1_xz(:, :, :), dst_xz(:, :, :), dst_xz_full(:, :, :)
    complex(C_DOUBLE_COMPLEX), allocatable :: src0_y(:, :, :), src1_y(:, :, :), dst_y(:, :, :)
    complex(C_DOUBLE_COMPLEX) :: line(-1:ny + 1)
    real(C_DOUBLE) :: a(1:ny + 1, -2:2), eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    integer(C_INT) :: ix, iz, ix_local, iz_local, ix_global, iz_global

    allocate (src0_xz(ylB, 2*nz + 1, nxB), src1_xz(ylB, 2*nz + 1, nxB), dst_xz(ylB, 2*nz + 1, nxB))
    allocate (dst_xz_full(ny + 3, 2*nz + 1, nxB))
    allocate (src0_y(ny + 3, zpyB, nxB), src1_y(ny + 3, zpyB, nxB), dst_y(ny + 3, zpyB, nxB))

    do ix = nx0, nxN
      do iz = -nz, nz
        src0_xz(:, iz + nz + 1, ix - nx0 + 1) = src0(yl0:ylN, iz + nz + 1, ix - nx0 + 1)
        src1_xz(:, iz + nz + 1, ix - nx0 + 1) = src1(yl0:ylN, iz + nz + 1, ix - nx0 + 1)
      end do
    end do

    call transpose_xz_to_y_pencil(src0_xz, src0_y)
    call transpose_xz_to_y_pencil(src1_xz, src1_y)

    do ix_local = 1, nxB
      ix_global = nx0 + ix_local - 1
      do iz_local = 1, zpyB
        iz_global = zpy0 + iz_local - 1 - (nz + 1)
        call build_line(ix_global, iz_global, src0_y(:, iz_local, ix_local), src1_y(:, iz_local, ix_local), &
                        line, a, eqm1, eq0, eqn, eqnp1, ny)
        call ys_solve_ghost_system(line, a, eqm1, eq0, eqn, eqnp1, ny)
        dst_y(:, iz_local, ix_local) = line
      end do
    end do

    call transpose_y_pencil_to_xz(dst_y, dst_xz)
    call allgather_y_blocks_to_xz_full(dst_xz, dst_xz_full)

    do ix = nx0, nxN
      do iz = -nz, nz
        dst(:, iz + nz + 1, ix - nx0 + 1) = dst_xz_full(:, iz + nz + 1, ix - nx0 + 1)
      end do
    end do

    deallocate (src0_xz, src1_xz, dst_xz, dst_xz_full, src0_y, src1_y, dst_y)
  end subroutine ys_solve_ghost_field_with_y_pencil

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
    complex(C_DOUBLE_COMPLEX) :: factored_a(ylB, ylB)
    complex(C_DOUBLE_COMPLEX) :: rhs_local(0:ylB - 1), left_couple(0:ylB - 1, 2), right_couple(0:ylB - 1, 2)
    complex(C_DOUBLE_COMPLEX), allocatable :: packed(:), gathered(:), mat(:, :), rhs(:), sol(:)
    integer(C_INT) :: iblock, row0

    if (npy_grid == 1 .or. ylB < 4) then
      call ys_solve_ghost_block_local(line, a, eqm1, eq0, eqn, eqnp1, ny, local_block)
      return
    end if

    call ys_build_local_reduction(line, a, eqm1, eq0, eqn, eqnp1, ny, c_local, tl_local, tr_local, &
                                  factored_a, rhs_local, left_couple, right_couple)

#ifdef HAVE_MPI
    allocate (packed(20), gathered(20*npy_grid))
    call ys_pack_reduction(c_local, tl_local, tr_local, packed)
    call MPI_Allgather(packed, int(size(packed), kind=4), MPI_DOUBLE_COMPLEX, &
                       gathered, int(size(packed), kind=4), MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr)

    allocate (mat(niface*npy_grid, niface*npy_grid), rhs(niface*npy_grid), sol(niface*npy_grid))
    mat = (0.0d0, 0.0d0)
    rhs = (0.0d0, 0.0d0)
    do iblock = 0, npy_grid - 1
      row0 = iblock*niface
      call ys_unpack_reduction(gathered(20*iblock + 1:20*(iblock + 1)), c_local, tl_local, tr_local)
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
    call ys_reconstruct_local_block(factored_a, rhs_local, left_couple, right_couple, u_left, u_right, local_block)
    deallocate (packed, gathered, mat, rhs, sol)
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

  subroutine ys_build_local_reduction(line, a, eqm1, eq0, eqn, eqnp1, ny, c_local, tl_local, tr_local, factored_a, &
                                      rhs_local, left_couple, right_couple)
    integer(C_INT), intent(in) :: ny
    complex(C_DOUBLE_COMPLEX), intent(in) :: line(-1:ny + 1)
    real(C_DOUBLE), intent(in) :: a(1:ny + 1, -2:2)
    real(C_DOUBLE), intent(in) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: c_local(4), tl_local(4, 2), tr_local(4, 2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: factored_a(ylB, ylB)
    complex(C_DOUBLE_COMPLEX), intent(out) :: rhs_local(0:ylB - 1), left_couple(0:ylB - 1, 2), right_couple(0:ylB - 1, 2)
    complex(C_DOUBLE_COMPLEX) :: rhs_work(ylB)
    complex(C_DOUBLE_COMPLEX) :: particular(0:ylB - 1), left_resp(0:ylB - 1, 2), right_resp(0:ylB - 1, 2)
    complex(C_DOUBLE_COMPLEX) :: local_a(ylB, ylB)
    integer(C_INT) :: row, idx, cols(5), k, col
    complex(C_DOUBLE_COMPLEX) :: row_left(2), row_right(2)
    real(C_DOUBLE) :: coeffs(5)
    integer(C_INT), parameter :: local_half_bw = 4

    if (ylB < 4) error stop "ys_build_local_reduction requires at least 4 local y rows"

    local_a = (0.0d0, 0.0d0)
    rhs_local = (0.0d0, 0.0d0)
    left_couple = (0.0d0, 0.0d0)
    right_couple = (0.0d0, 0.0d0)

    do row = yl0 - 2, ylN - 2
      idx = row - (yl0 - 2)
      call ys_get_ghost_row(row, ny, a, eqm1, eq0, eqn, eqnp1, cols, coeffs)
      rhs_local(idx) = line(row)
      call ys_scatter_row(cols, coeffs, row_left, row_right, ny)
      left_couple(idx, :) = row_left
      right_couple(idx, :) = row_right
      do k = 1, 5
        col = cols(k)
        if (col >= yl0 - 2 .and. col <= ylN - 2) then
          local_a(idx + 1, col - (yl0 - 2) + 1) = cmplx(coeffs(k), 0.0d0, kind=C_DOUBLE)
        end if
      end do
    end do

    factored_a = local_a
    call ys_factor_banded_dense(factored_a, local_half_bw)

    rhs_work = rhs_local
    call ys_solve_factored_banded_dense(factored_a, rhs_work, local_half_bw)
    particular = rhs_work

    rhs_work = -left_couple(:, 1)
    call ys_solve_factored_banded_dense(factored_a, rhs_work, local_half_bw)
    left_resp(:, 1) = rhs_work

    rhs_work = -left_couple(:, 2)
    call ys_solve_factored_banded_dense(factored_a, rhs_work, local_half_bw)
    left_resp(:, 2) = rhs_work

    rhs_work = -right_couple(:, 1)
    call ys_solve_factored_banded_dense(factored_a, rhs_work, local_half_bw)
    right_resp(:, 1) = rhs_work

    rhs_work = -right_couple(:, 2)
    call ys_solve_factored_banded_dense(factored_a, rhs_work, local_half_bw)
    right_resp(:, 2) = rhs_work

    c_local = [particular(0), particular(1), particular(ylB - 2), particular(ylB - 1)]
    tl_local(:, 1) = [left_resp(0, 1), left_resp(1, 1), left_resp(ylB - 2, 1), left_resp(ylB - 1, 1)]
    tl_local(:, 2) = [left_resp(0, 2), left_resp(1, 2), left_resp(ylB - 2, 2), left_resp(ylB - 1, 2)]
    tr_local(:, 1) = [right_resp(0, 1), right_resp(1, 1), right_resp(ylB - 2, 1), right_resp(ylB - 1, 1)]
    tr_local(:, 2) = [right_resp(0, 2), right_resp(1, 2), right_resp(ylB - 2, 2), right_resp(ylB - 1, 2)]
  end subroutine ys_build_local_reduction

  subroutine ys_reconstruct_local_block(factored_a, rhs_local, left_couple, right_couple, u_left, u_right, local_block)
    complex(C_DOUBLE_COMPLEX), intent(in) :: factored_a(ylB, ylB)
    complex(C_DOUBLE_COMPLEX), intent(in) :: rhs_local(0:ylB - 1), left_couple(0:ylB - 1, 2), right_couple(0:ylB - 1, 2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: u_left(2), u_right(2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: local_block(ylB)
    complex(C_DOUBLE_COMPLEX) :: rhs_work(ylB)
    integer(C_INT), parameter :: local_half_bw = 4

    rhs_work = rhs_local - left_couple(:, 1)*u_left(1) - left_couple(:, 2)*u_left(2) - &
               right_couple(:, 1)*u_right(1) - right_couple(:, 2)*u_right(2)
    call ys_solve_factored_banded_dense(factored_a, rhs_work, local_half_bw)
    local_block = rhs_work
  end subroutine ys_reconstruct_local_block

  subroutine ys_get_ghost_row(row, ny, a, eqm1, eq0, eqn, eqnp1, cols, coeffs)
    integer(C_INT), intent(in) :: row, ny
    real(C_DOUBLE), intent(in) :: a(1:ny + 1, -2:2)
    real(C_DOUBLE), intent(in) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    integer(C_INT), intent(out) :: cols(5)
    real(C_DOUBLE), intent(out) :: coeffs(5)
    integer(C_INT) :: i

    if (row == -1) then
      cols = [(-1_C_INT + i, i=0_C_INT, 4_C_INT)]
      coeffs = eqm1(-2:2)
    else if (row == 0) then
      cols = [(-1_C_INT + i, i=0_C_INT, 4_C_INT)]
      coeffs = eq0(-2:2)
    else if (row == ny) then
      cols = [(ny - 3_C_INT + i, i=0_C_INT, 4_C_INT)]
      coeffs = eqn(-2:2)
    else if (row == ny + 1) then
      cols = [(ny - 3_C_INT + i, i=0_C_INT, 4_C_INT)]
      coeffs = eqnp1(-2:2)
    else
      cols = [(row - 2_C_INT + i, i=0_C_INT, 4_C_INT)]
      coeffs = a(row, -2:2)
    end if
  end subroutine ys_get_ghost_row

  subroutine ys_scatter_row(cols, coeffs, left_row, right_row, ny)
    integer(C_INT), intent(in) :: cols(5), ny
    real(C_DOUBLE), intent(in) :: coeffs(5)
    complex(C_DOUBLE_COMPLEX), intent(out) :: left_row(2), right_row(2)
    integer(C_INT) :: k, col

    left_row = (0.0d0, 0.0d0)
    right_row = (0.0d0, 0.0d0)
    do k = 1, 5
      col = cols(k)
      if (col == yl0 - 2 - 2 .and. col >= -1) then
        left_row(1) = cmplx(coeffs(k), 0.0d0, kind=C_DOUBLE)
      else if (col == yl0 - 2 - 1 .and. col >= -1) then
        left_row(2) = cmplx(coeffs(k), 0.0d0, kind=C_DOUBLE)
      else if (col == ylN - 2 + 1 .and. col <= ny + 1) then
        right_row(1) = cmplx(coeffs(k), 0.0d0, kind=C_DOUBLE)
      else if (col == ylN - 2 + 2 .and. col <= ny + 1) then
        right_row(2) = cmplx(coeffs(k), 0.0d0, kind=C_DOUBLE)
      end if
    end do
  end subroutine ys_scatter_row

  subroutine ys_pack_reduction(c_local, tl_local, tr_local, packed)
    complex(C_DOUBLE_COMPLEX), intent(in) :: c_local(4), tl_local(4, 2), tr_local(4, 2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: packed(20)

    packed(1:4) = c_local
    packed(5:12) = reshape(tl_local, [8])
    packed(13:20) = reshape(tr_local, [8])
  end subroutine ys_pack_reduction

  subroutine ys_unpack_reduction(packed, c_local, tl_local, tr_local)
    complex(C_DOUBLE_COMPLEX), intent(in) :: packed(20)
    complex(C_DOUBLE_COMPLEX), intent(out) :: c_local(4), tl_local(4, 2), tr_local(4, 2)

    c_local = packed(1:4)
    tl_local = reshape(packed(5:12), [4, 2])
    tr_local = reshape(packed(13:20), [4, 2])
  end subroutine ys_unpack_reduction

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
