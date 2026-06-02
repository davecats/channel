#include "header.h"

module y_line_solvers

  use, intrinsic :: iso_c_binding
  use mpi_transpose, only: ny0, nx0, nxN, nxB, yl0, ylN, ylB, zpy0, zpyB, &
                           transpose_xz_to_y_pencil, transpose_y_pencil_to_xz, allgather_y_blocks_to_xz_full

  implicit none
  private

  public :: ys_lu5decomp, ys_leftlu5div
  public :: ys_solve_compact_derivative, ys_solve_compact_system, ys_solve_ghost_system
  public :: ys_solve_ghost_field_with_y_pencil

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

end module y_line_solvers
