#include "header.h"

module y_line_solvers

  use, intrinsic :: iso_c_binding
  use mpi_transpose, only: ny0, nyN, nx0, nxN, nxB, npy_grid, ierr, ipy, MPI_COMM_Y
#ifdef HAVE_MPI
  use mpi_f08
#endif

  implicit none
  private

  public :: ys_lu5decomp, ys_leftlu5div
  public :: ys_solve_compact_derivative, ys_solve_compact_system, ys_solve_ghost_system
  public :: ys_prepare_ghost_field_workspace, ys_solve_ghost_field, ys_solve_ghost_field_reduced, ys_fill_ghost_padded_field
  public :: ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store

  integer(C_INT), save :: ys_workspace_ny = -1
  integer(C_INT), save :: ys_workspace_nz = -1
  integer(C_INT), save :: ys_workspace_nx = -1
  integer(C_INT), save :: ys_workspace_nlines = 0
  integer(C_INT), save :: ys_workspace_active_n = 0
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_rhs_store(:, :)
  real(C_DOUBLE), allocatable, save :: ys_matrix_store(:, :, :)
  real(C_DOUBLE), allocatable, save :: ys_eqm1_store(:, :), ys_eq0_store(:, :), ys_eqn_store(:, :), ys_eqnp1_store(:, :)
  real(C_DOUBLE), allocatable, save :: ys_factored_store(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_store(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_packed_send(:, :), ys_left_u(:, :), ys_right_u(:, :)

contains

  subroutine ys_prepare_ghost_field_workspace(ny, nz, nx_lines)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, nx_lines
    integer(C_INT) :: row_start, row_end, active_n, nlines

    ! The reduced interface solve keeps two boundary-adjacent unknowns per side,
    ! so each rank must own at least four physical/ghost rows.
    if (npy_grid > 1 .and. nyN - ny0 + 1 < 4) error stop "ys_solve_ghost_field requires at least four y rows per rank"
    row_start = ny0
    row_end = nyN
    active_n = row_end - row_start + 1
    if (active_n < 4) error stop "ys_solve_ghost_field requires at least four active y rows"

    nlines = nx_lines*(2*nz + 1)
    if (allocated(ys_rhs_store)) then
      if (ys_workspace_ny /= ny .or. ys_workspace_nz /= nz .or. ys_workspace_nx /= nx_lines .or. &
          ys_workspace_active_n /= active_n) then
        deallocate (ys_rhs_store, ys_matrix_store, ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store)
        deallocate (ys_factored_store, ys_reduced_store, ys_packed_send, ys_left_u, ys_right_u)
      end if
    end if

    if (.not. allocated(ys_rhs_store)) then
      allocate (ys_rhs_store(-1:ny + 1, nlines), ys_matrix_store(1:ny + 1, -2:2, nlines))
      allocate (ys_eqm1_store(-2:2, nlines), ys_eq0_store(-2:2, nlines), ys_eqn_store(-2:2, nlines), ys_eqnp1_store(-2:2, nlines))
      allocate (ys_factored_store(0:active_n - 1, -2:2, nlines), ys_reduced_store(0:active_n - 1, 5, nlines))
      allocate (ys_packed_send(20, nlines), ys_left_u(2, nlines), ys_right_u(2, nlines))
    end if

    ys_workspace_ny = ny
    ys_workspace_nz = nz
    ys_workspace_nx = nx_lines
    ys_workspace_nlines = nlines
    ys_workspace_active_n = active_n

    ys_rhs_store = (0.0d0, 0.0d0)
    ys_matrix_store = 0.0d0
    ys_eqm1_store = 0.0d0
    ys_eq0_store = 0.0d0
    ys_eqn_store = 0.0d0
    ys_eqnp1_store = 0.0d0
    ys_factored_store = 0.0d0
    ys_reduced_store = (0.0d0, 0.0d0)
    ys_packed_send = (0.0d0, 0.0d0)
    ys_left_u = (0.0d0, 0.0d0)
    ys_right_u = (0.0d0, 0.0d0)
  end subroutine ys_prepare_ghost_field_workspace

  subroutine ys_solve_ghost_field(dst, ny, nz)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(out) :: dst(:, :, :)

    call ys_solve_ghost_field_reduced(dst, ny, nz)
  end subroutine ys_solve_ghost_field

  subroutine ys_fill_ghost_padded_field(dst, ny, nz)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX) :: lower_pair(2), upper_pair(2)
    complex(C_DOUBLE_COMPLEX) :: lower_rhs0, upper_rhsn
    real(C_DOUBLE) :: lower_eq0(-2:2), upper_eqn(-2:2)
    integer(C_INT) :: ix, iz, iline, nlines_z

   if (.not. allocated(ys_rhs_store)) error stop "ys_prepare_ghost_field_workspace must be called before ys_fill_ghost_padded_field"

    nlines_z = 2*nz + 1

    do ix = nx0, nxN
      do iz = -nz, nz
        iline = (ix - nx0)*nlines_z + (iz + nz + 1)
        lower_pair = ys_left_u(:, iline)
        upper_pair = ys_right_u(:, iline)

        if (ny0 == 1) then
          lower_rhs0 = ys_rhs_store(0, iline) - ys_rhs_store(-1, iline)*ys_eq0_store(-2, iline)/ys_eqm1_store(-2, iline)
          lower_eq0 = ys_eq0_store(:, iline) - ys_eqm1_store(:, iline)*ys_eq0_store(-2, iline)/ys_eqm1_store(-2, iline)
          lower_eq0(-2) = 0.0d0

          dst(0, iz, ix) = (lower_rhs0 - sum(lower_eq0(0:2)*dst(1:3, iz, ix)))/lower_eq0(-1)
          dst(-1, iz, ix) = (ys_rhs_store(-1, iline) - sum(ys_eqm1_store(-1:2, iline)*dst(0:3, iz, ix)))/ys_eqm1_store(-2, iline)
        else
          dst(ny0 - 2, iz, ix) = lower_pair(1)
          dst(ny0 - 1, iz, ix) = lower_pair(2)
        end if

        if (nyN == ny - 1) then
          upper_rhsn = ys_rhs_store(ny, iline) - ys_rhs_store(ny + 1, iline)*ys_eqn_store(2, iline)/ys_eqnp1_store(2, iline)
          upper_eqn = ys_eqn_store(:, iline) - ys_eqnp1_store(:, iline)*ys_eqn_store(2, iline)/ys_eqnp1_store(2, iline)
          upper_eqn(2) = 0.0d0

          dst(ny, iz, ix) = (upper_rhsn - sum(upper_eqn(-2:0)*dst(ny - 3:ny - 1, iz, ix)))/upper_eqn(1)
          dst(ny + 1, iz, ix) = (ys_rhs_store(ny + 1, iline) - sum(ys_eqnp1_store(-2:1, iline)*dst(ny - 3:ny, iz, ix)))/ys_eqnp1_store(2, iline)
        else
          dst(nyN + 1, iz, ix) = upper_pair(1)
          dst(nyN + 2, iz, ix) = upper_pair(2)
        end if
      end do
    end do
  end subroutine ys_fill_ghost_padded_field

  !$omp begin declare target
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
  !$omp end declare target

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

  subroutine ys_solve_ghost_field_reduced(dst, ny, nz)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(out) :: dst(:, :, :)
    complex(C_DOUBLE_COMPLEX) :: local_block(size(dst, 1))
    complex(C_DOUBLE_COMPLEX) :: u_left(2), u_right(2)
    complex(C_DOUBLE_COMPLEX) :: lower_rhs0, upper_rhsn
    real(C_DOUBLE) :: row_coeffs(-2:2), lower_eq0(-2:2), upper_eqn(-2:2)
    integer(C_INT) :: ix, iz, iline, nlines, nlines_z, row_start, row_end, active_n, row, idx, col, global_col
    logical :: has_lower_boundary, has_upper_boundary

    row_start = ny0
    row_end = nyN
    active_n = row_end - row_start + 1
    has_lower_boundary = (row_start == 1)
    has_upper_boundary = (row_end == ny - 1)
    if (.not. allocated(ys_rhs_store)) error stop "ys_prepare_ghost_field_workspace must be called before ys_solve_ghost_field"
    if (ys_workspace_ny /= ny .or. ys_workspace_nz /= nz .or. ys_workspace_nx /= size(dst, 3) .or. &
        ys_workspace_active_n /= size(dst, 1)) then
      error stop "ys_solve_ghost_field workspace does not match requested local solve dimensions"
    end if

    nlines_z = 2*nz + 1
    nlines = ys_workspace_nlines

    do iline = 1, nlines
      ys_factored_store(:, :, iline) = 0.0d0
      ys_reduced_store(:, :, iline) = (0.0d0, 0.0d0)

      lower_rhs0 = (0.0d0, 0.0d0)
      lower_eq0 = 0.0d0
      if (has_lower_boundary) then
        lower_rhs0 = ys_rhs_store(0, iline) - ys_rhs_store(-1, iline)*ys_eq0_store(-2, iline)/ys_eqm1_store(-2, iline)
        lower_eq0 = ys_eq0_store(:, iline) - ys_eqm1_store(:, iline)*ys_eq0_store(-2, iline)/ys_eqm1_store(-2, iline)
        lower_eq0(-2) = 0.0d0
      end if
      upper_rhsn = (0.0d0, 0.0d0)
      upper_eqn = 0.0d0
      if (has_upper_boundary) then
        upper_rhsn = ys_rhs_store(ny, iline) - ys_rhs_store(ny + 1, iline)*ys_eqn_store(2, iline)/ys_eqnp1_store(2, iline)
        upper_eqn = ys_eqn_store(:, iline) - ys_eqnp1_store(:, iline)*ys_eqn_store(2, iline)/ys_eqnp1_store(2, iline)
        upper_eqn(2) = 0.0d0
      end if

      do row = row_start, row_end
        ys_reduced_store(row - row_start, 1, iline) = ys_rhs_store(row, iline)
      end do

      if (has_lower_boundary) then
        ys_reduced_store(0, 1, iline) = ys_reduced_store(0, 1, iline) - ys_rhs_store(-1, iline)*ys_matrix_store(1, -2, iline)/ys_eqm1_store(-2, iline)
       row_coeffs = ys_matrix_store(1, -2:2, iline) - ys_eqm1_store(:, iline)*ys_matrix_store(1, -2, iline)/ys_eqm1_store(-2, iline)
        row_coeffs(-2) = 0.0d0
        ys_reduced_store(0, 1, iline) = ys_reduced_store(0, 1, iline) - lower_rhs0*row_coeffs(-1)/lower_eq0(-1)
        row_coeffs = row_coeffs - lower_eq0*row_coeffs(-1)/lower_eq0(-1)
        row_coeffs(-1) = 0.0d0
        ys_factored_store(0, 0:2, iline) = row_coeffs(0:2)

        ys_reduced_store(1, 1, iline) = ys_reduced_store(1, 1, iline) - lower_rhs0*ys_matrix_store(2, -2, iline)/lower_eq0(-1)
        row_coeffs = ys_matrix_store(2, -2:2, iline)
        row_coeffs(-2:1) = row_coeffs(-2:1) - lower_eq0(-1:2)*ys_matrix_store(2, -2, iline)/lower_eq0(-1)
        row_coeffs(-2) = 0.0d0
        ys_factored_store(1, -1:2, iline) = row_coeffs(-1:2)
      end if

      if (has_upper_boundary) then
        ys_reduced_store(active_n - 1, 1, iline) = ys_reduced_store(active_n - 1, 1, iline) - &
                                              ys_rhs_store(ny + 1, iline)*ys_matrix_store(ny - 1, 2, iline)/ys_eqnp1_store(2, iline)
        row_coeffs = ys_matrix_store(ny - 1, -2:2, iline) - ys_eqnp1_store(:, iline)*ys_matrix_store(ny - 1, 2, iline)/ys_eqnp1_store(2, iline)
        row_coeffs(2) = 0.0d0
        ys_reduced_store(active_n - 1, 1, iline) = ys_reduced_store(active_n - 1, 1, iline) - upper_rhsn*row_coeffs(1)/upper_eqn(1)
        row_coeffs = row_coeffs - upper_eqn*row_coeffs(1)/upper_eqn(1)
        row_coeffs(1) = 0.0d0
        ys_factored_store(active_n - 1, -2:0, iline) = row_coeffs(-2:0)

        ys_reduced_store(active_n - 2, 1, iline) = ys_reduced_store(active_n - 2, 1, iline) - upper_rhsn*ys_matrix_store(ny - 2, 2, iline)/upper_eqn(1)
        row_coeffs = ys_matrix_store(ny - 2, -2:2, iline)
        row_coeffs(-1:2) = row_coeffs(-1:2) - upper_eqn(-2:1)*ys_matrix_store(ny - 2, 2, iline)/upper_eqn(1)
        row_coeffs(2) = 0.0d0
        ys_factored_store(active_n - 2, -2:1, iline) = row_coeffs(-2:1)
      end if

      do row = row_start, row_end
        idx = row - row_start
        if ((has_lower_boundary .and. row <= 2) .or. (has_upper_boundary .and. row >= ny - 2)) cycle
        row_coeffs = ys_matrix_store(row, -2:2, iline)
        do col = -2, 2
          global_col = row + col
          if (global_col < row_start) then
            if (global_col == row_start - 2) ys_reduced_store(idx, 2, iline) = -cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
            if (global_col == row_start - 1) ys_reduced_store(idx, 3, iline) = -cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
          else if (global_col > row_end) then
            if (global_col == row_end + 1) ys_reduced_store(idx, 4, iline) = -cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
            if (global_col == row_end + 2) ys_reduced_store(idx, 5, iline) = -cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
          else
            ys_factored_store(idx, col, iline) = row_coeffs(col)
          end if
        end do
      end do

      call ys_factor_penta(ys_factored_store(:, :, iline))
      call ys_solve_factored_penta_multi(ys_reduced_store(:, :, iline), ys_factored_store(:, :, iline))

      ys_packed_send(1:4, iline) = [ys_reduced_store(0, 1, iline), ys_reduced_store(1, 1, iline), &
                                    ys_reduced_store(active_n - 2, 1, iline), ys_reduced_store(active_n - 1, 1, iline)]
      ys_packed_send(5:12, iline) = [ys_reduced_store(0, 2, iline), ys_reduced_store(1, 2, iline), &
                                     ys_reduced_store(active_n - 2, 2, iline), ys_reduced_store(active_n - 1, 2, iline), &
                                     ys_reduced_store(0, 3, iline), ys_reduced_store(1, 3, iline), &
                                     ys_reduced_store(active_n - 2, 3, iline), ys_reduced_store(active_n - 1, 3, iline)]
      ys_packed_send(13:20, iline) = [ys_reduced_store(0, 4, iline), ys_reduced_store(1, 4, iline), &
                                      ys_reduced_store(active_n - 2, 4, iline), ys_reduced_store(active_n - 1, 4, iline), &
                                      ys_reduced_store(0, 5, iline), ys_reduced_store(1, 5, iline), &
                                      ys_reduced_store(active_n - 2, 5, iline), ys_reduced_store(active_n - 1, 5, iline)]
    end do

    call ys_solve_reduced_interfaces(ys_packed_send, ys_left_u, ys_right_u)

    do ix = nx0, nxN
      do iz = -nz, nz
        iline = (ix - nx0)*nlines_z + (iz + nz + 1)
        u_left = ys_left_u(:, iline)
        u_right = ys_right_u(:, iline)

        ys_reduced_store(:, 1, iline) = ys_reduced_store(:, 1, iline) + ys_reduced_store(:, 2, iline)*u_left(1) + &
                                        ys_reduced_store(:, 3, iline)*u_left(2) + ys_reduced_store(:, 4, iline)*u_right(1) + &
                                        ys_reduced_store(:, 5, iline)*u_right(2)

        local_block = (0.0d0, 0.0d0)
        do row = row_start, row_end
          local_block(row - row_start + 1) = ys_reduced_store(row - row_start, 1, iline)
        end do

        dst(:, iz + nz + 1, ix - nx0 + 1) = local_block
      end do
    end do

  end subroutine ys_solve_ghost_field_reduced

  subroutine ys_solve_reduced_interfaces(packed_send, left_u, right_u)
    complex(C_DOUBLE_COMPLEX), intent(in) :: packed_send(:, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: left_u(:, :), right_u(:, :)
    complex(C_DOUBLE_COMPLEX), allocatable :: packed_recv(:, :, :), band_a(:, :), rhs(:)
    complex(C_DOUBLE_COMPLEX) :: packed_remote(20)
    integer(C_INT), parameter :: bw = 5
    integer(C_INT) :: nlines, niface, iline, iblock, row0

    nlines = size(packed_send, 2)
    niface = 4*npy_grid
    allocate (packed_recv(20, nlines, npy_grid), band_a(1:niface, 1:2*bw + 1), rhs(niface))

#ifdef HAVE_MPI
    call MPI_Allgather(packed_send, 20*nlines, MPI_DOUBLE_COMPLEX, packed_recv, 20*nlines, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr)
#else
    packed_recv(:, :, 1) = packed_send
#endif

    do iline = 1, nlines
      band_a = (0.0d0, 0.0d0)
      rhs = (0.0d0, 0.0d0)
      do iblock = 0, npy_grid - 1
        row0 = 4*iblock
        packed_remote = packed_recv(:, iline, iblock + 1)
        band_a(row0 + 1, bw + 1) = (1.0d0, 0.0d0)
        band_a(row0 + 2, bw + 1) = (1.0d0, 0.0d0)
        band_a(row0 + 3, bw + 1) = (1.0d0, 0.0d0)
        band_a(row0 + 4, bw + 1) = (1.0d0, 0.0d0)
        rhs(row0 + 1:row0 + 4) = packed_remote(1:4)
        if (iblock > 0) then
          band_a(row0 + 1, bw - 1) = -packed_remote(5)
          band_a(row0 + 2, bw - 2) = -packed_remote(6)
          band_a(row0 + 3, bw - 3) = -packed_remote(7)
          band_a(row0 + 4, bw - 4) = -packed_remote(8)
          band_a(row0 + 1, bw) = -packed_remote(9)
          band_a(row0 + 2, bw - 1) = -packed_remote(10)
          band_a(row0 + 3, bw - 2) = -packed_remote(11)
          band_a(row0 + 4, bw - 3) = -packed_remote(12)
        end if
        if (iblock < npy_grid - 1) then
          band_a(row0 + 1, bw + 5) = -packed_remote(13)
          band_a(row0 + 2, bw + 4) = -packed_remote(14)
          band_a(row0 + 3, bw + 3) = -packed_remote(15)
          band_a(row0 + 4, bw + 2) = -packed_remote(16)
          band_a(row0 + 1, bw + 6) = -packed_remote(17)
          band_a(row0 + 2, bw + 5) = -packed_remote(18)
          band_a(row0 + 3, bw + 4) = -packed_remote(19)
          band_a(row0 + 4, bw + 3) = -packed_remote(20)
        end if
      end do

      call ys_factor_banded_complex(band_a)
      call ys_solve_factored_banded_complex(rhs, band_a)
      left_u(:, iline) = (0.0d0, 0.0d0)
      right_u(:, iline) = (0.0d0, 0.0d0)
      row0 = 4*ipy
      if (ipy > 0) then
        left_u(:, iline) = rhs(row0 - 1:row0)
      end if
      if (ipy < npy_grid - 1) then
        right_u(:, iline) = rhs(row0 + 5:row0 + 6)
      end if
    end do

    deallocate (packed_recv, band_a, rhs)
  end subroutine ys_solve_reduced_interfaces

  !$omp begin declare target
  subroutine ys_factor_banded_complex(a)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: a(:, :)
    integer(C_INT), parameter :: bw = 5
    integer(C_INT) :: n, i, j, t
    complex(C_DOUBLE_COMPLEX) :: piv, factor

    n = size(a, 1)
    do i = 1, n
      piv = a(i, bw + 1)
      do j = 1, min(bw, n - i)
        factor = a(i + j, bw + 1 - j)/piv
        a(i + j, bw + 1 - j) = factor
        do t = 1, min(bw, n - i)
          if (t - j > bw) cycle
          a(i + j, bw + 1 + t - j) = a(i + j, bw + 1 + t - j) - factor*a(i, bw + 1 + t)
        end do
      end do
    end do
  end subroutine ys_factor_banded_complex

  subroutine ys_solve_factored_banded_complex(rhs, a)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: a(:, :)
    integer(C_INT), parameter :: bw = 5
    integer(C_INT) :: n, i, j

    n = size(a, 1)
    do i = 1, n
      do j = max(1_C_INT, i - bw), i - 1
        rhs(i) = rhs(i) - a(i, bw + 1 + j - i)*rhs(j)
      end do
    end do

    do i = n, 1, -1
      do j = i + 1, min(n, i + bw)
        rhs(i) = rhs(i) - a(i, bw + 1 + j - i)*rhs(j)
      end do
      rhs(i) = rhs(i)/a(i, bw + 1)
    end do
  end subroutine ys_solve_factored_banded_complex

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
  !$omp end declare target

end module y_line_solvers
