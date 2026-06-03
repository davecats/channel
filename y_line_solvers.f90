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
  public :: ys_prepare_ghost_field_workspace, ys_release_ghost_field_workspace, ys_solve_ghost_field_reduced
  public :: ys_local_rhs, ys_local_operator
  public :: ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs
  public :: ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row

  integer(C_INT), save :: ys_workspace_ny = -1
  integer(C_INT), save :: ys_workspace_nz = -1
  integer(C_INT), save :: ys_workspace_nx = -1
  integer(C_INT), save :: ys_workspace_nlines = 0
  integer(C_INT), save :: ys_workspace_active_n = 0
  integer(C_INT), save :: ys_workspace_npy = -1
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_local_rhs(:, :)
  real(C_DOUBLE), allocatable, save :: ys_local_operator(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_lower_ghost_rhs(:), ys_lower_boundary_rhs(:), ys_upper_boundary_rhs(:), ys_upper_ghost_rhs(:)
  real(C_DOUBLE), allocatable, save :: ys_lower_ghost_row(:, :), ys_lower_boundary_row(:, :), ys_upper_boundary_row(:, :), ys_upper_ghost_row(:, :)
  real(C_DOUBLE), allocatable, save :: ys_interior_lu(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_interior_response_columns(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rows_send(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_left_interface_values(:, :), ys_right_interface_values(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rows_recv(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_matrix_lu(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rhs(:, :)

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
    if (allocated(ys_local_rhs)) then
      if (ys_workspace_ny /= ny .or. ys_workspace_nz /= nz .or. ys_workspace_nx /= nx_lines .or. &
          ys_workspace_active_n /= active_n .or. ys_workspace_npy /= npy_grid) then
        !$omp target exit data map(delete: ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
        !$omp& ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ys_interior_lu, ys_interior_response_columns, &
        !$omp& ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
  deallocate (ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs)
        deallocate (ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row)
deallocate (ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values)
        deallocate (ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
      end if
    end if

    if (.not. allocated(ys_local_rhs)) then
      allocate (ys_local_rhs(ny0:nyN, nlines), ys_local_operator(ny0:nyN, -2:2, nlines))
     allocate (ys_lower_ghost_rhs(nlines), ys_lower_boundary_rhs(nlines), ys_upper_boundary_rhs(nlines), ys_upper_ghost_rhs(nlines))
      allocate (ys_lower_ghost_row(-2:2, nlines), ys_lower_boundary_row(-2:2, nlines), ys_upper_boundary_row(-2:2, nlines), ys_upper_ghost_row(-2:2, nlines))
      allocate (ys_interior_lu(0:active_n - 1, -2:2, nlines), ys_interior_response_columns(0:active_n - 1, 5, nlines))
      allocate (ys_reduced_rows_send(20, nlines), ys_left_interface_values(2, nlines), ys_right_interface_values(2, nlines))
      allocate (ys_reduced_rows_recv(20, nlines, npy_grid), ys_reduced_matrix_lu(4*npy_grid, 11, nlines), ys_reduced_rhs(4*npy_grid, nlines))
      !$omp target enter data map(alloc: ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
      !$omp& ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ys_interior_lu, ys_interior_response_columns, &
      !$omp& ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
    end if

    ys_workspace_ny = ny
    ys_workspace_nz = nz
    ys_workspace_nx = nx_lines
    ys_workspace_nlines = nlines
    ys_workspace_active_n = active_n
    ys_workspace_npy = npy_grid

    ys_local_rhs = (0.0d0, 0.0d0)
    ys_local_operator = 0.0d0
    ys_lower_ghost_rhs = (0.0d0, 0.0d0)
    ys_lower_boundary_rhs = (0.0d0, 0.0d0)
    ys_upper_boundary_rhs = (0.0d0, 0.0d0)
    ys_upper_ghost_rhs = (0.0d0, 0.0d0)
    ys_lower_ghost_row = 0.0d0
    ys_lower_boundary_row = 0.0d0
    ys_upper_boundary_row = 0.0d0
    ys_upper_ghost_row = 0.0d0
    ys_interior_lu = 0.0d0
    ys_interior_response_columns = (0.0d0, 0.0d0)
    ys_reduced_rows_send = (0.0d0, 0.0d0)
    ys_left_interface_values = (0.0d0, 0.0d0)
    ys_right_interface_values = (0.0d0, 0.0d0)
    ys_reduced_rows_recv = (0.0d0, 0.0d0)
    ys_reduced_matrix_lu = (0.0d0, 0.0d0)
    ys_reduced_rhs = (0.0d0, 0.0d0)
  end subroutine ys_prepare_ghost_field_workspace

  subroutine ys_release_ghost_field_workspace()
    implicit none

    if (.not. allocated(ys_local_rhs)) return

    !$omp target exit data map(delete: ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
    !$omp& ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ys_interior_lu, ys_interior_response_columns, &
    !$omp& ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
  deallocate (ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs)
    deallocate (ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row)
deallocate (ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values)
    deallocate (ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)

    ys_workspace_ny = -1
    ys_workspace_nz = -1
    ys_workspace_nx = -1
    ys_workspace_nlines = 0
    ys_workspace_active_n = 0
    ys_workspace_npy = -1
  end subroutine ys_release_ghost_field_workspace

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

  subroutine ys_solve_ghost_field_reduced(dst, ny, nz)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)
    complex(C_DOUBLE_COMPLEX) :: lower_rhs0, upper_rhsn
    real(C_DOUBLE) :: row_coeffs(-2:2), lower_eq0(-2:2), upper_eqn(-2:2)
    integer(C_INT) :: ix, iz, iline, nlines, nlines_z, row_start, row_end, active_n, row, local_idx, col, coupled_row
    integer(C_INT) :: dst_row_base, lower_inner0, lower_inner2, upper_inner0, upper_inner2, upper_inner3
    logical :: has_lower_boundary, has_upper_boundary, has_padded_dst

    row_start = ny0
    row_end = nyN
    active_n = row_end - row_start + 1
    has_lower_boundary = (row_start == 1)
    has_upper_boundary = (row_end == ny - 1)
    has_padded_dst = (size(dst, 1) == active_n + 4)
    if (.not. has_padded_dst .and. size(dst, 1) /= active_n) then
      error stop "ys_solve_ghost_field_reduced expected either local-only or ghost-padded dst"
    end if
    dst_row_base = 1
    if (has_padded_dst) dst_row_base = 3
    if (.not. allocated(ys_local_rhs)) error stop "ys_prepare_ghost_field_workspace must be called before ys_solve_ghost_field"
    if (ys_workspace_ny /= ny .or. ys_workspace_nz /= nz .or. ys_workspace_nx /= size(dst, 3) .or. &
        ys_workspace_active_n /= active_n) then
      error stop "ys_solve_ghost_field workspace does not match requested local solve dimensions"
    end if

    nlines_z = 2*nz + 1
    nlines = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
    !$omp& ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ys_interior_lu, ys_interior_response_columns, &
    !$omp& ys_reduced_rows_send, row_start, row_end, active_n, has_lower_boundary, has_upper_boundary, nlines) &
    !$omp private(iline, row, local_idx, col, coupled_row, lower_rhs0, upper_rhsn, row_coeffs, lower_eq0, upper_eqn)
    do iline = 1, nlines
      ys_interior_lu(:, :, iline) = 0.0d0
      ys_interior_response_columns(:, :, iline) = (0.0d0, 0.0d0)

      lower_rhs0 = (0.0d0, 0.0d0)
      lower_eq0 = 0.0d0
      if (has_lower_boundary) then
lower_rhs0 = ys_lower_boundary_rhs(iline) - ys_lower_ghost_rhs(iline)*ys_lower_boundary_row(-2, iline)/ys_lower_ghost_row(-2, iline)
        lower_eq0 = ys_lower_boundary_row(:, iline) - ys_lower_ghost_row(:, iline)*ys_lower_boundary_row(-2, iline)/ys_lower_ghost_row(-2, iline)
        lower_eq0(-2) = 0.0d0
      end if
      upper_rhsn = (0.0d0, 0.0d0)
      upper_eqn = 0.0d0
      if (has_upper_boundary) then
  upper_rhsn = ys_upper_boundary_rhs(iline) - ys_upper_ghost_rhs(iline)*ys_upper_boundary_row(2, iline)/ys_upper_ghost_row(2, iline)
        upper_eqn = ys_upper_boundary_row(:, iline) - ys_upper_ghost_row(:, iline)*ys_upper_boundary_row(2, iline)/ys_upper_ghost_row(2, iline)
        upper_eqn(2) = 0.0d0
      end if

      do row = row_start, row_end
        local_idx = row - row_start
        ys_interior_response_columns(local_idx, 1, iline) = ys_local_rhs(row, iline)
      end do

      if (has_lower_boundary) then
        ys_interior_response_columns(0, 1, iline) = ys_interior_response_columns(0, 1, iline) - &
                                     ys_lower_ghost_rhs(iline)*ys_local_operator(row_start, -2, iline)/ys_lower_ghost_row(-2, iline)
        row_coeffs = ys_local_operator(row_start, -2:2, iline) - ys_lower_ghost_row(:, iline)*ys_local_operator(row_start, -2, iline)/ys_lower_ghost_row(-2, iline)
        row_coeffs(-2) = 0.0d0
     ys_interior_response_columns(0, 1, iline) = ys_interior_response_columns(0, 1, iline) - lower_rhs0*row_coeffs(-1)/lower_eq0(-1)
        row_coeffs = row_coeffs - lower_eq0*row_coeffs(-1)/lower_eq0(-1)
        row_coeffs(-1) = 0.0d0
        ys_interior_lu(0, 0:2, iline) = row_coeffs(0:2)

        ys_interior_response_columns(1, 1, iline) = ys_interior_response_columns(1, 1, iline) - lower_rhs0*ys_local_operator(row_start + 1, -2, iline)/lower_eq0(-1)
        row_coeffs = ys_local_operator(row_start + 1, -2:2, iline)
        row_coeffs(-2:1) = row_coeffs(-2:1) - lower_eq0(-1:2)*ys_local_operator(row_start + 1, -2, iline)/lower_eq0(-1)
        row_coeffs(-2) = 0.0d0
        ys_interior_lu(1, -1:2, iline) = row_coeffs(-1:2)
      end if

      if (has_upper_boundary) then
        ys_interior_response_columns(active_n - 1, 1, iline) = ys_interior_response_columns(active_n - 1, 1, iline) - &
                                         ys_upper_ghost_rhs(iline)*ys_local_operator(row_end, 2, iline)/ys_upper_ghost_row(2, iline)
        row_coeffs = ys_local_operator(row_end, -2:2, iline) - ys_upper_ghost_row(:, iline)*ys_local_operator(row_end, 2, iline)/ys_upper_ghost_row(2, iline)
        row_coeffs(2) = 0.0d0
        ys_interior_response_columns(active_n - 1, 1, iline) = ys_interior_response_columns(active_n - 1, 1, iline) - upper_rhsn*row_coeffs(1)/upper_eqn(1)
        row_coeffs = row_coeffs - upper_eqn*row_coeffs(1)/upper_eqn(1)
        row_coeffs(1) = 0.0d0
        ys_interior_lu(active_n - 1, -2:0, iline) = row_coeffs(-2:0)

        ys_interior_response_columns(active_n - 2, 1, iline) = ys_interior_response_columns(active_n - 2, 1, iline) - &
                                                               upper_rhsn*ys_local_operator(row_end - 1, 2, iline)/upper_eqn(1)
        row_coeffs = ys_local_operator(row_end - 1, -2:2, iline)
        row_coeffs(-1:2) = row_coeffs(-1:2) - upper_eqn(-2:1)*ys_local_operator(row_end - 1, 2, iline)/upper_eqn(1)
        row_coeffs(2) = 0.0d0
        ys_interior_lu(active_n - 2, -2:1, iline) = row_coeffs(-2:1)
      end if

      do row = row_start, row_end
        local_idx = row - row_start
        if ((has_lower_boundary .and. row <= row_start + 1) .or. (has_upper_boundary .and. row >= row_end - 1)) cycle
        row_coeffs = ys_local_operator(row, -2:2, iline)
        do col = -2, 2
          coupled_row = row + col
          if (coupled_row < row_start) then
 if (coupled_row == row_start - 2) ys_interior_response_columns(local_idx, 2, iline) = -cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
 if (coupled_row == row_start - 1) ys_interior_response_columns(local_idx, 3, iline) = -cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
          else if (coupled_row > row_end) then
   if (coupled_row == row_end + 1) ys_interior_response_columns(local_idx, 4, iline) = -cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
   if (coupled_row == row_end + 2) ys_interior_response_columns(local_idx, 5, iline) = -cmplx(row_coeffs(col), 0.0d0, kind=C_DOUBLE)
          else
            ys_interior_lu(local_idx, col, iline) = row_coeffs(col)
          end if
        end do
      end do

      call ys_factor_penta(ys_interior_lu(:, :, iline))
      call ys_solve_factored_penta_multi(ys_interior_response_columns(:, :, iline), ys_interior_lu(:, :, iline))

      ys_reduced_rows_send(1:4, iline) = [ys_interior_response_columns(0, 1, iline), ys_interior_response_columns(1, 1, iline), &
                         ys_interior_response_columns(active_n - 2, 1, iline), ys_interior_response_columns(active_n - 1, 1, iline)]
      ys_reduced_rows_send(5:12, iline) = [ys_interior_response_columns(0, 2, iline), ys_interior_response_columns(1, 2, iline), &
                       ys_interior_response_columns(active_n - 2, 2, iline), ys_interior_response_columns(active_n - 1, 2, iline), &
                                           ys_interior_response_columns(0, 3, iline), ys_interior_response_columns(1, 3, iline), &
                         ys_interior_response_columns(active_n - 2, 3, iline), ys_interior_response_columns(active_n - 1, 3, iline)]
      ys_reduced_rows_send(13:20, iline) = [ys_interior_response_columns(0, 4, iline), ys_interior_response_columns(1, 4, iline), &
                       ys_interior_response_columns(active_n - 2, 4, iline), ys_interior_response_columns(active_n - 1, 4, iline), &
                                            ys_interior_response_columns(0, 5, iline), ys_interior_response_columns(1, 5, iline), &
                         ys_interior_response_columns(active_n - 2, 5, iline), ys_interior_response_columns(active_n - 1, 5, iline)]
    end do
    !$omp end target teams distribute parallel do

    lower_inner0 = dst_row_base
    lower_inner2 = dst_row_base + 2
    upper_inner0 = active_n + dst_row_base - 3
    upper_inner2 = upper_inner0 + 2
    upper_inner3 = upper_inner0 + 3

    call ys_solve_reduced_interfaces()
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(dst, ys_interior_response_columns, ys_left_interface_values, ys_right_interface_values, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, &
    !$omp& ys_upper_ghost_rhs, ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, nx0, nxN, nz, nlines_z, dst_row_base, &
    !$omp& has_padded_dst, has_lower_boundary, has_upper_boundary, active_n, &
    !$omp& lower_inner0, lower_inner2, upper_inner0, upper_inner2, upper_inner3) &
    !$omp private(ix, iz, iline, row, lower_rhs0, upper_rhsn, lower_eq0, upper_eqn)
    do ix = nx0, nxN
      do iz = -nz, nz
        iline = (ix - nx0)*nlines_z + (iz + nz + 1)

        ys_interior_response_columns(:, 1, iline) = ys_interior_response_columns(:, 1, iline) + &
                                                    ys_interior_response_columns(:, 2, iline)*ys_left_interface_values(1, iline) + &
                                                    ys_interior_response_columns(:, 3, iline)*ys_left_interface_values(2, iline) + &
                                                   ys_interior_response_columns(:, 4, iline)*ys_right_interface_values(1, iline) + &
                                                    ys_interior_response_columns(:, 5, iline)*ys_right_interface_values(2, iline)

        do row = 0, active_n - 1
          dst(row + dst_row_base, iz + nz + 1, ix - nx0 + 1) = ys_interior_response_columns(row, 1, iline)
        end do

        if (has_padded_dst) then
          if (has_lower_boundary) then
lower_rhs0 = ys_lower_boundary_rhs(iline) - ys_lower_ghost_rhs(iline)*ys_lower_boundary_row(-2, iline)/ys_lower_ghost_row(-2, iline)
            lower_eq0 = ys_lower_boundary_row(:, iline) - ys_lower_ghost_row(:, iline)*ys_lower_boundary_row(-2, iline)/ys_lower_ghost_row(-2, iline)
            lower_eq0(-2) = 0.0d0

            dst(2, iz + nz + 1, ix - nx0 + 1) = (lower_rhs0 - sum(lower_eq0(0:2)*dst(lower_inner0:lower_inner2, iz + nz + 1, ix - nx0 + 1)))/lower_eq0(-1)
            dst(1, iz + nz + 1, ix - nx0 + 1) = (ys_lower_ghost_rhs(iline) - &
                             sum(ys_lower_ghost_row(-1:2, iline)*dst(2:5, iz + nz + 1, ix - nx0 + 1)))/ys_lower_ghost_row(-2, iline)
          else
            dst(1, iz + nz + 1, ix - nx0 + 1) = ys_left_interface_values(1, iline)
            dst(2, iz + nz + 1, ix - nx0 + 1) = ys_left_interface_values(2, iline)
          end if

          if (has_upper_boundary) then
  upper_rhsn = ys_upper_boundary_rhs(iline) - ys_upper_ghost_rhs(iline)*ys_upper_boundary_row(2, iline)/ys_upper_ghost_row(2, iline)
            upper_eqn = ys_upper_boundary_row(:, iline) - ys_upper_ghost_row(:, iline)*ys_upper_boundary_row(2, iline)/ys_upper_ghost_row(2, iline)
            upper_eqn(2) = 0.0d0

            dst(active_n + 3, iz + nz + 1, ix - nx0 + 1) = (upper_rhsn - &
                                        sum(upper_eqn(-2:0)*dst(upper_inner0:upper_inner2, iz + nz + 1, ix - nx0 + 1)))/upper_eqn(1)
            dst(active_n + 4, iz + nz + 1, ix - nx0 + 1) = (ys_upper_ghost_rhs(iline) - &
        sum(ys_upper_ghost_row(-2:1, iline)*dst(upper_inner0:upper_inner3, iz + nz + 1, ix - nx0 + 1)))/ys_upper_ghost_row(2, iline)
          else
            dst(active_n + 3, iz + nz + 1, ix - nx0 + 1) = ys_right_interface_values(1, iline)
            dst(active_n + 4, iz + nz + 1, ix - nx0 + 1) = ys_right_interface_values(2, iline)
          end if
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    !$omp target update from(dst)

  end subroutine ys_solve_ghost_field_reduced

  subroutine ys_solve_reduced_interfaces()
    integer(C_INT), parameter :: bw = 5
    complex(C_DOUBLE_COMPLEX) :: solve_piv, solve_factor, packed_remote(20)
    integer(C_INT) :: nlines, niface, iline, iblock, row0, i, j, t

    nlines = size(ys_reduced_rows_send, 2)
    niface = 4*npy_grid

#ifdef HAVE_MPI
#ifndef HAVE_HIP
    !$omp target data use_device_ptr(ys_reduced_rows_send, ys_reduced_rows_recv)
#endif
    call MPI_Allgather(ys_reduced_rows_send, 20*nlines, MPI_DOUBLE_COMPLEX, ys_reduced_rows_recv, 20*nlines, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr)
#ifndef HAVE_HIP
    !$omp end target data
#endif
#else
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_reduced_rows_send, ys_reduced_rows_recv, nlines) private(iline, i)
    do iline = 1, nlines
      do i = 1, 20
        ys_reduced_rows_recv(i, iline, 1) = ys_reduced_rows_send(i, iline)
      end do
    end do
    !$omp end target teams distribute parallel do
#endif

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs, ys_left_interface_values, ys_right_interface_values, nlines, niface, npy_grid, ipy) &
    !$omp private(iline, iblock, row0, packed_remote, i, j, t, solve_piv, solve_factor)
    do iline = 1, nlines
      ys_reduced_matrix_lu(:, :, iline) = (0.0d0, 0.0d0)
      ys_reduced_rhs(:, iline) = (0.0d0, 0.0d0)
      do iblock = 0, npy_grid - 1
        row0 = 4*iblock
        packed_remote = ys_reduced_rows_recv(:, iline, iblock + 1)
        ys_reduced_matrix_lu(row0 + 1, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 2, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 3, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 4, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_rhs(row0 + 1:row0 + 4, iline) = packed_remote(1:4)
        if (iblock > 0) then
          ys_reduced_matrix_lu(row0 + 1, bw - 1, iline) = -packed_remote(5)
          ys_reduced_matrix_lu(row0 + 2, bw - 2, iline) = -packed_remote(6)
          ys_reduced_matrix_lu(row0 + 3, bw - 3, iline) = -packed_remote(7)
          ys_reduced_matrix_lu(row0 + 4, bw - 4, iline) = -packed_remote(8)
          ys_reduced_matrix_lu(row0 + 1, bw, iline) = -packed_remote(9)
          ys_reduced_matrix_lu(row0 + 2, bw - 1, iline) = -packed_remote(10)
          ys_reduced_matrix_lu(row0 + 3, bw - 2, iline) = -packed_remote(11)
          ys_reduced_matrix_lu(row0 + 4, bw - 3, iline) = -packed_remote(12)
        end if
        if (iblock < npy_grid - 1) then
          ys_reduced_matrix_lu(row0 + 1, bw + 5, iline) = -packed_remote(13)
          ys_reduced_matrix_lu(row0 + 2, bw + 4, iline) = -packed_remote(14)
          ys_reduced_matrix_lu(row0 + 3, bw + 3, iline) = -packed_remote(15)
          ys_reduced_matrix_lu(row0 + 4, bw + 2, iline) = -packed_remote(16)
          ys_reduced_matrix_lu(row0 + 1, bw + 6, iline) = -packed_remote(17)
          ys_reduced_matrix_lu(row0 + 2, bw + 5, iline) = -packed_remote(18)
          ys_reduced_matrix_lu(row0 + 3, bw + 4, iline) = -packed_remote(19)
          ys_reduced_matrix_lu(row0 + 4, bw + 3, iline) = -packed_remote(20)
        end if
      end do

      call ys_factor_banded_complex(ys_reduced_matrix_lu(:, :, iline))
      call ys_solve_factored_banded_complex(ys_reduced_rhs(:, iline), ys_reduced_matrix_lu(:, :, iline))

      ys_left_interface_values(:, iline) = (0.0d0, 0.0d0)
      ys_right_interface_values(:, iline) = (0.0d0, 0.0d0)
      row0 = 4*ipy
      if (ipy > 0) then
        ys_left_interface_values(:, iline) = ys_reduced_rhs(row0 - 1:row0, iline)
      end if
      if (ipy < npy_grid - 1) then
        ys_right_interface_values(:, iline) = ys_reduced_rhs(row0 + 5:row0 + 6, iline)
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_solve_reduced_interfaces

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_factor_banded_complex)
#endif
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

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_solve_factored_banded_complex)
#endif
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

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_factor_penta)
#endif
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

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_solve_factored_penta_multi)
#endif
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
