#include "header.h"

module y_line_solvers

  use, intrinsic :: iso_c_binding
  use mpi_transpose, only: ny0, nyN, nx0, nxN, nxB, npy_grid, ipy, allgather_y_device_complex_rows
  use roctx, only: roctxPush, roctxPop
#ifdef HAVE_CUDA
  use cusparse
#endif

  implicit none
  private
  integer(C_INT), parameter :: YS_ENDPOINT_RESPONSE_CONST = 1_C_INT
  integer(C_INT), parameter :: YS_ENDPOINT_RESPONSE_EVEN_Z = 2_C_INT

  public :: ys_prepare_assembled_workspace, ys_release_workspace
  public :: ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs
  public :: ys_eqm1, ys_eq0, ys_eqn, ys_eqnp1
  public :: ys_boundary_lower_rhs0, ys_boundary_upper_rhsn, ys_boundary_lower_eq, ys_boundary_upper_eq
  public :: ys_boundary_lower_rhs0_owner, ys_boundary_upper_rhsn_owner, ys_boundary_lower_eq_owner, ys_boundary_upper_eq_owner
  public :: ys_lower_ghost_owner, ys_lower_boundary_owner, ys_upper_boundary_owner, ys_upper_ghost_owner
  public :: ys_eqm1_owner, ys_eq0_owner, ys_eqn_owner, ys_eqnp1_owner
  public :: ys_gpsv_matrix, ys_gpsv_rhs, ys_gpsv_line_matrix, ys_gpsv_line_rhs
  public :: ys_gpsv_owner_matrix, ys_gpsv_owner_rhs, ys_owner_nz, ys_owner_nx, ys_owner_ix0, ys_owner_ixN
  public :: ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x
  public :: ys_solve_packed_pentadiagonal, ys_solve_endpoint_schur

  integer(C_INT), save :: ys_workspace_ny = -1
  integer(C_INT), save :: ys_workspace_nz = -1
  integer(C_INT), save :: ys_workspace_nx = -1
  integer(C_INT), save :: ys_workspace_nlines = 0
  integer(C_INT), save :: ys_workspace_active_n = 0
  integer(C_INT), save :: ys_workspace_npy = -1
  integer(C_INT), save :: ys_workspace_line_start = 1
  integer(C_INT), save :: ys_workspace_line_end = 0
  integer(C_INT), save :: ys_owner_nz = 0
  integer(C_INT), save :: ys_owner_nx = 0
  integer(C_INT), save :: ys_owner_ix0 = 1
  integer(C_INT), save :: ys_owner_ixN = 0
  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: ys_gpsv_matrix_store(:), ys_gpsv_rhs_store(:)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_gpsv_ds(:), ys_gpsv_dl(:), ys_gpsv_d(:), ys_gpsv_du(:), ys_gpsv_dw(:), ys_gpsv_x(:)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_gpsv_matrix(:, :, :, :), ys_gpsv_rhs(:, :, :)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_gpsv_line_matrix(:, :, :), ys_gpsv_line_rhs(:, :)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_gpsv_owner_matrix(:, :, :, :), ys_gpsv_owner_rhs(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: ys_lower_ghost_store(:), ys_lower_boundary_store(:), ys_upper_boundary_store(:), ys_upper_ghost_store(:)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_lower_ghost_rhs(:), ys_lower_boundary_rhs(:), ys_upper_boundary_rhs(:), ys_upper_ghost_rhs(:)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_lower_ghost_owner(:, :), ys_lower_boundary_owner(:, :), ys_upper_boundary_owner(:, :), ys_upper_ghost_owner(:, :)
  real(C_DOUBLE), allocatable, target, save :: ys_eqm1_store(:), ys_eq0_store(:), ys_eqn_store(:), ys_eqnp1_store(:)
  real(C_DOUBLE), pointer, save :: ys_eqm1(:, :), ys_eq0(:, :), ys_eqn(:, :), ys_eqnp1(:, :)
  real(C_DOUBLE), pointer, save :: ys_eqm1_owner(:, :, :), ys_eq0_owner(:, :, :), ys_eqn_owner(:, :, :), ys_eqnp1_owner(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: ys_boundary_lower_rhs0_store(:), ys_boundary_upper_rhsn_store(:)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_boundary_lower_rhs0(:), ys_boundary_upper_rhsn(:)
  real(C_DOUBLE), allocatable, target, save :: ys_boundary_lower_eq_store(:), ys_boundary_upper_eq_store(:)
  real(C_DOUBLE), pointer, save :: ys_boundary_lower_eq(:, :), ys_boundary_upper_eq(:, :)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_boundary_lower_rhs0_owner(:, :), ys_boundary_upper_rhsn_owner(:, :)
  real(C_DOUBLE), pointer, save :: ys_boundary_lower_eq_owner(:, :, :), ys_boundary_upper_eq_owner(:, :, :)
  real(C_DOUBLE), allocatable, save :: ys_interior_lu(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_interior_response_columns(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rows_send(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_left_interface_values(:, :), ys_right_interface_values(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rows_recv(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_matrix_lu(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rhs(:, :)
#ifdef HAVE_CUDA
  type(cusparseHandle), save :: ys_gpsv_handle
  logical, save :: ys_gpsv_handle_created = .false.
#endif
  integer(C_INT), save :: ys_gpsv_n = -1, ys_gpsv_batch = -1
  integer(C_INT), save :: ys_batch_n = -1, ys_batch_count = -1
#ifdef HAVE_CUDA
  integer(8), save :: ys_gpsv_buffer_size = 0_8
#endif
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_batch_ds(:), ys_batch_dl(:), ys_batch_d(:), ys_batch_du(:), ys_batch_dw(:), ys_batch_x(:)
#ifdef HAVE_CUDA
  character(c_char), allocatable, save :: ys_gpsv_buffer(:)
#endif

contains

  subroutine ys_release_core_workspace()
    implicit none

    if (.not. allocated(ys_gpsv_rhs_store)) return

    !$omp target exit data map(delete: ys_gpsv_matrix_store, ys_gpsv_rhs_store, &
    !$omp& ys_lower_ghost_store, ys_lower_boundary_store, ys_upper_boundary_store, ys_upper_ghost_store, &
    !$omp& ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store, ys_boundary_lower_rhs0_store, ys_boundary_upper_rhsn_store, &
    !$omp& ys_boundary_lower_eq_store, ys_boundary_upper_eq_store)
    nullify (ys_gpsv_matrix, ys_gpsv_rhs, ys_gpsv_line_matrix, ys_gpsv_line_rhs, ys_gpsv_owner_matrix, ys_gpsv_owner_rhs, &
             ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, ys_lower_ghost_rhs, ys_lower_boundary_rhs, &
             ys_upper_boundary_rhs, ys_upper_ghost_rhs, ys_lower_ghost_owner, ys_lower_boundary_owner, ys_upper_boundary_owner, ys_upper_ghost_owner, &
             ys_eqm1, ys_eq0, ys_eqn, ys_eqnp1, ys_eqm1_owner, ys_eq0_owner, ys_eqn_owner, ys_eqnp1_owner, ys_boundary_lower_rhs0, &
   ys_boundary_upper_rhsn, ys_boundary_lower_eq, ys_boundary_upper_eq, ys_boundary_lower_rhs0_owner, ys_boundary_upper_rhsn_owner, &
             ys_boundary_lower_eq_owner, ys_boundary_upper_eq_owner)
    deallocate (ys_gpsv_matrix_store, ys_gpsv_rhs_store)
    deallocate (ys_lower_ghost_store, ys_lower_boundary_store, ys_upper_boundary_store, ys_upper_ghost_store)
    deallocate (ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store)
    deallocate (ys_boundary_lower_rhs0_store, ys_boundary_upper_rhsn_store, ys_boundary_lower_eq_store, ys_boundary_upper_eq_store)
  end subroutine ys_release_core_workspace

  subroutine ys_release_reduced_workspace()
    implicit none

    if (.not. allocated(ys_interior_lu)) return

    !$omp target exit data map(delete: ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, &
    !$omp& ys_right_interface_values, ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
deallocate (ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values)
    deallocate (ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
  end subroutine ys_release_reduced_workspace

  subroutine ys_allocate_core_workspace(row_start, row_end, line_start, nlines)
    implicit none
    integer(C_INT), intent(in) :: row_start, row_end, line_start, nlines
    integer(C_INT) :: active_n, line_end

    active_n = row_end - row_start + 1
    line_end = line_start + nlines - 1
    allocate (ys_gpsv_matrix_store(5*nlines*active_n), ys_gpsv_rhs_store(nlines*active_n))
    allocate (ys_lower_ghost_store(nlines), ys_lower_boundary_store(nlines), ys_upper_boundary_store(nlines), ys_upper_ghost_store(nlines))
    allocate (ys_eqm1_store(5*nlines), ys_eq0_store(5*nlines), ys_eqn_store(5*nlines), ys_eqnp1_store(5*nlines))
    allocate (ys_boundary_lower_rhs0_store(nlines), ys_boundary_upper_rhsn_store(nlines))
    allocate (ys_boundary_lower_eq_store(4*nlines), ys_boundary_upper_eq_store(4*nlines))
    ys_gpsv_ds(1:nlines*active_n) => ys_gpsv_matrix_store(1:nlines*active_n)
    ys_gpsv_dl(1:nlines*active_n) => ys_gpsv_matrix_store(nlines*active_n + 1:2*nlines*active_n)
    ys_gpsv_d(1:nlines*active_n) => ys_gpsv_matrix_store(2*nlines*active_n + 1:3*nlines*active_n)
    ys_gpsv_du(1:nlines*active_n) => ys_gpsv_matrix_store(3*nlines*active_n + 1:4*nlines*active_n)
    ys_gpsv_dw(1:nlines*active_n) => ys_gpsv_matrix_store(4*nlines*active_n + 1:5*nlines*active_n)
    ys_gpsv_x(1:nlines*active_n) => ys_gpsv_rhs_store(1:nlines*active_n)
    ys_gpsv_line_matrix(line_start:line_end, row_start:row_end, -2:2) => ys_gpsv_matrix_store
    ys_gpsv_line_rhs(line_start:line_end, row_start:row_end) => ys_gpsv_rhs_store
  ys_gpsv_owner_matrix(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN, row_start:row_end, -2:2) => ys_gpsv_matrix_store
    ys_gpsv_owner_rhs(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN, row_start:row_end) => ys_gpsv_rhs_store
    ys_lower_ghost_rhs(line_start:line_end) => ys_lower_ghost_store
    ys_lower_boundary_rhs(line_start:line_end) => ys_lower_boundary_store
    ys_upper_boundary_rhs(line_start:line_end) => ys_upper_boundary_store
    ys_upper_ghost_rhs(line_start:line_end) => ys_upper_ghost_store
    ys_lower_ghost_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_lower_ghost_store
    ys_lower_boundary_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_lower_boundary_store
    ys_upper_boundary_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_upper_boundary_store
    ys_upper_ghost_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_upper_ghost_store
    ys_eqm1(-2:2, line_start:line_end) => ys_eqm1_store
    ys_eq0(-2:2, line_start:line_end) => ys_eq0_store
    ys_eqn(-2:2, line_start:line_end) => ys_eqn_store
    ys_eqnp1(-2:2, line_start:line_end) => ys_eqnp1_store
    ys_eqm1_owner(-2:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_eqm1_store
    ys_eq0_owner(-2:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_eq0_store
    ys_eqn_owner(-2:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_eqn_store
    ys_eqnp1_owner(-2:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_eqnp1_store
    ys_boundary_lower_rhs0(line_start:line_end) => ys_boundary_lower_rhs0_store
    ys_boundary_upper_rhsn(line_start:line_end) => ys_boundary_upper_rhsn_store
    ys_boundary_lower_eq(-1:2, line_start:line_end) => ys_boundary_lower_eq_store
    ys_boundary_upper_eq(-2:1, line_start:line_end) => ys_boundary_upper_eq_store
    ys_boundary_lower_rhs0_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_boundary_lower_rhs0_store
    ys_boundary_upper_rhsn_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_boundary_upper_rhsn_store
    ys_boundary_lower_eq_owner(-1:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_boundary_lower_eq_store
    ys_boundary_upper_eq_owner(-2:1, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_boundary_upper_eq_store
    !$omp target enter data map(alloc: ys_gpsv_matrix_store, ys_gpsv_rhs_store, &
    !$omp& ys_lower_ghost_store, ys_lower_boundary_store, ys_upper_boundary_store, ys_upper_ghost_store, &
    !$omp& ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store, ys_boundary_lower_rhs0_store, ys_boundary_upper_rhsn_store, &
    !$omp& ys_boundary_lower_eq_store, ys_boundary_upper_eq_store)
  end subroutine ys_allocate_core_workspace

  subroutine ys_allocate_reduced_workspace(active_n, nlines, npy_count)
    implicit none
    integer(C_INT), intent(in) :: active_n, nlines, npy_count

    allocate (ys_interior_lu(0:active_n - 1, -2:2, nlines), ys_interior_response_columns(0:active_n - 1, 5, nlines))
    allocate (ys_reduced_rows_send(20, nlines), ys_left_interface_values(2, nlines), ys_right_interface_values(2, nlines))
    allocate (ys_reduced_rows_recv(20, nlines, npy_count), ys_reduced_matrix_lu(4*npy_count, 11, nlines), ys_reduced_rhs(4*npy_count, nlines))
    !$omp target enter data map(alloc: ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, &
    !$omp& ys_right_interface_values, ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
  end subroutine ys_allocate_reduced_workspace

  subroutine ys_prepare_assembled_workspace(ny, nz, row_start, row_end, line_start, nlines, use_reduced_backend)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, row_start, row_end, line_start, nlines
    logical, intent(in) :: use_reduced_backend
    integer(C_INT) :: active_n, line_end, nlines_z

    active_n = row_end - row_start + 1
    line_end = line_start + nlines - 1
    nlines_z = 2*nz + 1
    if (active_n < 1) error stop "ys_prepare_assembled_workspace requires at least one row"
    if (mod(line_start - 1, nlines_z) /= 0) error stop "ys_prepare_assembled_workspace requires ix-aligned line_start"
    if (mod(nlines, nlines_z) /= 0) error stop "ys_prepare_assembled_workspace requires full ix columns"
    if (use_reduced_backend) then
      if (npy_grid > 1 .and. active_n < 4) error stop "ys_solve_ghost_field requires at least four y rows per rank"
      if (active_n < 4) error stop "ys_solve_ghost_field requires at least four active y rows"
    end if

    ys_owner_nz = nlines_z
    ys_owner_nx = nlines/nlines_z
    ys_owner_ix0 = nx0 + (line_start - 1)/nlines_z
    ys_owner_ixN = ys_owner_ix0 + ys_owner_nx - 1
    if (allocated(ys_gpsv_rhs_store)) then
      if (ys_workspace_ny /= ny .or. ys_workspace_nz /= nz .or. ys_workspace_nlines /= nlines .or. &
          ys_workspace_active_n /= active_n .or. ys_workspace_npy /= merge(npy_grid, -1_C_INT, use_reduced_backend) .or. &
          ys_workspace_line_start /= line_start) then
        call ys_release_core_workspace()
        call ys_release_reduced_workspace()
      end if
    end if

    if (.not. allocated(ys_gpsv_rhs_store)) then
      call ys_allocate_core_workspace(row_start, row_end, line_start, nlines)
    end if
    if (use_reduced_backend) then
      if (.not. allocated(ys_interior_lu)) call ys_allocate_reduced_workspace(active_n, nlines, npy_grid)
    else if (allocated(ys_interior_lu)) then
      call ys_release_reduced_workspace()
    end if

    ys_workspace_ny = ny
    ys_workspace_nz = nz
    ys_workspace_nx = nlines/(2*nz + 1)
    ys_workspace_nlines = nlines
    ys_workspace_active_n = active_n
    ys_workspace_npy = merge(npy_grid, -1_C_INT, use_reduced_backend)
    ys_workspace_line_start = line_start
    ys_workspace_line_end = line_end
    ys_gpsv_line_matrix(line_start:line_end, row_start:row_end, -2:2) => ys_gpsv_matrix_store
    ys_gpsv_line_rhs(line_start:line_end, row_start:row_end) => ys_gpsv_rhs_store
  ys_gpsv_owner_matrix(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN, row_start:row_end, -2:2) => ys_gpsv_matrix_store
    ys_gpsv_owner_rhs(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN, row_start:row_end) => ys_gpsv_rhs_store
    ys_lower_ghost_rhs(line_start:line_end) => ys_lower_ghost_store
    ys_lower_boundary_rhs(line_start:line_end) => ys_lower_boundary_store
    ys_upper_boundary_rhs(line_start:line_end) => ys_upper_boundary_store
    ys_upper_ghost_rhs(line_start:line_end) => ys_upper_ghost_store
    ys_lower_ghost_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_lower_ghost_store
    ys_lower_boundary_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_lower_boundary_store
    ys_upper_boundary_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_upper_boundary_store
    ys_upper_ghost_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_upper_ghost_store
    ys_eqm1(-2:2, line_start:line_end) => ys_eqm1_store
    ys_eq0(-2:2, line_start:line_end) => ys_eq0_store
    ys_eqn(-2:2, line_start:line_end) => ys_eqn_store
    ys_eqnp1(-2:2, line_start:line_end) => ys_eqnp1_store
    ys_eqm1_owner(-2:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_eqm1_store
    ys_eq0_owner(-2:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_eq0_store
    ys_eqn_owner(-2:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_eqn_store
    ys_eqnp1_owner(-2:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_eqnp1_store
    ys_boundary_lower_rhs0(line_start:line_end) => ys_boundary_lower_rhs0_store
    ys_boundary_upper_rhsn(line_start:line_end) => ys_boundary_upper_rhsn_store
    ys_boundary_lower_eq(-1:2, line_start:line_end) => ys_boundary_lower_eq_store
    ys_boundary_upper_eq(-2:1, line_start:line_end) => ys_boundary_upper_eq_store
    ys_boundary_lower_rhs0_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_boundary_lower_rhs0_store
    ys_boundary_upper_rhsn_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_boundary_upper_rhsn_store
    ys_boundary_lower_eq_owner(-1:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_boundary_lower_eq_store
    ys_boundary_upper_eq_owner(-2:1, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => ys_boundary_upper_eq_store
    if (line_start == 1_C_INT .and. nlines == (nxN - nx0 + 1)*(2*nz + 1)) then
      ys_gpsv_matrix(-nz:nz, nx0:nxN, row_start:row_end, -2:2) => ys_gpsv_matrix_store
      ys_gpsv_rhs(-nz:nz, nx0:nxN, row_start:row_end) => ys_gpsv_rhs_store
    else
      nullify (ys_gpsv_matrix, ys_gpsv_rhs)
    end if
  end subroutine ys_prepare_assembled_workspace

  subroutine ys_release_workspace()
    implicit none

    if (.not. allocated(ys_gpsv_rhs_store)) return

    call ys_release_core_workspace()
    call ys_release_reduced_workspace()
    call ys_release_gpsv_workspace()

    ys_workspace_ny = -1
    ys_workspace_nz = -1
    ys_workspace_nx = -1
    ys_workspace_nlines = 0
    ys_workspace_active_n = 0
    ys_workspace_npy = -1
    ys_workspace_line_start = 1
    ys_workspace_line_end = 0
    ys_owner_nz = 0
    ys_owner_nx = 0
    ys_owner_ix0 = 1
    ys_owner_ixN = 0
  end subroutine ys_release_workspace

#ifdef HAVE_CUDA
  subroutine ys_check_cusparse(status, where)
    integer(C_INT), intent(in) :: status
    character(*), intent(in) :: where

    if (status /= CUSPARSE_STATUS_SUCCESS) then
      print *, "cuSPARSE error in ", trim(where), ": status=", status
      error stop
    end if
  end subroutine ys_check_cusparse
#endif

  subroutine ys_release_gpsv_workspace()
    implicit none
#ifdef HAVE_CUDA
    integer(C_INT) :: status
#endif

    if (allocated(ys_batch_ds)) then
      !$omp target exit data map(delete: ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
      deallocate (ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
    end if
#ifdef HAVE_CUDA
    if (allocated(ys_gpsv_buffer)) then
      !$omp target exit data map(delete: ys_gpsv_buffer)
      deallocate (ys_gpsv_buffer)
    end if
    if (ys_gpsv_handle_created) then
      status = cusparseDestroy(ys_gpsv_handle)
      call ys_check_cusparse(status, "cusparseDestroy")
      ys_gpsv_handle_created = .false.
    end if
#endif
    ys_gpsv_n = -1
    ys_gpsv_batch = -1
    ys_batch_n = -1
    ys_batch_count = -1
#ifdef HAVE_CUDA
    ys_gpsv_buffer_size = 0_8
#endif
  end subroutine ys_release_gpsv_workspace

  subroutine ys_prepare_cusparse_workspace(n, batch_count, ds, dl, d, du, dw, x)
    implicit none
    integer(C_INT), intent(in) :: n, batch_count
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:), dw(:), x(:)
#ifdef HAVE_CUDA
    integer(C_INT) :: status
    integer(8) :: buffer_size
#endif

    if (ys_gpsv_n == n .and. ys_gpsv_batch == batch_count .and. allocated(ys_gpsv_buffer)) return

#ifdef HAVE_CUDA
    if (allocated(ys_gpsv_buffer)) then
      !$omp target exit data map(delete: ys_gpsv_buffer)
      deallocate (ys_gpsv_buffer)
    end if

    if (.not. ys_gpsv_handle_created) then
      status = cusparseCreate(ys_gpsv_handle)
      call ys_check_cusparse(status, "cusparseCreate")
      ys_gpsv_handle_created = .true.
    end if
    !$omp target data use_device_addr(ds, dl, d, du, dw, x)
    status = cusparseZgpsvInterleavedBatch_bufferSize(ys_gpsv_handle, 0_C_INT, n, ds, dl, d, du, dw, x, batch_count, buffer_size)
    !$omp end target data
    call ys_check_cusparse(status, "cusparseZgpsvInterleavedBatch_bufferSize")
    ys_gpsv_buffer_size = buffer_size
    allocate (ys_gpsv_buffer(max(1, int(buffer_size))))
    !$omp target enter data map(alloc: ys_gpsv_buffer)
#endif

    ys_gpsv_n = n
    ys_gpsv_batch = batch_count
  end subroutine ys_prepare_cusparse_workspace

  subroutine ys_prepare_batch_workspace(n, batch_count)
    implicit none
    integer(C_INT), intent(in) :: n, batch_count

    if (ys_batch_n /= n .or. ys_batch_count /= batch_count) then
      if (allocated(ys_batch_ds)) then
        !$omp target exit data map(delete: ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
        deallocate (ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
      end if
      allocate (ys_batch_ds(n*batch_count), ys_batch_dl(n*batch_count), ys_batch_d(n*batch_count), &
                ys_batch_du(n*batch_count), ys_batch_dw(n*batch_count), ys_batch_x(n*batch_count))
      !$omp target enter data map(alloc: ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
      ys_batch_n = n
      ys_batch_count = batch_count
    end if

    call ys_prepare_cusparse_workspace(n, batch_count, ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
  end subroutine ys_prepare_batch_workspace

  subroutine ys_solve_interleaved_pentadiagonal(ds, dl, d, du, dw, x, n, batch_count, label)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:), dw(:), x(:)
    integer(C_INT), intent(in) :: n, batch_count
    character(*), intent(in) :: label
#ifdef HAVE_CUDA
    integer(C_INT) :: status
#else
    integer(C_INT) :: iline
#endif

    call roctxPush(label)
#ifdef HAVE_CUDA
    call ys_prepare_cusparse_workspace(n, batch_count, ds, dl, d, du, dw, x)
    !$omp target data use_device_addr(ds, dl, d, du, dw, x, ys_gpsv_buffer)
    status = cusparseZgpsvInterleavedBatch(ys_gpsv_handle, 0_C_INT, n, ds, dl, d, du, dw, x, batch_count, ys_gpsv_buffer)
    !$omp end target data
    call ys_check_cusparse(status, "cusparseZgpsvInterleavedBatch "//trim(label))
#else
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ds, dl, d, du, dw, x, n, batch_count) &
    !$omp private(iline)
    do iline = 1, batch_count
      call ys_factor_penta_interleaved(ds, dl, d, du, dw, batch_count, iline, n)
      call ys_solve_factored_penta_interleaved(x, ds, dl, d, du, dw, batch_count, iline, n)
    end do
    !$omp end target teams distribute parallel do
#endif
    call roctxPop(label)
  end subroutine ys_solve_interleaved_pentadiagonal

  subroutine ys_solve_packed_pentadiagonal(n, batch_count, label)
    implicit none
    integer(C_INT), intent(in) :: n, batch_count
    character(*), intent(in) :: label

    call ys_solve_interleaved_pentadiagonal(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, &
                                            n, batch_count, label)
  end subroutine ys_solve_packed_pentadiagonal

  subroutine ys_solve_endpoint_schur(dst, symmetric_operator)
    implicit none
    logical, intent(in) :: symmetric_operator
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)
    complex(C_DOUBLE_COMPLEX) :: rhs_value
    complex(C_DOUBLE_COMPLEX) :: c0, c1, y00, y01, y10, y11, g0, g1
    complex(C_DOUBLE_COMPLEX) :: s00, s01, s10, s11, t00, t01, t10, t11, det
    complex(C_DOUBLE_COMPLEX) :: vleft1, vleft2, vright1, vright2
    complex(C_DOUBLE_COMPLEX) :: s4(4, 4), rhs4(4, 5), pivot4, factor4
    complex(C_DOUBLE_COMPLEX) :: row_coeffs(-2:2), coeff
    integer(C_INT) :: row_start, row_end, active_n, nlines, nlines_z, nz, dst_row_base, response_mode
    integer(C_INT) :: nI, nresp, batch_count, exposed_n, interior_base
    integer(C_INT) :: sys, iline, ref_iline, resp, resp_index
    integer(C_INT) :: local_i, local_idx, row, col, coupled_row, p, j, offset, actual_p
    integer(C_INT) :: exposed_slot, response_slot, rhs_col, iface, k, m
    integer(C_INT) :: ix, iz, abs_iz, ix_local
    logical :: is_actual
    logical :: has_left_interface, has_right_interface, has_padded_dst

    row_start = ny0
    row_end = nyN
    active_n = ys_workspace_active_n
    nlines = ys_workspace_nlines
    nlines_z = size(dst, 2)
    nz = (nlines_z - 1)/2
    has_padded_dst = (size(dst, 1) == active_n + 4)
    dst_row_base = merge(3_C_INT, 1_C_INT, has_padded_dst)
    response_mode = merge(YS_ENDPOINT_RESPONSE_EVEN_Z, YS_ENDPOINT_RESPONSE_CONST, symmetric_operator)

    has_left_interface = (ipy > 0)
    has_right_interface = (ipy < npy_grid - 1)
    exposed_n = 0
    if (has_left_interface) exposed_n = exposed_n + 2
    if (has_right_interface) exposed_n = exposed_n + 2
    interior_base = 0
    if (has_left_interface) interior_base = 2
    nI = active_n - exposed_n
    if (nI <= 0) error stop "endpoint Schur solve needs at least one interior row"
    select case (response_mode)
    case (YS_ENDPOINT_RESPONSE_CONST)
      nresp = 1_C_INT
    case (YS_ENDPOINT_RESPONSE_EVEN_Z)
      nresp = (nxN - nx0 + 1)*(nz + 1)
    case default
      error stop "unknown endpoint Schur response mode"
    end select
    batch_count = nlines + exposed_n*nresp

    call ys_prepare_batch_workspace(nI, batch_count)

    call roctxPush("ys_endpoint_pack_plus_response")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x, &
    !$omp& row_start, active_n, nI, nlines, nlines_z, nresp, batch_count, nz, nx0, &
    !$omp& response_mode, ipy, has_left_interface, has_right_interface, exposed_n, interior_base) &
    !$omp private(sys, iline, ref_iline, resp, local_i, local_idx, row, col, coupled_row, p, j, offset, actual_p, is_actual, ix_local, abs_iz, &
    !$omp& rhs_value, row_coeffs, coeff, exposed_slot, response_slot)
    do local_i = 0, nI - 1
      do sys = 1, batch_count
        is_actual = (sys <= nlines)
        if (is_actual) then
          iline = sys
          ref_iline = iline
        else
          resp = mod(sys - nlines - 1, nresp) + 1
          select case (response_mode)
          case (YS_ENDPOINT_RESPONSE_CONST)
            ref_iline = 1_C_INT
          case default
            ix_local = (resp - 1)/(nz + 1)
            abs_iz = mod(resp - 1, nz + 1)
            ref_iline = ix_local*nlines_z + (abs_iz + nz + 1)
          end select
          iline = ref_iline
        end if
        local_idx = interior_base + local_i
        row = row_start + local_idx
        rhs_value = (0.0d0, 0.0d0)
        if (is_actual) then
          actual_p = local_idx*nlines + iline
          rhs_value = ys_gpsv_x(actual_p)
        end if
        actual_p = local_idx*nlines + ref_iline
        row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)

        p = local_i*batch_count + sys
        ys_batch_ds(p) = (0.0d0, 0.0d0)
        ys_batch_dl(p) = (0.0d0, 0.0d0)
        ys_batch_d(p) = (0.0d0, 0.0d0)
        ys_batch_du(p) = (0.0d0, 0.0d0)
        ys_batch_dw(p) = (0.0d0, 0.0d0)

        do col = -2, 2
          coeff = row_coeffs(col)
          if (coeff == 0.0d0) cycle
          coupled_row = local_idx + col
          if (coupled_row >= interior_base .and. coupled_row <= interior_base + nI - 1) then
            j = coupled_row - interior_base
            offset = j - local_i
          else
            exposed_slot = 0
            if (has_left_interface) then
              if (coupled_row == 0) exposed_slot = 1
              if (coupled_row == 1) exposed_slot = 2
            end if
            if (has_right_interface) then
              if (coupled_row == active_n - 2) exposed_slot = exposed_n - 1
              if (coupled_row == active_n - 1) exposed_slot = exposed_n
            end if
            if (exposed_slot > 0) then
              if (sys > nlines) then
                response_slot = (sys - nlines - 1)/nresp + 1
                if (response_slot == exposed_slot) rhs_value = coeff
              end if
            end if
            cycle
          end if

          select case (offset)
          case (-2)
            ys_batch_ds(p) = coeff
          case (-1)
            ys_batch_dl(p) = coeff
          case (0)
            ys_batch_d(p) = coeff
          case (1)
            ys_batch_du(p) = coeff
          case (2)
            ys_batch_dw(p) = coeff
          end select
        end do
        ys_batch_x(p) = rhs_value
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_endpoint_pack_plus_response")

    call ys_solve_interleaved_pentadiagonal(ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x, &
                                            nI, batch_count, "ys_endpoint_gpsv_plus_response")

    call roctxPush("ys_endpoint_pack_schur")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, ys_batch_x, ys_reduced_rows_send, nI, nlines, nlines_z, nresp, batch_count, &
    !$omp& nz, nx0, row_start, active_n, response_mode, ipy, has_left_interface, has_right_interface, exposed_n, interior_base) &
    !$omp private(iline, row, ix, iz, abs_iz, resp_index, c0, c1, y00, y01, y10, y11, g0, g1, s00, s01, s10, s11, &
    !$omp& t00, t01, t10, t11, det, row_coeffs, s4, rhs4, pivot4, factor4, local_idx, col, coeff, coupled_row, j, &
    !$omp& exposed_slot, iface, k, m, rhs_col, actual_p)
    do iline = 1, nlines
      ix = (iline - 1)/nlines_z + nx0
      iz = mod(iline - 1, nlines_z) - nz
      select case (response_mode)
      case (YS_ENDPOINT_RESPONSE_CONST)
        resp_index = 1_C_INT
      case default
        abs_iz = abs(iz)
        resp_index = (ix - nx0)*(nz + 1) + abs_iz + 1
      end select
      ys_reduced_rows_send(:, iline) = (0.0d0, 0.0d0)
      if (.not. has_left_interface) then
        row = row_start + nI
        actual_p = (row - row_start)*nlines + iline
        row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)
        c0 = ys_batch_x((nI - 2)*batch_count + iline)
        c1 = ys_batch_x((nI - 1)*batch_count + iline)
        y00 = ys_batch_x((nI - 2)*batch_count + nlines + resp_index)
        y10 = ys_batch_x((nI - 1)*batch_count + nlines + resp_index)
        y01 = ys_batch_x((nI - 2)*batch_count + nlines + nresp + resp_index)
        y11 = ys_batch_x((nI - 1)*batch_count + nlines + nresp + resp_index)
        g0 = ys_gpsv_x(actual_p) - row_coeffs(-2)*c0 - row_coeffs(-1)*c1
        s00 = row_coeffs(0) - row_coeffs(-2)*y00 - row_coeffs(-1)*y10
        s01 = row_coeffs(1) - row_coeffs(-2)*y01 - row_coeffs(-1)*y11
        t00 = row_coeffs(2)
        t01 = (0.0d0, 0.0d0)

        row = row_start + nI + 1
        actual_p = (row - row_start)*nlines + iline
        row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)
        c1 = ys_batch_x((nI - 1)*batch_count + iline)
        y10 = ys_batch_x((nI - 1)*batch_count + nlines + resp_index)
        y11 = ys_batch_x((nI - 1)*batch_count + nlines + nresp + resp_index)
        g1 = ys_gpsv_x(actual_p) - row_coeffs(-2)*c1
        s10 = row_coeffs(-1) - row_coeffs(-2)*y10
        s11 = row_coeffs(0) - row_coeffs(-2)*y11
        t10 = row_coeffs(1)
        t11 = row_coeffs(2)

        det = s00*s11 - s01*s10
        ys_reduced_rows_send(3, iline) = (s11*g0 - s01*g1)/det
        ys_reduced_rows_send(4, iline) = (-s10*g0 + s00*g1)/det
        ys_reduced_rows_send(15, iline) = -((s11*t00 - s01*t10)/det)
        ys_reduced_rows_send(16, iline) = -((-s10*t00 + s00*t10)/det)
        ys_reduced_rows_send(19, iline) = -((s11*t01 - s01*t11)/det)
        ys_reduced_rows_send(20, iline) = -((-s10*t01 + s00*t11)/det)
      else if (.not. has_right_interface) then
        row = row_start
        actual_p = iline
        row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)
        c0 = ys_batch_x(iline)
        y00 = ys_batch_x(nlines + resp_index)
        y01 = ys_batch_x(nlines + nresp + resp_index)
        g0 = ys_gpsv_x(actual_p) - row_coeffs(2)*c0
        s00 = row_coeffs(0) - row_coeffs(2)*y00
        s01 = row_coeffs(1) - row_coeffs(2)*y01
        t00 = row_coeffs(-2)
        t01 = row_coeffs(-1)

        row = row_start + 1
        actual_p = nlines + iline
        row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)
        c0 = ys_batch_x(iline)
        c1 = ys_batch_x(batch_count + iline)
        y00 = ys_batch_x(nlines + resp_index)
        y10 = ys_batch_x(batch_count + nlines + resp_index)
        y01 = ys_batch_x(nlines + nresp + resp_index)
        y11 = ys_batch_x(batch_count + nlines + nresp + resp_index)
        g1 = ys_gpsv_x(actual_p) - row_coeffs(1)*c0 - row_coeffs(2)*c1
        s10 = row_coeffs(-1) - row_coeffs(1)*y00 - row_coeffs(2)*y10
        s11 = row_coeffs(0) - row_coeffs(1)*y01 - row_coeffs(2)*y11
        t10 = (0.0d0, 0.0d0)
        t11 = row_coeffs(-2)

        det = s00*s11 - s01*s10
        ys_reduced_rows_send(1, iline) = (s11*g0 - s01*g1)/det
        ys_reduced_rows_send(2, iline) = (-s10*g0 + s00*g1)/det
        ys_reduced_rows_send(5, iline) = -((s11*t00 - s01*t10)/det)
        ys_reduced_rows_send(6, iline) = -((-s10*t00 + s00*t10)/det)
        ys_reduced_rows_send(9, iline) = -((s11*t01 - s01*t11)/det)
        ys_reduced_rows_send(10, iline) = -((-s10*t01 + s00*t11)/det)
      else
        do k = 1, 4
          do m = 1, 4
            s4(k, m) = (0.0d0, 0.0d0)
          end do
          do rhs_col = 1, 5
            rhs4(k, rhs_col) = (0.0d0, 0.0d0)
          end do
        end do
        do iface = 1, 4
          select case (iface)
          case (1)
            local_idx = 0
          case (2)
            local_idx = 1
          case (3)
            local_idx = active_n - 2
          case default
            local_idx = active_n - 1
          end select
          row = row_start + local_idx
          actual_p = local_idx*nlines + iline
        row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)
          rhs4(iface, 1) = ys_gpsv_x(actual_p)
          do col = -2, 2
            coeff = row_coeffs(col)
            if (coeff == 0.0d0) cycle
            coupled_row = local_idx + col
            if (coupled_row >= interior_base .and. coupled_row <= interior_base + nI - 1) then
              j = coupled_row - interior_base
              rhs4(iface, 1) = rhs4(iface, 1) - coeff*ys_batch_x(j*batch_count + iline)
              do exposed_slot = 1, exposed_n
                s4(iface, exposed_slot) = s4(iface, exposed_slot) - &
                                          coeff*ys_batch_x(j*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)
              end do
            else if (coupled_row == -2) then
              rhs4(iface, 2) = rhs4(iface, 2) + coeff
            else if (coupled_row == -1) then
              rhs4(iface, 3) = rhs4(iface, 3) + coeff
            else if (coupled_row == active_n) then
              rhs4(iface, 4) = rhs4(iface, 4) + coeff
            else if (coupled_row == active_n + 1) then
              rhs4(iface, 5) = rhs4(iface, 5) + coeff
            else
              exposed_slot = 0
              if (coupled_row == 0) exposed_slot = 1
              if (coupled_row == 1) exposed_slot = 2
              if (coupled_row == active_n - 2) exposed_slot = 3
              if (coupled_row == active_n - 1) exposed_slot = 4
              if (exposed_slot > 0) s4(iface, exposed_slot) = s4(iface, exposed_slot) + coeff
            end if
          end do
        end do

        do k = 1, 4
          pivot4 = s4(k, k)
          do m = k + 1, 4
            factor4 = s4(m, k)/pivot4
            s4(m, k) = factor4
            do j = k + 1, 4
              s4(m, j) = s4(m, j) - factor4*s4(k, j)
            end do
            do rhs_col = 1, 5
              rhs4(m, rhs_col) = rhs4(m, rhs_col) - factor4*rhs4(k, rhs_col)
            end do
          end do
        end do
        do rhs_col = 1, 5
          do k = 4, 1, -1
            do j = k + 1, 4
              rhs4(k, rhs_col) = rhs4(k, rhs_col) - s4(k, j)*rhs4(j, rhs_col)
            end do
            rhs4(k, rhs_col) = rhs4(k, rhs_col)/s4(k, k)
          end do
        end do

        do iface = 1, 4
          ys_reduced_rows_send(iface, iline) = rhs4(iface, 1)
          ys_reduced_rows_send(4 + iface, iline) = -rhs4(iface, 2)
          ys_reduced_rows_send(8 + iface, iline) = -rhs4(iface, 3)
          ys_reduced_rows_send(12 + iface, iline) = -rhs4(iface, 4)
          ys_reduced_rows_send(16 + iface, iline) = -rhs4(iface, 5)
        end do
      end if
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_endpoint_pack_schur")

    call ys_solve_reduced_interfaces()

    call roctxPush("ys_endpoint_reconstruct")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(dst, ys_batch_x, ys_reduced_rhs, ys_left_interface_values, ys_right_interface_values, nlines, nlines_z, nresp, nx0, nz, &
    !$omp& active_n, dst_row_base, has_padded_dst, nI, batch_count, response_mode, ipy, has_left_interface, has_right_interface, interior_base) &
    !$omp private(iline, local_i, local_idx, p, ix, iz, abs_iz, resp_index, vleft1, vleft2, vright1, vright2, &
    !$omp& exposed_slot)
    do iline = 1, nlines
      ix = (iline - 1)/nlines_z + nx0
      iz = mod(iline - 1, nlines_z) - nz
      select case (response_mode)
      case (YS_ENDPOINT_RESPONSE_CONST)
        resp_index = 1_C_INT
      case default
        abs_iz = abs(iz)
        resp_index = (ix - nx0)*(nz + 1) + abs_iz + 1
      end select
      vleft1 = (0.0d0, 0.0d0)
      vleft2 = (0.0d0, 0.0d0)
      vright1 = (0.0d0, 0.0d0)
      vright2 = (0.0d0, 0.0d0)
      if (has_left_interface) then
        vleft1 = ys_reduced_rhs(4*ipy + 1, iline)
        vleft2 = ys_reduced_rhs(4*ipy + 2, iline)
        dst(dst_row_base, iz + nz + 1, ix - nx0 + 1) = vleft1
        dst(dst_row_base + 1, iz + nz + 1, ix - nx0 + 1) = vleft2
      end if
      if (has_right_interface) then
        vright1 = ys_reduced_rhs(4*ipy + 3, iline)
        vright2 = ys_reduced_rhs(4*ipy + 4, iline)
        dst(active_n - 2 + dst_row_base, iz + nz + 1, ix - nx0 + 1) = vright1
        dst(active_n - 1 + dst_row_base, iz + nz + 1, ix - nx0 + 1) = vright2
      end if
      do local_i = 0, nI - 1
        p = local_i*batch_count + iline
        local_idx = interior_base + local_i
        dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = ys_batch_x(p)
        exposed_slot = 0
        if (has_left_interface) then
          exposed_slot = exposed_slot + 1
          dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = &
            dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) - &
            ys_batch_x(local_i*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)*vleft1
          exposed_slot = exposed_slot + 1
          dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = &
            dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) - &
            ys_batch_x(local_i*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)*vleft2
        end if
        if (has_right_interface) then
          exposed_slot = exposed_slot + 1
          dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = &
            dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) - &
            ys_batch_x(local_i*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)*vright1
          exposed_slot = exposed_slot + 1
          dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = &
            dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) - &
            ys_batch_x(local_i*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)*vright2
        end if
      end do
      if (has_padded_dst) then
        if (has_left_interface) then
          dst(1, iz + nz + 1, ix - nx0 + 1) = ys_left_interface_values(1, iline)
          dst(2, iz + nz + 1, ix - nx0 + 1) = ys_left_interface_values(2, iline)
        end if
        if (has_right_interface) then
          dst(active_n + 3, iz + nz + 1, ix - nx0 + 1) = ys_right_interface_values(1, iline)
          dst(active_n + 4, iz + nz + 1, ix - nx0 + 1) = ys_right_interface_values(2, iline)
        end if
      end if
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_endpoint_reconstruct")
  end subroutine ys_solve_endpoint_schur

  subroutine ys_solve_reduced_interfaces()
    integer(C_INT), parameter :: bw = 5
    complex(C_DOUBLE_COMPLEX) :: solve_piv, solve_factor, packed_remote(20)
    integer(C_INT) :: nlines, niface, iline, iblock, row0, i, j, t

    nlines = size(ys_reduced_rows_send, 2)
    niface = 4*npy_grid

    call allgather_y_device_complex_rows(ys_reduced_rows_send, ys_reduced_rows_recv, 20_C_INT, nlines, &
                                         "MPI_Allgather reduced_y_interfaces")

    call roctxPush("ys_reduced_interfaces_solve")
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
    call roctxPop("ys_reduced_interfaces_solve")
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
  !$omp declare target(ys_factor_penta_interleaved)
#endif
  subroutine ys_factor_penta_interleaved(ds, dl, d, du, dw, stride, first, n)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i, p, p1, p2
    complex(C_DOUBLE_COMPLEX) :: factor

    do i = 0, n - 1
      p = first + i*stride
      d(p) = 1.0d0/d(p)

      if (i + 1 < n) then
        p1 = first + (i + 1)*stride
        factor = dl(p1)*d(p)
        dl(p1) = factor
        d(p1) = d(p1) - factor*du(p)
        if (i + 2 < n) du(p1) = du(p1) - factor*dw(p)
      end if

      if (i + 2 < n) then
        p2 = first + (i + 2)*stride
        factor = ds(p2)*d(p)
        ds(p2) = factor
        dl(p2) = dl(p2) - factor*du(p)
        d(p2) = d(p2) - factor*dw(p)
      end if
    end do
  end subroutine ys_factor_penta_interleaved

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_solve_factored_penta_interleaved)
#endif
  subroutine ys_solve_factored_penta_interleaved(rhs, ds, dl, d, du, dw, stride, first, n)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i, p

    do i = 0, n - 1
      p = first + i*stride
      if (i >= 2) rhs(p) = rhs(p) - ds(p)*rhs(first + (i - 2)*stride)
      if (i >= 1) rhs(p) = rhs(p) - dl(p)*rhs(first + (i - 1)*stride)
    end do

    do i = n - 1, 0, -1
      p = first + i*stride
      if (i + 1 < n) rhs(p) = rhs(p) - du(p)*rhs(first + (i + 1)*stride)
      if (i + 2 < n) rhs(p) = rhs(p) - dw(p)*rhs(first + (i + 2)*stride)
      rhs(p) = rhs(p)*d(p)
    end do
  end subroutine ys_solve_factored_penta_interleaved
end module y_line_solvers
