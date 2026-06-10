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
  integer(C_INT), parameter :: YS_REDUCED_BW = 5_C_INT

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

  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: ys_lower_ghost_store(:), ys_lower_boundary_store(:)
  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: ys_upper_boundary_store(:), ys_upper_ghost_store(:)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_lower_ghost_rhs(:), ys_lower_boundary_rhs(:)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_upper_boundary_rhs(:), ys_upper_ghost_rhs(:)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_lower_ghost_owner(:, :), ys_lower_boundary_owner(:, :)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_upper_boundary_owner(:, :), ys_upper_ghost_owner(:, :)

  real(C_DOUBLE), allocatable, target, save :: ys_eqm1_store(:), ys_eq0_store(:), ys_eqn_store(:), ys_eqnp1_store(:)
  real(C_DOUBLE), pointer, save :: ys_eqm1(:, :), ys_eq0(:, :), ys_eqn(:, :), ys_eqnp1(:, :)
  real(C_DOUBLE), pointer, save :: ys_eqm1_owner(:, :, :), ys_eq0_owner(:, :, :), ys_eqn_owner(:, :, :), ys_eqnp1_owner(:, :, :)

  complex(C_DOUBLE_COMPLEX), allocatable, target, save :: ys_boundary_lower_rhs0_store(:), ys_boundary_upper_rhsn_store(:)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_boundary_lower_rhs0(:), ys_boundary_upper_rhsn(:)
  real(C_DOUBLE), allocatable, target, save :: ys_boundary_lower_eq_store(:), ys_boundary_upper_eq_store(:)
  real(C_DOUBLE), pointer, save :: ys_boundary_lower_eq(:, :), ys_boundary_upper_eq(:, :)
  complex(C_DOUBLE_COMPLEX), pointer, save :: ys_boundary_lower_rhs0_owner(:, :), ys_boundary_upper_rhsn_owner(:, :)
  real(C_DOUBLE), pointer, save :: ys_boundary_lower_eq_owner(:, :, :), ys_boundary_upper_eq_owner(:, :, :)

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
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_batch_ds(:), ys_batch_dl(:), ys_batch_d(:)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_batch_du(:), ys_batch_dw(:), ys_batch_x(:)
#ifdef HAVE_CUDA
  character(c_char), allocatable, save :: ys_gpsv_buffer(:)
#endif

contains

  subroutine ys_reset_workspace_state()
    implicit none

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
  end subroutine ys_reset_workspace_state

  logical function ys_core_shape_matches(ny, nz, active_n, nlines, line_start, wanted_npy) result(matches)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, active_n, nlines, line_start, wanted_npy

    matches = allocated(ys_gpsv_rhs_store) .and. &
              ys_workspace_ny == ny .and. ys_workspace_nz == nz .and. &
              ys_workspace_nlines == nlines .and. ys_workspace_active_n == active_n .and. &
              ys_workspace_npy == wanted_npy .and. ys_workspace_line_start == line_start
  end function ys_core_shape_matches

  subroutine ys_nullify_core_views()
    implicit none

    nullify (ys_gpsv_matrix, ys_gpsv_rhs, ys_gpsv_line_matrix, ys_gpsv_line_rhs, &
             ys_gpsv_owner_matrix, ys_gpsv_owner_rhs, ys_gpsv_ds, ys_gpsv_dl, &
             ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, ys_lower_ghost_rhs, &
             ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
             ys_lower_ghost_owner, ys_lower_boundary_owner, ys_upper_boundary_owner, &
             ys_upper_ghost_owner, ys_eqm1, ys_eq0, ys_eqn, ys_eqnp1, &
             ys_eqm1_owner, ys_eq0_owner, ys_eqn_owner, ys_eqnp1_owner, &
             ys_boundary_lower_rhs0, ys_boundary_upper_rhsn, ys_boundary_lower_eq, &
             ys_boundary_upper_eq, ys_boundary_lower_rhs0_owner, &
             ys_boundary_upper_rhsn_owner, ys_boundary_lower_eq_owner, &
             ys_boundary_upper_eq_owner)
  end subroutine ys_nullify_core_views

  subroutine ys_bind_core_views(row_start, row_end, line_start, nlines)
    implicit none
    integer(C_INT), intent(in) :: row_start, row_end, line_start, nlines
    integer(C_INT) :: active_n, line_end, nall

    active_n = row_end - row_start + 1
    line_end = line_start + nlines - 1
    nall = nlines*active_n

    ys_gpsv_ds(1:nall) => ys_gpsv_matrix_store(1:nall)
    ys_gpsv_dl(1:nall) => ys_gpsv_matrix_store(nall + 1:2*nall)
    ys_gpsv_d(1:nall) => ys_gpsv_matrix_store(2*nall + 1:3*nall)
    ys_gpsv_du(1:nall) => ys_gpsv_matrix_store(3*nall + 1:4*nall)
    ys_gpsv_dw(1:nall) => ys_gpsv_matrix_store(4*nall + 1:5*nall)
    ys_gpsv_x(1:nall) => ys_gpsv_rhs_store(1:nall)

    ys_gpsv_line_matrix(line_start:line_end, row_start:row_end, -2:2) => ys_gpsv_matrix_store
    ys_gpsv_line_rhs(line_start:line_end, row_start:row_end) => ys_gpsv_rhs_store
    ys_gpsv_owner_matrix(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN, &
                         row_start:row_end, -2:2) => ys_gpsv_matrix_store
    ys_gpsv_owner_rhs(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN, &
                      row_start:row_end) => ys_gpsv_rhs_store

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
    ys_boundary_lower_rhs0_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => &
      ys_boundary_lower_rhs0_store
    ys_boundary_upper_rhsn_owner(-ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => &
      ys_boundary_upper_rhsn_store
    ys_boundary_lower_eq_owner(-1:2, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => &
      ys_boundary_lower_eq_store
    ys_boundary_upper_eq_owner(-2:1, -ys_workspace_nz:ys_workspace_nz, ys_owner_ix0:ys_owner_ixN) => &
      ys_boundary_upper_eq_store

    if (line_start == 1_C_INT .and. nlines == (nxN - nx0 + 1)*(2*ys_workspace_nz + 1)) then
      ys_gpsv_matrix(-ys_workspace_nz:ys_workspace_nz, nx0:nxN, row_start:row_end, -2:2) => &
        ys_gpsv_matrix_store
      ys_gpsv_rhs(-ys_workspace_nz:ys_workspace_nz, nx0:nxN, row_start:row_end) => ys_gpsv_rhs_store
    else
      nullify (ys_gpsv_matrix, ys_gpsv_rhs)
    end if
  end subroutine ys_bind_core_views

  subroutine ys_release_core_workspace()
    implicit none

    if (.not. allocated(ys_gpsv_rhs_store)) return

    !$omp target exit data map(delete: ys_gpsv_matrix_store, ys_gpsv_rhs_store, &
    !$omp& ys_lower_ghost_store, ys_lower_boundary_store, ys_upper_boundary_store, ys_upper_ghost_store, &
    !$omp& ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store, &
    !$omp& ys_boundary_lower_rhs0_store, ys_boundary_upper_rhsn_store, &
    !$omp& ys_boundary_lower_eq_store, ys_boundary_upper_eq_store)

    call ys_nullify_core_views()

    deallocate (ys_gpsv_matrix_store, ys_gpsv_rhs_store)
    deallocate (ys_lower_ghost_store, ys_lower_boundary_store, ys_upper_boundary_store, ys_upper_ghost_store)
    deallocate (ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store)
    deallocate (ys_boundary_lower_rhs0_store, ys_boundary_upper_rhsn_store, &
                ys_boundary_lower_eq_store, ys_boundary_upper_eq_store)
  end subroutine ys_release_core_workspace

  subroutine ys_release_reduced_workspace()
    implicit none

    if (allocated(ys_reduced_rows_send)) then
      !$omp target exit data map(delete: ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, &
      !$omp& ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
      deallocate (ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values)
      deallocate (ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
    end if

  end subroutine ys_release_reduced_workspace

  subroutine ys_allocate_core_workspace(row_start, row_end, line_start, nlines)
    implicit none
    integer(C_INT), intent(in) :: row_start, row_end, line_start, nlines
    integer(C_INT) :: active_n

    active_n = row_end - row_start + 1

    allocate (ys_gpsv_matrix_store(5*nlines*active_n), ys_gpsv_rhs_store(nlines*active_n))
    allocate (ys_lower_ghost_store(nlines), ys_lower_boundary_store(nlines), &
              ys_upper_boundary_store(nlines), ys_upper_ghost_store(nlines))
    allocate (ys_eqm1_store(5*nlines), ys_eq0_store(5*nlines), &
              ys_eqn_store(5*nlines), ys_eqnp1_store(5*nlines))
    allocate (ys_boundary_lower_rhs0_store(nlines), ys_boundary_upper_rhsn_store(nlines))
    allocate (ys_boundary_lower_eq_store(4*nlines), ys_boundary_upper_eq_store(4*nlines))

    call ys_bind_core_views(row_start, row_end, line_start, nlines)

    !$omp target enter data map(alloc: ys_gpsv_matrix_store, ys_gpsv_rhs_store, &
    !$omp& ys_lower_ghost_store, ys_lower_boundary_store, ys_upper_boundary_store, ys_upper_ghost_store, &
    !$omp& ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store, &
    !$omp& ys_boundary_lower_rhs0_store, ys_boundary_upper_rhsn_store, &
    !$omp& ys_boundary_lower_eq_store, ys_boundary_upper_eq_store)
  end subroutine ys_allocate_core_workspace

  subroutine ys_allocate_reduced_workspace(nlines, npy_count)
    implicit none
    integer(C_INT), intent(in) :: nlines, npy_count

    allocate (ys_reduced_rows_send(20, nlines), ys_left_interface_values(2, nlines), &
              ys_right_interface_values(2, nlines))
    allocate (ys_reduced_rows_recv(20, nlines, npy_count), &
              ys_reduced_matrix_lu(4*npy_count, 2*YS_REDUCED_BW + 1, nlines), &
              ys_reduced_rhs(4*npy_count, nlines))

    !$omp target enter data map(alloc: ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, &
    !$omp& ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
  end subroutine ys_allocate_reduced_workspace

  subroutine ys_prepare_assembled_workspace(ny, nz, row_start, row_end, line_start, nlines, use_reduced_backend)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, row_start, row_end, line_start, nlines
    logical, intent(in) :: use_reduced_backend
    integer(C_INT) :: active_n, line_end, nlines_z, wanted_npy

    active_n = row_end - row_start + 1
    line_end = line_start + nlines - 1
    nlines_z = 2*nz + 1
    wanted_npy = merge(npy_grid, -1_C_INT, use_reduced_backend)

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
      if (.not. ys_core_shape_matches(ny, nz, active_n, nlines, line_start, wanted_npy)) then
        call ys_release_core_workspace()
        call ys_release_reduced_workspace()
      end if
    end if

    ys_workspace_ny = ny
    ys_workspace_nz = nz
    ys_workspace_nx = nlines/nlines_z
    ys_workspace_nlines = nlines
    ys_workspace_active_n = active_n
    ys_workspace_npy = wanted_npy
    ys_workspace_line_start = line_start
    ys_workspace_line_end = line_end

    if (.not. allocated(ys_gpsv_rhs_store)) then
      call ys_allocate_core_workspace(row_start, row_end, line_start, nlines)
    else
      call ys_bind_core_views(row_start, row_end, line_start, nlines)
    end if

    if (use_reduced_backend) then
      if (.not. allocated(ys_reduced_rows_send)) call ys_allocate_reduced_workspace(nlines, npy_grid)
    else if (allocated(ys_reduced_rows_send)) then
      call ys_release_reduced_workspace()
    end if
  end subroutine ys_prepare_assembled_workspace

  subroutine ys_release_workspace()
    implicit none

    call ys_release_core_workspace()
    call ys_release_reduced_workspace()
    call ys_release_gpsv_workspace()
    call ys_reset_workspace_state()
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

#ifndef HAVE_CUDA
    ys_gpsv_n = n
    ys_gpsv_batch = batch_count
#else
    if (ys_gpsv_n == n .and. ys_gpsv_batch == batch_count .and. allocated(ys_gpsv_buffer)) return

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
    status = cusparseZgpsvInterleavedBatch_bufferSize(ys_gpsv_handle, 0_C_INT, n, ds, dl, d, du, dw, x, &
                                                      batch_count, buffer_size)
    !$omp end target data
    call ys_check_cusparse(status, "cusparseZgpsvInterleavedBatch_bufferSize")

    ys_gpsv_buffer_size = buffer_size
    allocate (ys_gpsv_buffer(max(1, int(buffer_size))))
    !$omp target enter data map(alloc: ys_gpsv_buffer)

    ys_gpsv_n = n
    ys_gpsv_batch = batch_count
#endif
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

    call ys_prepare_cusparse_workspace(n, batch_count, ys_batch_ds, ys_batch_dl, ys_batch_d, &
                                       ys_batch_du, ys_batch_dw, ys_batch_x)
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
    status = cusparseZgpsvInterleavedBatch(ys_gpsv_handle, 0_C_INT, n, ds, dl, d, du, dw, x, &
                                           batch_count, ys_gpsv_buffer)
    !$omp end target data
    call ys_check_cusparse(status, "cusparseZgpsvInterleavedBatch "//trim(label))
#else
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ds, dl, d, du, dw, x, n, batch_count) private(iline)
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
    complex(C_DOUBLE_COMPLEX) :: rhs_value, coeff, row_coeffs(-2:2), s4(4, 4), rhs4(4, 5)
    complex(C_DOUBLE_COMPLEX) :: pivot4, factor4, iface_value(4)
    integer(C_INT) :: row_start, active_n, nlines, nlines_z, nz, nx_count, dst_row_base, response_mode
    integer(C_INT) :: nI, nresp, batch_count, exposed_n, left_count, interior_base
    integer(C_INT) :: sys, iline, ref_iline, resp, local_i, local_idx, col, coupled_row
    integer(C_INT) :: p, j, offset, actual_p, response_slot, exposed_slot, ghost_col
    integer(C_INT) :: ix_local, iz_index, iz, abs_iz, resp_index, iface, slot, rhs_col, k, m
    logical :: is_actual, has_left_interface, has_right_interface, has_padded_dst

    row_start = ny0
    active_n = ys_workspace_active_n
    nlines = ys_workspace_nlines
    nlines_z = size(dst, 2)
    nz = (nlines_z - 1)/2
    nx_count = nlines/nlines_z
    has_padded_dst = (size(dst, 1) == active_n + 4)
    dst_row_base = merge(3_C_INT, 1_C_INT, has_padded_dst)
    response_mode = merge(YS_ENDPOINT_RESPONSE_EVEN_Z, YS_ENDPOINT_RESPONSE_CONST, symmetric_operator)

    has_left_interface = (ipy > 0)
    has_right_interface = (ipy < npy_grid - 1)
    left_count = merge(2_C_INT, 0_C_INT, has_left_interface)
    exposed_n = left_count + merge(2_C_INT, 0_C_INT, has_right_interface)
    interior_base = left_count
    nI = active_n - exposed_n
    if (nI <= 0) error stop "endpoint Schur solve needs at least one interior row"

    select case (response_mode)
    case (YS_ENDPOINT_RESPONSE_CONST)
      nresp = 1_C_INT
    case (YS_ENDPOINT_RESPONSE_EVEN_Z)
      nresp = nx_count*(nz + 1)
    case default
      error stop "unknown endpoint Schur response mode"
    end select
    batch_count = nlines + exposed_n*nresp

    call ys_prepare_batch_workspace(nI, batch_count)

    call roctxPush("ys_endpoint_pack_plus_response")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, &
    !$omp& ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x, &
    !$omp& active_n, nI, nlines, nlines_z, nz, nresp, batch_count, response_mode, &
    !$omp& has_left_interface, has_right_interface, left_count, exposed_n, interior_base) &
    !$omp private(sys, iline, ref_iline, resp, local_i, local_idx, col, coupled_row, p, j, offset, &
    !$omp& actual_p, response_slot, exposed_slot, ix_local, abs_iz, rhs_value, coeff, row_coeffs, is_actual)
    do local_i = 0, nI - 1
      do sys = 1, batch_count
        is_actual = (sys <= nlines)
        if (is_actual) then
          iline = sys
          ref_iline = iline
          response_slot = 0_C_INT
        else
          resp = mod(sys - nlines - 1, nresp) + 1
          response_slot = (sys - nlines - 1)/nresp + 1
          select case (response_mode)
          case (YS_ENDPOINT_RESPONSE_CONST)
            ref_iline = 1_C_INT
          case default
            ix_local = (resp - 1)/(nz + 1)
            abs_iz = mod(resp - 1, nz + 1)
            ref_iline = ix_local*nlines_z + abs_iz + nz + 1
          end select
          iline = ref_iline
        end if

        local_idx = interior_base + local_i
        actual_p = local_idx*nlines + ref_iline
        row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), &
                       ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)
        rhs_value = (0.0d0, 0.0d0)
        if (is_actual) rhs_value = ys_gpsv_x(local_idx*nlines + iline)

        p = local_i*batch_count + sys
        ys_batch_ds(p) = (0.0d0, 0.0d0)
        ys_batch_dl(p) = (0.0d0, 0.0d0)
        ys_batch_d(p) = (0.0d0, 0.0d0)
        ys_batch_du(p) = (0.0d0, 0.0d0)
        ys_batch_dw(p) = (0.0d0, 0.0d0)

        do col = -2, 2
          coeff = row_coeffs(col)
          if (coeff == (0.0d0, 0.0d0)) cycle
          coupled_row = local_idx + col
          if (coupled_row >= interior_base .and. coupled_row <= interior_base + nI - 1) then
            j = coupled_row - interior_base
            offset = j - local_i
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
          else if (.not. is_actual) then
            exposed_slot = 0_C_INT
            if (has_left_interface) then
              if (coupled_row == 0) exposed_slot = 1_C_INT
              if (coupled_row == 1) exposed_slot = 2_C_INT
            end if
            if (has_right_interface) then
              if (coupled_row == active_n - 2) exposed_slot = left_count + 1
              if (coupled_row == active_n - 1) exposed_slot = left_count + 2
            end if
            if (exposed_slot == response_slot) rhs_value = coeff
          end if
        end do
        ys_batch_x(p) = rhs_value
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_endpoint_pack_plus_response")

    call ys_solve_interleaved_pentadiagonal(ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x, &
                                            nI, batch_count, "ys_endpoint_gpsv_plus_response")

    if (exposed_n > 0) then
      call roctxPush("ys_endpoint_pack_schur")
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, ys_batch_x, &
      !$omp& ys_reduced_rows_send, active_n, nI, nlines, nlines_z, nz, nresp, batch_count, response_mode, &
      !$omp& has_left_interface, has_right_interface, left_count, exposed_n, interior_base) &
      !$omp private(iline, ix_local, iz, abs_iz, resp_index, iface, slot, local_idx, actual_p, &
      !$omp& row_coeffs, rhs4, s4, col, coeff, coupled_row, j, exposed_slot, ghost_col, &
      !$omp& k, m, rhs_col, pivot4, factor4)
      do iline = 1, nlines
        ix_local = (iline - 1)/nlines_z
        iz = mod(iline - 1, nlines_z) - nz
        select case (response_mode)
        case (YS_ENDPOINT_RESPONSE_CONST)
          resp_index = 1_C_INT
        case default
          abs_iz = abs(iz)
          resp_index = ix_local*(nz + 1) + abs_iz + 1
        end select

        ys_reduced_rows_send(:, iline) = (0.0d0, 0.0d0)
        s4(:, :) = (0.0d0, 0.0d0)
        rhs4(:, :) = (0.0d0, 0.0d0)

        do iface = 1, exposed_n
          slot = iface
          if (.not. has_left_interface) slot = iface + 2
          local_idx = merge(slot - 1, active_n + slot - 5, slot <= 2)
          actual_p = local_idx*nlines + iline
          row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), &
                         ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)
          rhs4(iface, 1) = ys_gpsv_x(actual_p)

          do col = -2, 2
            coeff = row_coeffs(col)
            if (coeff == (0.0d0, 0.0d0)) cycle
            coupled_row = local_idx + col
            if (coupled_row >= interior_base .and. coupled_row <= interior_base + nI - 1) then
              j = coupled_row - interior_base
              rhs4(iface, 1) = rhs4(iface, 1) - coeff*ys_batch_x(j*batch_count + iline)
              do exposed_slot = 1, exposed_n
                s4(iface, exposed_slot) = s4(iface, exposed_slot) - &
                                          coeff*ys_batch_x(j*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)
              end do
            else
              exposed_slot = 0_C_INT
              if (has_left_interface) then
                if (coupled_row == 0) exposed_slot = 1_C_INT
                if (coupled_row == 1) exposed_slot = 2_C_INT
              end if
              if (has_right_interface) then
                if (coupled_row == active_n - 2) exposed_slot = left_count + 1
                if (coupled_row == active_n - 1) exposed_slot = left_count + 2
              end if
              if (exposed_slot > 0) then
                s4(iface, exposed_slot) = s4(iface, exposed_slot) + coeff
              else
                ghost_col = 0_C_INT
                if (coupled_row == -2) ghost_col = 1_C_INT
                if (coupled_row == -1) ghost_col = 2_C_INT
                if (coupled_row == active_n) ghost_col = 3_C_INT
                if (coupled_row == active_n + 1) ghost_col = 4_C_INT
                if (ghost_col > 0) rhs4(iface, ghost_col + 1) = rhs4(iface, ghost_col + 1) + coeff
              end if
            end if
          end do
        end do

        do k = 1, exposed_n
          pivot4 = s4(k, k)
          do m = k + 1, exposed_n
            factor4 = s4(m, k)/pivot4
            s4(m, k) = factor4
            s4(m, k + 1:exposed_n) = s4(m, k + 1:exposed_n) - factor4*s4(k, k + 1:exposed_n)
            rhs4(m, 1:5) = rhs4(m, 1:5) - factor4*rhs4(k, 1:5)
          end do
        end do
        do rhs_col = 1, 5
          do k = exposed_n, 1, -1
            if (k < exposed_n) rhs4(k, rhs_col) = rhs4(k, rhs_col) - &
                                                  sum(s4(k, k + 1:exposed_n)*rhs4(k + 1:exposed_n, rhs_col))
            rhs4(k, rhs_col) = rhs4(k, rhs_col)/s4(k, k)
          end do
        end do

        do iface = 1, exposed_n
          slot = iface
          if (.not. has_left_interface) slot = iface + 2
          ys_reduced_rows_send(slot, iline) = rhs4(iface, 1)
          do ghost_col = 1, 4
            ys_reduced_rows_send(4*ghost_col + slot, iline) = -rhs4(iface, ghost_col + 1)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
      call roctxPop("ys_endpoint_pack_schur")

      call ys_solve_reduced_interfaces()
    end if

    call roctxPush("ys_endpoint_reconstruct")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(dst, ys_batch_x, ys_reduced_rhs, ys_left_interface_values, ys_right_interface_values, active_n, nI, nlines, &
    !$omp& nlines_z, nz, nresp, batch_count, response_mode, dst_row_base, has_padded_dst, &
    !$omp& has_left_interface, has_right_interface, left_count, exposed_n, interior_base, ipy) &
    !$omp private(iline, ix_local, iz_index, iz, abs_iz, resp_index, iface_value, local_i, local_idx, p, exposed_slot)
    do iline = 1, nlines
      ix_local = (iline - 1)/nlines_z
      iz_index = mod(iline - 1, nlines_z) + 1
      iz = iz_index - nz - 1
      select case (response_mode)
      case (YS_ENDPOINT_RESPONSE_CONST)
        resp_index = 1_C_INT
      case default
        abs_iz = abs(iz)
        resp_index = ix_local*(nz + 1) + abs_iz + 1
      end select

      iface_value(:) = (0.0d0, 0.0d0)
      if (has_left_interface) then
        iface_value(1) = ys_reduced_rhs(4*ipy + 1, iline)
        iface_value(2) = ys_reduced_rhs(4*ipy + 2, iline)
        dst(dst_row_base, iz_index, ix_local + 1) = iface_value(1)
        dst(dst_row_base + 1, iz_index, ix_local + 1) = iface_value(2)
        if (has_padded_dst) dst(1:2, iz_index, ix_local + 1) = ys_left_interface_values(:, iline)
      end if
      if (has_right_interface) then
        iface_value(left_count + 1) = ys_reduced_rhs(4*ipy + 3, iline)
        iface_value(left_count + 2) = ys_reduced_rhs(4*ipy + 4, iline)
        dst(active_n - 2 + dst_row_base, iz_index, ix_local + 1) = iface_value(left_count + 1)
        dst(active_n - 1 + dst_row_base, iz_index, ix_local + 1) = iface_value(left_count + 2)
        if (has_padded_dst) dst(active_n + 3:active_n + 4, iz_index, ix_local + 1) = ys_right_interface_values(:, iline)
      end if

      do local_i = 0, nI - 1
        local_idx = interior_base + local_i
        p = local_i*batch_count + iline
        dst(local_idx + dst_row_base, iz_index, ix_local + 1) = ys_batch_x(p)
        do exposed_slot = 1, exposed_n
          dst(local_idx + dst_row_base, iz_index, ix_local + 1) = &
            dst(local_idx + dst_row_base, iz_index, ix_local + 1) - &
            ys_batch_x(local_i*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)*iface_value(exposed_slot)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_endpoint_reconstruct")
  end subroutine ys_solve_endpoint_schur

  subroutine ys_solve_reduced_interfaces()
    implicit none
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    complex(C_DOUBLE_COMPLEX) :: packed_remote(20)
    integer(C_INT) :: nlines, niface, iline, iblock, row0

    nlines = size(ys_reduced_rows_send, 2)
    niface = 4*npy_grid

    call allgather_y_device_complex_rows(ys_reduced_rows_send, ys_reduced_rows_recv, 20_C_INT, nlines, &
                                         "MPI_Allgather reduced_y_interfaces")

    call roctxPush("ys_reduced_interfaces_solve")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs, ys_left_interface_values, &
    !$omp& ys_right_interface_values, nlines, niface, npy_grid, ipy) &
    !$omp private(iline, iblock, row0, packed_remote)
    do iline = 1, nlines
      ys_reduced_matrix_lu(:, :, iline) = (0.0d0, 0.0d0)
      ys_reduced_rhs(:, iline) = (0.0d0, 0.0d0)

      do iblock = 0, npy_grid - 1
        row0 = 4*iblock
        packed_remote = ys_reduced_rows_recv(:, iline, iblock + 1)
        ys_reduced_matrix_lu(row0 + 1:row0 + 4, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_rhs(row0 + 1:row0 + 4, iline) = packed_remote(1:4)
      end do

      do iblock = 1, npy_grid - 1
        row0 = 4*iblock
        packed_remote = ys_reduced_rows_recv(:, iline, iblock + 1)
        ys_reduced_matrix_lu(row0 + 1, bw - 1, iline) = -packed_remote(5)
        ys_reduced_matrix_lu(row0 + 2, bw - 2, iline) = -packed_remote(6)
        ys_reduced_matrix_lu(row0 + 3, bw - 3, iline) = -packed_remote(7)
        ys_reduced_matrix_lu(row0 + 4, bw - 4, iline) = -packed_remote(8)
        ys_reduced_matrix_lu(row0 + 1, bw, iline) = -packed_remote(9)
        ys_reduced_matrix_lu(row0 + 2, bw - 1, iline) = -packed_remote(10)
        ys_reduced_matrix_lu(row0 + 3, bw - 2, iline) = -packed_remote(11)
        ys_reduced_matrix_lu(row0 + 4, bw - 3, iline) = -packed_remote(12)
      end do

      do iblock = 0, npy_grid - 2
        row0 = 4*iblock
        packed_remote = ys_reduced_rows_recv(:, iline, iblock + 1)
        ys_reduced_matrix_lu(row0 + 1, bw + 5, iline) = -packed_remote(13)
        ys_reduced_matrix_lu(row0 + 2, bw + 4, iline) = -packed_remote(14)
        ys_reduced_matrix_lu(row0 + 3, bw + 3, iline) = -packed_remote(15)
        ys_reduced_matrix_lu(row0 + 4, bw + 2, iline) = -packed_remote(16)
        ys_reduced_matrix_lu(row0 + 1, bw + 6, iline) = -packed_remote(17)
        ys_reduced_matrix_lu(row0 + 2, bw + 5, iline) = -packed_remote(18)
        ys_reduced_matrix_lu(row0 + 3, bw + 4, iline) = -packed_remote(19)
        ys_reduced_matrix_lu(row0 + 4, bw + 3, iline) = -packed_remote(20)
      end do

      call ys_factor_banded_complex(ys_reduced_matrix_lu(:, :, iline))
      call ys_solve_factored_banded_complex(ys_reduced_rhs(:, iline), ys_reduced_matrix_lu(:, :, iline))

      ys_left_interface_values(:, iline) = (0.0d0, 0.0d0)
      ys_right_interface_values(:, iline) = (0.0d0, 0.0d0)
      row0 = 4*ipy
      if (ipy > 0) ys_left_interface_values(:, iline) = ys_reduced_rhs(row0 - 1:row0, iline)
      if (ipy < npy_grid - 1) ys_right_interface_values(:, iline) = ys_reduced_rhs(row0 + 5:row0 + 6, iline)
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_reduced_interfaces_solve")
  end subroutine ys_solve_reduced_interfaces

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_factor_banded_complex)
#endif
  subroutine ys_factor_banded_complex(a)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: a(:, :)
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: n, i, j, last
    complex(C_DOUBLE_COMPLEX) :: piv, factor

    n = size(a, 1)
    do i = 1, n
      piv = a(i, bw + 1)
      last = min(bw, n - i)
      do j = 1, last
        factor = a(i + j, bw + 1 - j)/piv
        a(i + j, bw + 1 - j) = factor
        a(i + j, bw + 2 - j:bw + 1 + last - j) = &
          a(i + j, bw + 2 - j:bw + 1 + last - j) - factor*a(i, bw + 2:bw + 1 + last)
      end do
    end do
  end subroutine ys_factor_banded_complex

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_solve_factored_banded_complex)
#endif
  subroutine ys_solve_factored_banded_complex(rhs, a)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: a(:, :)
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: n, i, lo, hi

    n = size(a, 1)
    if (n <= 0) return

    do i = 2, n
      lo = max(1_C_INT, i - bw)
      rhs(i) = rhs(i) - sum(a(i, bw + 1 + lo - i:bw)*rhs(lo:i - 1))
    end do

    rhs(n) = rhs(n)/a(n, bw + 1)
    do i = n - 1, 1, -1
      hi = min(n, i + bw)
      rhs(i) = rhs(i) - sum(a(i, bw + 2:bw + 1 + hi - i)*rhs(i + 1:hi))
      rhs(i) = rhs(i)/a(i, bw + 1)
    end do
  end subroutine ys_solve_factored_banded_complex

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_factor_penta_interleaved)
#endif
  subroutine ys_factor_penta_interleaved(ds, dl, d, du, dw, stride, first, n)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i, p, p1, p2
    complex(C_DOUBLE_COMPLEX) :: factor

    if (n <= 0) return
    if (n == 1) then
      d(first) = 1.0d0/d(first)
      return
    end if

    do i = 0, n - 3
      p = first + i*stride
      d(p) = 1.0d0/d(p)

      p1 = p + stride
      factor = dl(p1)*d(p)
      dl(p1) = factor
      d(p1) = d(p1) - factor*du(p)
      du(p1) = du(p1) - factor*dw(p)

      p2 = p + 2*stride
      factor = ds(p2)*d(p)
      ds(p2) = factor
      dl(p2) = dl(p2) - factor*du(p)
      d(p2) = d(p2) - factor*dw(p)
    end do

    p = first + (n - 2)*stride
    d(p) = 1.0d0/d(p)
    p1 = p + stride
    factor = dl(p1)*d(p)
    dl(p1) = factor
    d(p1) = d(p1) - factor*du(p)

    d(p1) = 1.0d0/d(p1)
  end subroutine ys_factor_penta_interleaved

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_solve_factored_penta_interleaved)
#endif
  subroutine ys_solve_factored_penta_interleaved(rhs, ds, dl, d, du, dw, stride, first, n)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i, p

    if (n <= 0) return

    if (n >= 2) then
      p = first + stride
      rhs(p) = rhs(p) - dl(p)*rhs(first)
    end if
    do i = 2, n - 1
      p = first + i*stride
      rhs(p) = rhs(p) - ds(p)*rhs(first + (i - 2)*stride) - dl(p)*rhs(first + (i - 1)*stride)
    end do

    p = first + (n - 1)*stride
    rhs(p) = rhs(p)*d(p)
    if (n >= 2) then
      p = first + (n - 2)*stride
      rhs(p) = (rhs(p) - du(p)*rhs(p + stride))*d(p)
    end if
    do i = n - 3, 0, -1
      p = first + i*stride
      rhs(p) = (rhs(p) - du(p)*rhs(p + stride) - dw(p)*rhs(p + 2*stride))*d(p)
    end do
  end subroutine ys_solve_factored_penta_interleaved

end module y_line_solvers
