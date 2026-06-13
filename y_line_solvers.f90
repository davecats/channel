#include "header.h"

module y_line_solvers

  use, intrinsic :: iso_c_binding
  use mpi_transpose, only: ny0, nyN, nx0, nxN, nxB, npy_grid, ipy
#ifdef HAVE_MPI
  use mpi_transpose, only: MPI_COMM_Y, ensure_ycomm_buffers, ycomm_sendbuf, ycomm_recvbuf
#endif
  use roctx, only: roctxPush, roctxPop
#ifdef HAVE_MPI
  use mpi_f08
#endif
#ifdef HAVE_CUDA
  use cusparse
#elif defined(HAVE_HIP)
  use omp_lib
  use hipfort_hipsparse
#endif

  implicit none
  private

#ifdef HAVE_HIP
  interface hipsparseZgpsvInterleavedBatch_bufferSizeExt
    function hipsparseZgpsvInterleavedBatch_bufferSizeExt_(handle, algo, m, ds, dl, d, du, dw, x, &
                                                           batch_count, pBufferSizeInBytes) &
      bind(c, name="hipsparseZgpsvInterleavedBatch_bufferSizeExt")
      use, intrinsic :: iso_c_binding
      use hipfort_hipsparse_enums
      implicit none
      integer(kind(HIPSPARSE_STATUS_SUCCESS)) :: hipsparseZgpsvInterleavedBatch_bufferSizeExt_
      type(c_ptr), value :: handle
      integer(c_int), value :: algo, m, batch_count
      type(c_ptr), value :: ds, dl, d, du, dw, x
      integer(c_size_t) :: pBufferSizeInBytes
    end function hipsparseZgpsvInterleavedBatch_bufferSizeExt_
  end interface

  interface hipsparseZgpsvInterleavedBatch
    function hipsparseZgpsvInterleavedBatch_(handle, algo, m, ds, dl, d, du, dw, x, batch_count, pBuffer) &
      bind(c, name="hipsparseZgpsvInterleavedBatch")
      use, intrinsic :: iso_c_binding
      use hipfort_hipsparse_enums
      implicit none
      integer(kind(HIPSPARSE_STATUS_SUCCESS)) :: hipsparseZgpsvInterleavedBatch_
      type(c_ptr), value :: handle
      integer(c_int), value :: algo, m, batch_count
      type(c_ptr), value :: ds, dl, d, du, dw, x, pBuffer
    end function hipsparseZgpsvInterleavedBatch_
  end interface

  interface
    function hipDeviceSynchronize() bind(c, name="hipDeviceSynchronize")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: hipDeviceSynchronize
    end function hipDeviceSynchronize
  end interface

#endif

  integer(C_INT), parameter :: YS_ENDPOINT_RESPONSE_CONST = 1_C_INT
  integer(C_INT), parameter :: YS_ENDPOINT_RESPONSE_EVEN_Z = 2_C_INT
  integer(C_INT), parameter :: YS_REDUCED_BW = 5_C_INT
  integer(C_INT), parameter :: YS_REDUCED_RETURN_WIDTH = 8_C_INT
  integer(C_INT), parameter :: YS_HIER_EXCHANGE_ALLGATHER = 1_C_INT
  integer(C_INT), parameter :: YS_HIER_EXCHANGE_ALLTOALL = 2_C_INT
  integer(C_INT), parameter :: YS_LEAF_MODE_COMPOSE = 1_C_INT
  integer(C_INT), parameter :: YS_LEAF_MODE_NODE_SOURCE = 2_C_INT

  public :: ys_prepare_assembled_workspace, ys_release_workspace
  public :: ys_set_node_leaf_source_mode, ys_get_node_leaf_layout
  public :: ys_node_leaf_transpose_source_to_group, ys_node_leaf_transpose_solution_from_group
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
  integer(C_INT), save :: ys_workspace_reduced_node_size = -1
  integer(C_INT), save :: ys_workspace_reduced_leaf_mode = YS_LEAF_MODE_COMPOSE
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
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_matrix_lu(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rhs(:, :)
  integer(C_INT), save :: ys_hier_node_size = 1
  integer(C_INT), save :: ys_hier_ngroups = 1
  integer(C_INT), save :: ys_hier_group_id = 0
  integer(C_INT), save :: ys_hier_group_first = 0
  integer(C_INT), save :: ys_hier_group_count = 1
  integer(C_INT), save :: ys_hier_group_rank = 0
  integer(C_INT), save :: ys_hier_column_rank = 0
  integer(C_INT), save :: ys_hier_column_count = 1
  integer(C_INT), save :: ys_hier_node_first_line = 1
  integer(C_INT), save :: ys_hier_node_owned_nlines = 0
  integer(C_INT), save :: ys_hier_global_first_line = 1
  integer(C_INT), save :: ys_hier_global_owned_nlines = 0
  integer(C_INT), save :: ys_hier_column_solve_first_line = 1
  integer(C_INT), save :: ys_hier_column_solve_nlines = 0
  integer(C_INT), save :: ys_hier_column_send_elems = 0
  integer(C_INT), save :: ys_hier_column_recv_elems = 0
  integer(C_INT), save :: ys_hier_column_a2a_send_elems = 0
  integer(C_INT), save :: ys_hier_column_a2a_recv_elems = 0
  integer(C_INT), save :: ys_hier_column_return_send_elems = 0
  integer(C_INT), save :: ys_hier_column_return_recv_elems = 0
  integer(C_INT), save :: ys_hier_global_exchange = YS_HIER_EXCHANGE_ALLGATHER
  integer(C_INT), save :: ys_reduced_leaf_mode = YS_LEAF_MODE_COMPOSE
  logical, save :: ys_hier_use_column_global = .false.
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_hier_node_rows(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_hier_recover_basis(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_hier_group_values(:, :)
  integer(C_INT), allocatable, save :: ys_hier_group_first_by_group(:), ys_hier_group_count_by_group(:)
  integer(C_INT), allocatable, save :: ys_hier_rank_group(:), ys_hier_rank_group_pos(:)
  integer(C_INT), allocatable, save :: ys_hier_rank_line_first(:), ys_hier_rank_line_count(:)
  integer(C_INT), allocatable, save :: ys_hier_global_line_first(:), ys_hier_global_line_count(:)
  integer, allocatable, save :: ys_hier_leaf_send_counts(:), ys_hier_leaf_send_displs(:)
  integer, allocatable, save :: ys_hier_leaf_recv_counts(:), ys_hier_leaf_recv_displs(:)
  integer, allocatable, save :: ys_hier_global_send_counts(:), ys_hier_global_send_displs(:)
  integer, allocatable, save :: ys_hier_global_recv_counts(:), ys_hier_global_recv_displs(:)
  integer, allocatable, save :: ys_hier_column_global_send_counts(:), ys_hier_column_global_send_displs(:)
  integer, allocatable, save :: ys_hier_column_global_recv_counts(:), ys_hier_column_global_recv_displs(:)
  integer, allocatable, save :: ys_hier_column_return_send_counts(:), ys_hier_column_return_send_displs(:)
  integer, allocatable, save :: ys_hier_column_return_recv_counts(:), ys_hier_column_return_recv_displs(:)
  integer, allocatable, save :: ys_hier_group_return_send_counts(:), ys_hier_group_return_send_displs(:)
  integer, allocatable, save :: ys_hier_group_return_recv_counts(:), ys_hier_group_return_recv_displs(:)
  integer, allocatable, save :: ys_hier_leaf_return_send_counts(:), ys_hier_leaf_return_send_displs(:)
  integer, allocatable, save :: ys_hier_leaf_return_recv_counts(:), ys_hier_leaf_return_recv_displs(:)
  integer(C_INT), save :: ys_hier_leaf_send_elems = 0, ys_hier_leaf_recv_elems = 0
  integer(C_INT), save :: ys_hier_global_send_elems = 0, ys_hier_global_recv_elems = 0
  integer(C_INT), save :: ys_hier_group_return_send_elems = 0, ys_hier_group_return_recv_elems = 0
  integer(C_INT), save :: ys_hier_leaf_return_send_elems = 0, ys_hier_leaf_return_recv_elems = 0
  logical, save :: ys_hier_comm_stats_enabled = .false.
#ifdef HAVE_MPI
  type(MPI_Comm), save :: ys_hier_group_comm = MPI_COMM_NULL
  type(MPI_Comm), save :: ys_hier_column_comm = MPI_COMM_NULL
  logical, save :: ys_hier_group_comm_active = .false.
  logical, save :: ys_hier_column_comm_active = .false.
#endif

#ifdef HAVE_CUDA
  type(cusparseHandle), save :: ys_gpsv_handle
  logical, save :: ys_gpsv_handle_created = .false.
#elif defined(HAVE_HIP)
  type(c_ptr), save :: ys_gpsv_handle = c_null_ptr
  type(c_ptr), save :: ys_gpsv_buffer = c_null_ptr
  logical, save :: ys_gpsv_handle_created = .false.
#endif
  integer(C_INT), save :: ys_gpsv_n = -1, ys_gpsv_batch = -1
  integer(C_INT), save :: ys_batch_n = -1, ys_batch_count = -1
#if defined(HAVE_CUDA)
  integer(8), save :: ys_gpsv_buffer_size = 0_8
#elif defined(HAVE_HIP)
  integer(C_SIZE_T), save :: ys_gpsv_buffer_size = 0_C_SIZE_T
#endif
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_batch_ds(:), ys_batch_dl(:), ys_batch_d(:)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_batch_du(:), ys_batch_dw(:), ys_batch_x(:)
#if defined(HAVE_CUDA)
  character(c_char), allocatable, target, save :: ys_gpsv_buffer(:)
  !$omp declare target(ys_factor_banded_complex)
  !$omp declare target(ys_solve_factored_banded_complex)
  !$omp declare target(ys_solve_factored_banded_complex_multi)
  !$omp declare target(ys_factor_penta_interleaved)
  !$omp declare target(ys_solve_factored_penta_interleaved)
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
    ys_workspace_reduced_node_size = -1
    ys_workspace_reduced_leaf_mode = YS_LEAF_MODE_COMPOSE
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
              ys_workspace_npy == wanted_npy .and. ys_workspace_line_start == line_start .and. &
              ys_workspace_reduced_leaf_mode == ys_reduced_leaf_mode
  end function ys_core_shape_matches

  subroutine ys_set_node_leaf_source_mode(enabled)
    implicit none
    logical, intent(in) :: enabled

    ys_reduced_leaf_mode = merge(YS_LEAF_MODE_NODE_SOURCE, YS_LEAF_MODE_COMPOSE, enabled)
  end subroutine ys_set_node_leaf_source_mode

  subroutine ys_get_node_leaf_layout(ny, nlines, nlines_z, node_row_start, node_row_end, node_first_line, node_line_count)
    implicit none
    integer(C_INT), intent(in) :: ny, nlines, nlines_z
    integer(C_INT), intent(out) :: node_row_start, node_row_end, node_first_line, node_line_count
    integer(C_INT) :: group_last, first_x, x_count, nx_lines

    call ys_configure_reduced_solver(npy_grid)
    group_last = ys_hier_group_first + ys_hier_group_count - 1_C_INT
    node_row_start = 1_C_INT + ys_hier_group_first*(ny - 1_C_INT)/npy_grid
    node_row_end = (group_last + 1_C_INT)*(ny - 1_C_INT)/npy_grid
    nx_lines = nlines/nlines_z
    if (.not. ys_hier_use_column_global .or. nx_lines < ys_hier_group_count) then
      node_first_line = 1_C_INT
      node_line_count = -1_C_INT
      return
    end if
    call ys_split_line_range(ys_hier_group_rank, nx_lines, ys_hier_group_count, first_x, x_count)
    node_first_line = (first_x - 1_C_INT)*nlines_z + 1_C_INT
    node_line_count = x_count*nlines_z
  end subroutine ys_get_node_leaf_layout

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
#ifdef HAVE_MPI
    integer :: ierr_local
#endif

    if (allocated(ys_reduced_rows_send)) then
      !$omp target exit data map(delete: ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, &
      !$omp& ys_reduced_matrix_lu, ys_reduced_rhs)
      deallocate (ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values)
      deallocate (ys_reduced_matrix_lu, ys_reduced_rhs)
    end if
    if (allocated(ys_hier_node_rows)) then
      !$omp target exit data map(delete: ys_hier_node_rows, ys_hier_recover_basis, ys_hier_group_values)
      deallocate (ys_hier_node_rows, ys_hier_recover_basis, ys_hier_group_values)
    end if
    if (allocated(ys_hier_group_first_by_group)) then
      !$omp target exit data map(delete: ys_hier_group_first_by_group, ys_hier_group_count_by_group, &
      !$omp& ys_hier_rank_group, ys_hier_rank_group_pos, ys_hier_rank_line_first, ys_hier_rank_line_count, &
      !$omp& ys_hier_global_line_first, ys_hier_global_line_count, ys_hier_leaf_send_counts, &
      !$omp& ys_hier_leaf_send_displs, ys_hier_leaf_recv_counts, ys_hier_leaf_recv_displs, &
      !$omp& ys_hier_global_send_counts, ys_hier_global_send_displs, ys_hier_global_recv_counts, &
      !$omp& ys_hier_global_recv_displs, ys_hier_column_global_send_counts, &
      !$omp& ys_hier_column_global_send_displs, ys_hier_column_global_recv_counts, &
      !$omp& ys_hier_column_global_recv_displs, ys_hier_column_return_send_counts, &
      !$omp& ys_hier_column_return_send_displs, ys_hier_column_return_recv_counts, &
      !$omp& ys_hier_column_return_recv_displs, ys_hier_group_return_send_counts, &
      !$omp& ys_hier_group_return_send_displs, ys_hier_group_return_recv_counts, &
      !$omp& ys_hier_group_return_recv_displs, ys_hier_leaf_return_send_counts, &
      !$omp& ys_hier_leaf_return_send_displs, ys_hier_leaf_return_recv_counts, ys_hier_leaf_return_recv_displs)
      deallocate (ys_hier_group_first_by_group, ys_hier_group_count_by_group, ys_hier_rank_group, &
                  ys_hier_rank_group_pos, ys_hier_rank_line_first, ys_hier_rank_line_count, &
                  ys_hier_global_line_first, ys_hier_global_line_count)
      deallocate (ys_hier_leaf_send_counts, ys_hier_leaf_send_displs, ys_hier_leaf_recv_counts, ys_hier_leaf_recv_displs)
      deallocate (ys_hier_global_send_counts, ys_hier_global_send_displs, &
                  ys_hier_global_recv_counts, ys_hier_global_recv_displs)
      deallocate (ys_hier_column_global_send_counts, ys_hier_column_global_send_displs, &
                  ys_hier_column_global_recv_counts, ys_hier_column_global_recv_displs)
      deallocate (ys_hier_column_return_send_counts, ys_hier_column_return_send_displs, &
                  ys_hier_column_return_recv_counts, ys_hier_column_return_recv_displs)
      deallocate (ys_hier_group_return_send_counts, ys_hier_group_return_send_displs, &
                  ys_hier_group_return_recv_counts, ys_hier_group_return_recv_displs)
      deallocate (ys_hier_leaf_return_send_counts, ys_hier_leaf_return_send_displs, &
                  ys_hier_leaf_return_recv_counts, ys_hier_leaf_return_recv_displs)
    end if
    ys_hier_leaf_send_elems = 0
    ys_hier_leaf_recv_elems = 0
    ys_hier_global_send_elems = 0
    ys_hier_global_recv_elems = 0
    ys_hier_group_return_send_elems = 0
    ys_hier_group_return_recv_elems = 0
    ys_hier_leaf_return_send_elems = 0
    ys_hier_leaf_return_recv_elems = 0
    ys_hier_column_send_elems = 0
    ys_hier_column_recv_elems = 0
    ys_hier_column_a2a_send_elems = 0
    ys_hier_column_a2a_recv_elems = 0
    ys_hier_column_return_send_elems = 0
    ys_hier_column_return_recv_elems = 0
    ys_workspace_reduced_node_size = -1
    ys_workspace_reduced_leaf_mode = YS_LEAF_MODE_COMPOSE
#ifdef HAVE_MPI
    if (ys_hier_group_comm_active) then
      call MPI_Comm_free(ys_hier_group_comm, ierr_local)
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_free y-Schur group communicator failed"
      ys_hier_group_comm = MPI_COMM_NULL
      ys_hier_group_comm_active = .false.
    end if
    if (ys_hier_column_comm_active) then
      call MPI_Comm_free(ys_hier_column_comm, ierr_local)
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_free y-Schur column communicator failed"
      ys_hier_column_comm = MPI_COMM_NULL
      ys_hier_column_comm_active = .false.
    end if
#endif

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
    integer(C_INT) :: comm_send_elems, comm_recv_elems, node_alloc_nlines
#ifdef HAVE_MPI
    integer :: ierr_local
#endif

    call ys_configure_reduced_solver(npy_count)
#ifdef HAVE_MPI
    call MPI_Comm_split(MPI_COMM_Y, ys_hier_group_id, ys_hier_group_rank, ys_hier_group_comm, ierr_local)
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_split y-Schur group communicator failed"
    ys_hier_group_comm_active = .true.
    call MPI_Comm_split(MPI_COMM_Y, ys_hier_group_rank, ys_hier_group_id, ys_hier_column_comm, ierr_local)
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_split y-Schur column communicator failed"
    ys_hier_column_comm_active = .true.
    call MPI_Comm_rank(ys_hier_column_comm, ys_hier_column_rank, ierr_local)
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_rank y-Schur column communicator failed"
    call MPI_Comm_size(ys_hier_column_comm, ys_hier_column_count, ierr_local)
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_size y-Schur column communicator failed"
#endif
    allocate (ys_reduced_rows_send(20, nlines), ys_left_interface_values(2, nlines), &
              ys_right_interface_values(2, nlines))
    allocate (ys_reduced_matrix_lu(4*npy_count, 2*YS_REDUCED_BW + 1, nlines), &
              ys_reduced_rhs(4*npy_count, nlines))

    if (ys_reduced_leaf_mode == YS_LEAF_MODE_NODE_SOURCE) then
      ys_hier_node_first_line = ys_workspace_line_start
      ys_hier_node_owned_nlines = nlines
    else
      call ys_split_line_range(ys_hier_group_rank, nlines, ys_hier_group_count, &
                               ys_hier_node_first_line, ys_hier_node_owned_nlines)
    end if
    ys_hier_column_send_elems = 20*ys_hier_node_owned_nlines
    ys_hier_column_recv_elems = ys_hier_column_send_elems*ys_hier_column_count
    node_alloc_nlines = max(1_C_INT, ys_hier_node_owned_nlines)
    allocate (ys_hier_node_rows(20, node_alloc_nlines), &
              ys_hier_recover_basis(4*ys_hier_node_size, 5, node_alloc_nlines), &
              ys_hier_group_values(YS_REDUCED_RETURN_WIDTH, node_alloc_nlines))
    allocate (ys_hier_group_first_by_group(ys_hier_ngroups), ys_hier_group_count_by_group(ys_hier_ngroups))
    allocate (ys_hier_rank_group(npy_count), ys_hier_rank_group_pos(npy_count), &
              ys_hier_rank_line_first(npy_count), ys_hier_rank_line_count(npy_count), &
              ys_hier_global_line_first(npy_count), ys_hier_global_line_count(npy_count))
    allocate (ys_hier_leaf_send_counts(npy_count), ys_hier_leaf_send_displs(npy_count), &
              ys_hier_leaf_recv_counts(npy_count), ys_hier_leaf_recv_displs(npy_count))
    allocate (ys_hier_global_send_counts(npy_count), ys_hier_global_send_displs(npy_count), &
              ys_hier_global_recv_counts(npy_count), ys_hier_global_recv_displs(npy_count))
    allocate (ys_hier_column_global_send_counts(ys_hier_column_count), &
              ys_hier_column_global_send_displs(ys_hier_column_count), &
              ys_hier_column_global_recv_counts(ys_hier_column_count), &
              ys_hier_column_global_recv_displs(ys_hier_column_count))
    allocate (ys_hier_column_return_send_counts(ys_hier_column_count), &
              ys_hier_column_return_send_displs(ys_hier_column_count), &
              ys_hier_column_return_recv_counts(ys_hier_column_count), &
              ys_hier_column_return_recv_displs(ys_hier_column_count))
    allocate (ys_hier_group_return_send_counts(npy_count), ys_hier_group_return_send_displs(npy_count), &
              ys_hier_group_return_recv_counts(npy_count), ys_hier_group_return_recv_displs(npy_count))
    allocate (ys_hier_leaf_return_send_counts(npy_count), ys_hier_leaf_return_send_displs(npy_count), &
              ys_hier_leaf_return_recv_counts(npy_count), ys_hier_leaf_return_recv_displs(npy_count))
    call ys_prepare_hier_comm_plan(nlines, npy_count)
    call ys_prepare_hier_column_comm_plan()

    comm_send_elems = max(ys_hier_leaf_send_elems, ys_hier_global_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_column_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_column_a2a_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_column_return_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_group_return_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_leaf_return_send_elems)
    comm_recv_elems = max(ys_hier_leaf_recv_elems, ys_hier_global_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_column_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_column_a2a_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_column_return_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_group_return_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_leaf_return_recv_elems)
    call ensure_ycomm_buffers(comm_send_elems, comm_recv_elems)

    !$omp target enter data map(alloc: ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, &
    !$omp& ys_reduced_matrix_lu, ys_reduced_rhs)
    !$omp target enter data map(alloc: ys_hier_node_rows, ys_hier_recover_basis, ys_hier_group_values)
    !$omp target enter data map(to: ys_hier_group_first_by_group, ys_hier_group_count_by_group, &
    !$omp& ys_hier_rank_group, ys_hier_rank_group_pos, ys_hier_rank_line_first, ys_hier_rank_line_count, &
    !$omp& ys_hier_global_line_first, ys_hier_global_line_count, ys_hier_leaf_send_counts, &
    !$omp& ys_hier_leaf_send_displs, ys_hier_leaf_recv_counts, ys_hier_leaf_recv_displs, &
    !$omp& ys_hier_global_send_counts, ys_hier_global_send_displs, ys_hier_global_recv_counts, &
    !$omp& ys_hier_global_recv_displs, ys_hier_column_global_send_counts, &
    !$omp& ys_hier_column_global_send_displs, ys_hier_column_global_recv_counts, &
    !$omp& ys_hier_column_global_recv_displs, ys_hier_column_return_send_counts, &
    !$omp& ys_hier_column_return_send_displs, ys_hier_column_return_recv_counts, &
    !$omp& ys_hier_column_return_recv_displs, ys_hier_group_return_send_counts, &
    !$omp& ys_hier_group_return_send_displs, ys_hier_group_return_recv_counts, &
    !$omp& ys_hier_group_return_recv_displs, ys_hier_leaf_return_send_counts, &
    !$omp& ys_hier_leaf_return_send_displs, ys_hier_leaf_return_recv_counts, ys_hier_leaf_return_recv_displs)
    ys_workspace_reduced_node_size = ys_hier_node_size
    ys_workspace_reduced_leaf_mode = ys_reduced_leaf_mode
  end subroutine ys_allocate_reduced_workspace

  subroutine ys_configure_reduced_solver(npy_count)
    implicit none
    integer(C_INT), intent(in) :: npy_count
    character(len=32) :: env_value
    character(len=32) :: env_switch
    integer :: env_length, env_status, ios, requested_node_size

    if (npy_count < 2) error stop "Hierarchical y-Schur requires at least two y ranks"
    ys_hier_node_size = min(4_C_INT, npy_count)
    call get_environment_variable("CHANNEL_Y_SCHUR_NODE_SIZE", env_value, length=env_length, status=env_status)
    if (env_status == 0 .and. env_length > 0) then
      read (env_value(:env_length), *, iostat=ios) requested_node_size
      if (ios /= 0) error stop "Invalid CHANNEL_Y_SCHUR_NODE_SIZE"
      ys_hier_node_size = int(requested_node_size, C_INT)
    end if
    if (ys_hier_node_size < 2 .or. ys_hier_node_size > npy_count) error stop "Invalid CHANNEL_Y_SCHUR_NODE_SIZE"

    ys_hier_global_exchange = YS_HIER_EXCHANGE_ALLGATHER
    call get_environment_variable("CHANNEL_Y_SCHUR_GLOBAL_EXCHANGE", env_value, length=env_length, status=env_status)
    if (env_status == 0 .and. env_length > 0) then
      env_switch = trim(adjustl(env_value(:env_length)))
      select case (env_switch)
      case ("allgather", "ALLGATHER", "gather", "GATHER")
        ys_hier_global_exchange = YS_HIER_EXCHANGE_ALLGATHER
      case ("alltoall", "ALLTOALL", "a2a", "A2A")
        ys_hier_global_exchange = YS_HIER_EXCHANGE_ALLTOALL
      case default
        error stop "Invalid CHANNEL_Y_SCHUR_GLOBAL_EXCHANGE"
      end select
    end if

    ys_hier_comm_stats_enabled = .false.
    call get_environment_variable("CHANNEL_Y_SCHUR_COMM_STATS", env_value, length=env_length, status=env_status)
    if (env_status == 0 .and. env_length > 0) then
      env_switch = trim(adjustl(env_value(:env_length)))
      select case (env_switch)
      case ("1", "true", "TRUE", "yes", "YES", "on", "ON")
        ys_hier_comm_stats_enabled = .true.
      case ("0", "false", "FALSE", "no", "NO", "off", "OFF")
        ys_hier_comm_stats_enabled = .false.
      case default
        error stop "Invalid CHANNEL_Y_SCHUR_COMM_STATS"
      end select
    end if

    ys_hier_ngroups = (npy_count + ys_hier_node_size - 1_C_INT)/ys_hier_node_size
    ys_hier_group_id = ipy/ys_hier_node_size
    ys_hier_group_first = ys_hier_group_id*ys_hier_node_size
    ys_hier_group_count = min(ys_hier_node_size, npy_count - ys_hier_group_first)
    ys_hier_group_rank = ipy - ys_hier_group_first
    ys_hier_use_column_global = (ys_hier_ngroups*ys_hier_node_size == npy_count)
  end subroutine ys_configure_reduced_solver

  subroutine ys_split_line_range(rank, nitems, nranks, first_item, item_count)
    implicit none
    integer(C_INT), intent(in) :: rank, nitems, nranks
    integer(C_INT), intent(out) :: first_item, item_count
    integer(C_INT) :: base_count, remainder

    base_count = nitems/nranks
    remainder = mod(nitems, nranks)
    item_count = base_count
    if (rank < remainder) item_count = item_count + 1
    first_item = rank*base_count + min(rank, remainder) + 1
  end subroutine ys_split_line_range

  integer(C_INT) function ys_range_intersection_count(first_a, count_a, first_b, count_b) result(count)
    implicit none
    integer(C_INT), intent(in) :: first_a, count_a, first_b, count_b
    integer(C_INT) :: lo, hi

    lo = max(first_a, first_b)
    hi = min(first_a + count_a - 1_C_INT, first_b + count_b - 1_C_INT)
    count = max(0_C_INT, hi - lo + 1_C_INT)
  end function ys_range_intersection_count

  subroutine ys_prepare_hier_comm_plan(nlines, npy_count)
    implicit none
    integer(C_INT), intent(in) :: nlines, npy_count
    integer(C_INT) :: rank, group_id, group_first, group_count, group_rank
    integer(C_INT) :: my_group_id, first_line, line_count, global_first, global_count
    integer(C_INT) :: send_offset, recv_offset, intersection_count

    my_group_id = ipy/ys_hier_node_size
    call ys_split_line_range(ipy, nlines, npy_count, ys_hier_global_first_line, ys_hier_global_owned_nlines)

    do group_id = 0, ys_hier_ngroups - 1
      group_first = group_id*ys_hier_node_size
      group_count = min(ys_hier_node_size, npy_count - group_first)
      ys_hier_group_first_by_group(group_id + 1) = group_first
      ys_hier_group_count_by_group(group_id + 1) = group_count
    end do

    do rank = 0, npy_count - 1
      group_id = rank/ys_hier_node_size
      group_first = group_id*ys_hier_node_size
      group_count = min(ys_hier_node_size, npy_count - group_first)
      group_rank = rank - group_first
      ys_hier_rank_group(rank + 1) = group_id
      ys_hier_rank_group_pos(rank + 1) = group_rank
      call ys_split_line_range(group_rank, nlines, group_count, first_line, line_count)
      ys_hier_rank_line_first(rank + 1) = first_line
      ys_hier_rank_line_count(rank + 1) = line_count
      call ys_split_line_range(rank, nlines, npy_count, first_line, line_count)
      ys_hier_global_line_first(rank + 1) = first_line
      ys_hier_global_line_count(rank + 1) = line_count
    end do

    send_offset = 0_C_INT
    recv_offset = 0_C_INT
    do rank = 0, npy_count - 1
      if (ys_hier_rank_group(rank + 1) == my_group_id) then
        ys_hier_leaf_send_counts(rank + 1) = 20*ys_hier_rank_line_count(rank + 1)
        ys_hier_leaf_recv_counts(rank + 1) = 20*ys_hier_node_owned_nlines
        ys_hier_leaf_return_send_counts(rank + 1) = YS_REDUCED_RETURN_WIDTH*ys_hier_node_owned_nlines
        ys_hier_leaf_return_recv_counts(rank + 1) = YS_REDUCED_RETURN_WIDTH*ys_hier_rank_line_count(rank + 1)
      else
        ys_hier_leaf_send_counts(rank + 1) = 0
        ys_hier_leaf_recv_counts(rank + 1) = 0
        ys_hier_leaf_return_send_counts(rank + 1) = 0
        ys_hier_leaf_return_recv_counts(rank + 1) = 0
      end if
      ys_hier_leaf_send_displs(rank + 1) = send_offset
      ys_hier_leaf_recv_displs(rank + 1) = recv_offset
      send_offset = send_offset + ys_hier_leaf_send_counts(rank + 1)
      recv_offset = recv_offset + ys_hier_leaf_recv_counts(rank + 1)
    end do
    ys_hier_leaf_send_elems = send_offset
    ys_hier_leaf_recv_elems = recv_offset

    send_offset = 0_C_INT
    recv_offset = 0_C_INT
    do rank = 0, npy_count - 1
      ys_hier_leaf_return_send_displs(rank + 1) = send_offset
      ys_hier_leaf_return_recv_displs(rank + 1) = recv_offset
      send_offset = send_offset + ys_hier_leaf_return_send_counts(rank + 1)
      recv_offset = recv_offset + ys_hier_leaf_return_recv_counts(rank + 1)
    end do
    ys_hier_leaf_return_send_elems = send_offset
    ys_hier_leaf_return_recv_elems = recv_offset

    send_offset = 0_C_INT
    recv_offset = 0_C_INT
    do rank = 0, npy_count - 1
      global_first = ys_hier_global_line_first(rank + 1)
      global_count = ys_hier_global_line_count(rank + 1)
      first_line = ys_hier_rank_line_first(rank + 1)
      line_count = ys_hier_rank_line_count(rank + 1)
      intersection_count = ys_range_intersection_count(ys_hier_node_first_line, ys_hier_node_owned_nlines, &
                                                       global_first, global_count)
      ys_hier_global_send_counts(rank + 1) = 20*intersection_count
      intersection_count = ys_range_intersection_count(ys_hier_global_first_line, ys_hier_global_owned_nlines, &
                                                       first_line, line_count)
      ys_hier_global_recv_counts(rank + 1) = 20*intersection_count
      ys_hier_global_send_displs(rank + 1) = send_offset
      ys_hier_global_recv_displs(rank + 1) = recv_offset
      send_offset = send_offset + ys_hier_global_send_counts(rank + 1)
      recv_offset = recv_offset + ys_hier_global_recv_counts(rank + 1)
    end do
    ys_hier_global_send_elems = send_offset
    ys_hier_global_recv_elems = recv_offset

    send_offset = 0_C_INT
    recv_offset = 0_C_INT
    do rank = 0, npy_count - 1
      first_line = ys_hier_rank_line_first(rank + 1)
      line_count = ys_hier_rank_line_count(rank + 1)
      global_first = ys_hier_global_line_first(rank + 1)
      global_count = ys_hier_global_line_count(rank + 1)
      intersection_count = ys_range_intersection_count(ys_hier_global_first_line, ys_hier_global_owned_nlines, &
                                                       first_line, line_count)
      ys_hier_group_return_send_counts(rank + 1) = YS_REDUCED_RETURN_WIDTH*intersection_count
      intersection_count = ys_range_intersection_count(ys_hier_node_first_line, ys_hier_node_owned_nlines, &
                                                       global_first, global_count)
      ys_hier_group_return_recv_counts(rank + 1) = YS_REDUCED_RETURN_WIDTH*intersection_count
      ys_hier_group_return_send_displs(rank + 1) = send_offset
      ys_hier_group_return_recv_displs(rank + 1) = recv_offset
      send_offset = send_offset + ys_hier_group_return_send_counts(rank + 1)
      recv_offset = recv_offset + ys_hier_group_return_recv_counts(rank + 1)
    end do
    ys_hier_group_return_send_elems = send_offset
    ys_hier_group_return_recv_elems = recv_offset
  end subroutine ys_prepare_hier_comm_plan

  subroutine ys_prepare_hier_column_comm_plan()
    implicit none
    integer(C_INT) :: rank, first_line, line_count
    integer(C_INT) :: send_offset, recv_offset

    call ys_split_line_range(ys_hier_column_rank, ys_hier_node_owned_nlines, ys_hier_column_count, &
                             ys_hier_column_solve_first_line, ys_hier_column_solve_nlines)

    send_offset = 0_C_INT
    recv_offset = 0_C_INT
    do rank = 0, ys_hier_column_count - 1
      call ys_split_line_range(rank, ys_hier_node_owned_nlines, ys_hier_column_count, first_line, line_count)
      ys_hier_column_global_send_counts(rank + 1) = 20*line_count
      ys_hier_column_global_recv_counts(rank + 1) = 20*ys_hier_column_solve_nlines
      ys_hier_column_global_send_displs(rank + 1) = send_offset
      ys_hier_column_global_recv_displs(rank + 1) = recv_offset
      send_offset = send_offset + ys_hier_column_global_send_counts(rank + 1)
      recv_offset = recv_offset + ys_hier_column_global_recv_counts(rank + 1)
    end do
    ys_hier_column_a2a_send_elems = send_offset
    ys_hier_column_a2a_recv_elems = recv_offset

    send_offset = 0_C_INT
    recv_offset = 0_C_INT
    do rank = 0, ys_hier_column_count - 1
      call ys_split_line_range(rank, ys_hier_node_owned_nlines, ys_hier_column_count, first_line, line_count)
      ys_hier_column_return_send_counts(rank + 1) = YS_REDUCED_RETURN_WIDTH*ys_hier_column_solve_nlines
      ys_hier_column_return_recv_counts(rank + 1) = YS_REDUCED_RETURN_WIDTH*line_count
      ys_hier_column_return_send_displs(rank + 1) = send_offset
      ys_hier_column_return_recv_displs(rank + 1) = recv_offset
      send_offset = send_offset + ys_hier_column_return_send_counts(rank + 1)
      recv_offset = recv_offset + ys_hier_column_return_recv_counts(rank + 1)
    end do
    ys_hier_column_return_send_elems = send_offset
    ys_hier_column_return_recv_elems = recv_offset
  end subroutine ys_prepare_hier_column_comm_plan

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
    if (use_reduced_backend .and. npy_grid > 1) then
      if (npy_grid > 1 .and. active_n < 4) error stop "ys_solve_ghost_field requires at least four y rows per rank"
      if (active_n < 4) error stop "ys_solve_ghost_field requires at least four active y rows"
      call ys_configure_reduced_solver(npy_grid)
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

    if (use_reduced_backend .and. npy_grid > 1) then
      if (allocated(ys_reduced_rows_send)) then
        if (ys_workspace_reduced_node_size /= ys_hier_node_size .or. &
            ys_workspace_reduced_leaf_mode /= ys_reduced_leaf_mode) call ys_release_reduced_workspace()
      end if
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

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  subroutine ys_check_gpusparse(status, where)
    integer(C_INT), intent(in) :: status
    character(*), intent(in) :: where

#ifdef HAVE_CUDA
    if (status /= CUSPARSE_STATUS_SUCCESS) then
      print *, "cuSPARSE error in ", trim(where), ": status=", status
      error stop
    end if
#elif defined(HAVE_HIP)
    if (status /= HIPSPARSE_STATUS_SUCCESS) then
      print *, "hipSPARSE error in ", trim(where), ": status=", status
      error stop
    end if
#endif
  end subroutine ys_check_gpusparse

#ifdef HAVE_HIP
  subroutine ys_check_hip(status, where)
    integer(C_INT), intent(in) :: status
    character(*), intent(in) :: where

    if (status /= 0_C_INT) then
      print *, "HIP runtime error in ", trim(where), ": status=", status
      error stop
    end if
  end subroutine ys_check_hip
#endif

  subroutine ys_create_gpusparse_handle()
    implicit none
    integer(C_INT) :: status

    if (ys_gpsv_handle_created) return

#ifdef HAVE_CUDA
    status = cusparseCreate(ys_gpsv_handle)
    call ys_check_gpusparse(status, "cusparseCreate")
#elif defined(HAVE_HIP)
    status = hipsparseCreate(ys_gpsv_handle)
    call ys_check_gpusparse(status, "hipsparseCreate")
#endif
    ys_gpsv_handle_created = .true.
  end subroutine ys_create_gpusparse_handle

  subroutine ys_destroy_gpusparse_handle()
    implicit none
    integer(C_INT) :: status

    if (.not. ys_gpsv_handle_created) return

#ifdef HAVE_CUDA
    status = cusparseDestroy(ys_gpsv_handle)
    call ys_check_gpusparse(status, "cusparseDestroy")
#elif defined(HAVE_HIP)
    status = hipsparseDestroy(ys_gpsv_handle)
    call ys_check_gpusparse(status, "hipsparseDestroy")
    ys_gpsv_handle = c_null_ptr
#endif
    ys_gpsv_handle_created = .false.
  end subroutine ys_destroy_gpusparse_handle

  subroutine ys_query_gpsv_buffer_size(ds, dl, d, du, dw, x, n, batch_count, buffer_size)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: ds(:), dl(:), d(:), du(:), dw(:), x(:)
    integer(C_INT), intent(in) :: n, batch_count
#ifdef HAVE_CUDA
    integer(8), intent(out) :: buffer_size
#elif defined(HAVE_HIP)
    integer(C_SIZE_T), intent(out) :: buffer_size
#endif
    integer(C_INT) :: status

    !$omp target data use_device_addr(ds, dl, d, du, dw, x)
#ifdef HAVE_CUDA
    status = cusparseZgpsvInterleavedBatch_bufferSize(ys_gpsv_handle, 0_C_INT, n, ds, dl, d, du, dw, x, &
                                                      batch_count, buffer_size)
#elif defined(HAVE_HIP)
    call ys_check_hip(hipDeviceSynchronize(), "hipDeviceSynchronize before gpsvInterleavedBatch_bufferSize")
    status = hipsparseZgpsvInterleavedBatch_bufferSizeExt(ys_gpsv_handle, 0_C_INT, n, c_loc(ds(1)), c_loc(dl(1)), &
                                                          c_loc(d(1)), c_loc(du(1)), c_loc(dw(1)), c_loc(x(1)), &
                                                          batch_count, buffer_size)
    call ys_check_hip(hipDeviceSynchronize(), "hipDeviceSynchronize after gpsvInterleavedBatch_bufferSize")
#endif
    !$omp end target data
    call ys_check_gpusparse(status, "gpsvInterleavedBatch_bufferSize")
  end subroutine ys_query_gpsv_buffer_size

  subroutine ys_call_gpsv_interleaved(ds, dl, d, du, dw, x, n, batch_count, label)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: ds(:), dl(:), d(:), du(:), dw(:), x(:)
    integer(C_INT), intent(in) :: n, batch_count
    character(*), intent(in) :: label
    integer(C_INT) :: status

#ifdef HAVE_CUDA
    !$omp target data use_device_addr(ds, dl, d, du, dw, x, ys_gpsv_buffer)
    status = cusparseZgpsvInterleavedBatch(ys_gpsv_handle, 0_C_INT, n, ds, dl, d, du, dw, x, &
                                           batch_count, ys_gpsv_buffer)
    call ys_check_gpusparse(status, "cusparseZgpsvInterleavedBatch "//trim(label))
    !$omp end target data
#elif defined(HAVE_HIP)
    !$omp target data use_device_addr(ds, dl, d, du, dw, x)
    call ys_check_hip(hipDeviceSynchronize(), "hipDeviceSynchronize before "//trim(label))
    status = hipsparseZgpsvInterleavedBatch(ys_gpsv_handle, 0_C_INT, n, c_loc(ds(1)), c_loc(dl(1)), c_loc(d(1)), &
                                            c_loc(du(1)), c_loc(dw(1)), c_loc(x(1)), batch_count, ys_gpsv_buffer)
    call ys_check_gpusparse(status, "hipsparseZgpsvInterleavedBatch "//trim(label))
    call ys_check_hip(hipDeviceSynchronize(), "hipDeviceSynchronize after "//trim(label))
    !$omp end target data
#endif
  end subroutine ys_call_gpsv_interleaved
#endif

  subroutine ys_release_gpsv_workspace()
    implicit none

    if (allocated(ys_batch_ds)) then
      !$omp target exit data map(delete: ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
      deallocate (ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
    end if
#ifdef HAVE_CUDA
    if (allocated(ys_gpsv_buffer)) then
      !$omp target exit data map(delete: ys_gpsv_buffer)
      deallocate (ys_gpsv_buffer)
    end if
#elif defined(HAVE_HIP)
    if (c_associated(ys_gpsv_buffer)) then
      call omp_target_free(ys_gpsv_buffer, omp_get_default_device())
      ys_gpsv_buffer = c_null_ptr
    end if
#endif
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    call ys_destroy_gpusparse_handle()
#endif
    ys_gpsv_n = -1
    ys_gpsv_batch = -1
    ys_batch_n = -1
    ys_batch_count = -1
#ifdef HAVE_CUDA
    ys_gpsv_buffer_size = 0_8
#elif defined(HAVE_HIP)
    ys_gpsv_buffer_size = 0_C_SIZE_T
    ys_gpsv_buffer = c_null_ptr
#endif
  end subroutine ys_release_gpsv_workspace

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  subroutine ys_prepare_gpusparse_workspace(n, batch_count, ds, dl, d, du, dw, x)
    implicit none
    integer(C_INT), intent(in) :: n, batch_count
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: ds(:), dl(:), d(:), du(:), dw(:), x(:)
#ifdef HAVE_CUDA
    integer(8) :: buffer_size
#elif defined(HAVE_HIP)
    integer(C_SIZE_T) :: buffer_size
#endif

    if (ys_gpsv_n == n .and. ys_gpsv_batch == batch_count &
#ifdef HAVE_CUDA
        .and. allocated(ys_gpsv_buffer) &
#elif defined(HAVE_HIP)
        .and. c_associated(ys_gpsv_buffer) &
#endif
        ) return

#ifdef HAVE_CUDA
    if (allocated(ys_gpsv_buffer)) then
      !$omp target exit data map(delete: ys_gpsv_buffer)
      deallocate (ys_gpsv_buffer)
    end if
#elif defined(HAVE_HIP)
    if (c_associated(ys_gpsv_buffer)) then
      call omp_target_free(ys_gpsv_buffer, omp_get_default_device())
      ys_gpsv_buffer = c_null_ptr
    end if
#endif

    call ys_create_gpusparse_handle()
    call ys_query_gpsv_buffer_size(ds, dl, d, du, dw, x, n, batch_count, buffer_size)

    ys_gpsv_buffer_size = buffer_size
#ifdef HAVE_CUDA
    allocate (ys_gpsv_buffer(max(1, int(buffer_size))))
    !$omp target enter data map(alloc: ys_gpsv_buffer)
#elif defined(HAVE_HIP)
    ! On MI300A, hipSPARSE gpsv rejects a workspace buffer created by mapping a
    ! Fortran character array with OpenMP target data, even though the operand
    ! arrays from use_device_addr() are accepted. omp_target_alloc() and hipMalloc()
    ! both work for this buffer; use the OpenMP allocator here to match mpi_transpose.
    ys_gpsv_buffer = omp_target_alloc(max(1_C_SIZE_T, buffer_size), omp_get_default_device())
    if (.not. c_associated(ys_gpsv_buffer)) then
      print *, "OpenMP target allocation failed in ys_prepare_gpusparse_workspace"
      error stop
    end if
#endif

    ys_gpsv_n = n
    ys_gpsv_batch = batch_count
  end subroutine ys_prepare_gpusparse_workspace
#else
  subroutine ys_prepare_gpusparse_workspace(n, batch_count)
    implicit none
    integer(C_INT), intent(in) :: n, batch_count

    ys_gpsv_n = n
    ys_gpsv_batch = batch_count
  end subroutine ys_prepare_gpusparse_workspace
#endif

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

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    call ys_prepare_gpusparse_workspace(n, batch_count, ys_batch_ds, ys_batch_dl, ys_batch_d, &
                                        ys_batch_du, ys_batch_dw, ys_batch_x)
#else
    call ys_prepare_gpusparse_workspace(n, batch_count)
#endif
  end subroutine ys_prepare_batch_workspace

  subroutine ys_solve_interleaved_pentadiagonal(ds, dl, d, du, dw, x, n, batch_count, label)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: ds(:), dl(:), d(:), du(:), dw(:), x(:)
    integer(C_INT), intent(in) :: n, batch_count
    character(*), intent(in) :: label
    integer(C_INT) :: iline, p, p1
    complex(C_DOUBLE_COMPLEX) :: d0inv, d1inv, lfac, rhs0, rhs1

    call roctxPush(label)
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    ! ROCSparse fails if the system size is 1: https://github.com/ROCm/rocSPARSE/blob/develop_deprecated/library/src/precond/rocsparse_gtsv.cpp
    if (n < 3) then
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(ds, dl, d, du, dw, x, n, batch_count) private(iline, p, p1, d0inv, d1inv, lfac, rhs0, rhs1)
      do iline = 1, batch_count
        p = iline
        if (n == 1) then
          d0inv = 1.0d0/d(p)
          d(p) = d0inv
          x(p) = x(p)*d0inv
        else
          p1 = p + batch_count
          d0inv = 1.0d0/d(p)
          lfac = dl(p1)*d0inv
          d1inv = 1.0d0/(d(p1) - lfac*du(p))
          rhs1 = (x(p1) - lfac*x(p))*d1inv
          rhs0 = (x(p) - du(p)*rhs1)*d0inv

          d(p) = d0inv
          dl(p1) = lfac
          d(p1) = d1inv
          x(p) = rhs0
          x(p1) = rhs1
        end if
      end do
      !$omp end target teams distribute parallel do
    else
      call ys_prepare_gpusparse_workspace(n, batch_count, ds, dl, d, du, dw, x)
      call ys_call_gpsv_interleaved(ds, dl, d, du, dw, x, n, batch_count, label)
    end if
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
    integer(C_INT) :: reduced_leaf_id, reduced_leaf_count, reduced_row0
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

    if (ys_reduced_leaf_mode == YS_LEAF_MODE_NODE_SOURCE) then
      reduced_leaf_id = ys_hier_group_id
      reduced_leaf_count = ys_hier_ngroups
    else
      reduced_leaf_id = ipy
      reduced_leaf_count = npy_grid
    end if
    reduced_row0 = 4*reduced_leaf_id

    has_left_interface = (reduced_leaf_id > 0)
    has_right_interface = (reduced_leaf_id < reduced_leaf_count - 1)
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
    !$omp& has_left_interface, has_right_interface, left_count, exposed_n, interior_base, reduced_row0) &
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
        iface_value(1) = ys_reduced_rhs(reduced_row0 + 1, iline)
        iface_value(2) = ys_reduced_rhs(reduced_row0 + 2, iline)
        dst(dst_row_base, iz_index, ix_local + 1) = iface_value(1)
        dst(dst_row_base + 1, iz_index, ix_local + 1) = iface_value(2)
        if (has_padded_dst) dst(1:2, iz_index, ix_local + 1) = ys_left_interface_values(:, iline)
      end if
      if (has_right_interface) then
        iface_value(left_count + 1) = ys_reduced_rhs(reduced_row0 + 3, iline)
        iface_value(left_count + 2) = ys_reduced_rhs(reduced_row0 + 4, iline)
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

    if (ys_reduced_leaf_mode == YS_LEAF_MODE_NODE_SOURCE) then
      call ys_solve_reduced_interfaces_node_leaf()
    else
      call ys_solve_reduced_interfaces_hier()
    end if
  end subroutine ys_solve_reduced_interfaces

  subroutine ys_solve_reduced_interfaces_node_leaf()
    implicit none
    integer :: ierr_local
    integer(C_INT) :: nlines, comm_send_elems, comm_recv_elems
    real(C_DOUBLE) :: comm_t0 = 0.0_C_DOUBLE
    real(C_DOUBLE) :: comm_elapsed = 0.0_C_DOUBLE

    if (.not. ys_hier_use_column_global) error stop "node-source y-Schur requires even node groups"

    nlines = size(ys_reduced_rows_send, 2)
    comm_send_elems = max(ys_hier_column_send_elems, ys_hier_column_a2a_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_column_return_send_elems)
    comm_recv_elems = max(ys_hier_column_recv_elems, ys_hier_column_a2a_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_column_return_recv_elems)
    call ensure_ycomm_buffers(comm_send_elems, comm_recv_elems)

    select case (ys_hier_global_exchange)
    case (YS_HIER_EXCHANGE_ALLGATHER)
      call ys_node_leaf_pack_rows_column(nlines)
      call roctxPush("MPI_Allgather node_leaf_rows_to_global")
      if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
      !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
      call MPI_Allgather(ycomm_sendbuf(1:ys_hier_column_send_elems), ys_hier_column_send_elems, MPI_DOUBLE_COMPLEX, &
                         ycomm_recvbuf(1:ys_hier_column_recv_elems), ys_hier_column_send_elems, MPI_DOUBLE_COMPLEX, &
                         ys_hier_column_comm, ierr_local)
      !$omp end target data
      if (ys_hier_comm_stats_enabled) then
        comm_elapsed = MPI_Wtime() - comm_t0
        call ys_hier_report_comm_stats_elems("node_leaf_rows_to_global", ys_hier_column_comm, &
                                             ys_hier_column_send_elems, ys_hier_column_recv_elems, comm_elapsed)
      end if
      call roctxPop("MPI_Allgather node_leaf_rows_to_global")
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Allgather node_leaf_rows_to_global failed"

      call ys_node_leaf_solve_global_groups_column(nlines)

    case (YS_HIER_EXCHANGE_ALLTOALL)
      call ys_node_leaf_pack_rows_to_column_global(nlines)
      call roctxPush("MPI_Alltoallv node_leaf_rows_to_global")
      if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
      !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
      call MPI_Alltoallv(ycomm_sendbuf(1:ys_hier_column_a2a_send_elems), ys_hier_column_global_send_counts, &
                         ys_hier_column_global_send_displs, MPI_DOUBLE_COMPLEX, &
                         ycomm_recvbuf(1:ys_hier_column_a2a_recv_elems), ys_hier_column_global_recv_counts, &
                         ys_hier_column_global_recv_displs, MPI_DOUBLE_COMPLEX, ys_hier_column_comm, ierr_local)
      !$omp end target data
      if (ys_hier_comm_stats_enabled) then
        comm_elapsed = MPI_Wtime() - comm_t0
        call ys_hier_report_comm_stats("node_leaf_rows_to_global", ys_hier_column_comm, &
                                       ys_hier_column_global_send_counts, ys_hier_column_global_recv_counts, comm_elapsed)
      end if
      call roctxPop("MPI_Alltoallv node_leaf_rows_to_global")
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv node_leaf_rows_to_global failed"

      call ys_node_leaf_solve_global_groups_column_alltoall()
      call roctxPush("MPI_Alltoallv node_leaf_values_from_global")
      if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
      !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
      call MPI_Alltoallv(ycomm_sendbuf(1:ys_hier_column_return_send_elems), ys_hier_column_return_send_counts, &
                         ys_hier_column_return_send_displs, MPI_DOUBLE_COMPLEX, &
                         ycomm_recvbuf(1:ys_hier_column_return_recv_elems), ys_hier_column_return_recv_counts, &
                         ys_hier_column_return_recv_displs, MPI_DOUBLE_COMPLEX, ys_hier_column_comm, ierr_local)
      !$omp end target data
      if (ys_hier_comm_stats_enabled) then
        comm_elapsed = MPI_Wtime() - comm_t0
        call ys_hier_report_comm_stats("node_leaf_values_from_global", ys_hier_column_comm, &
                                       ys_hier_column_return_send_counts, ys_hier_column_return_recv_counts, comm_elapsed)
      end if
      call roctxPop("MPI_Alltoallv node_leaf_values_from_global")
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv node_leaf_values_from_global failed"
      call ys_node_leaf_unpack_column_returns(nlines)
    end select
  end subroutine ys_solve_reduced_interfaces_node_leaf

  subroutine ys_solve_reduced_interfaces_hier()
    implicit none
    integer :: ierr_local
    integer(C_INT) :: nlines, comm_send_elems, comm_recv_elems
    real(C_DOUBLE) :: comm_t0 = 0.0_C_DOUBLE
    real(C_DOUBLE) :: comm_elapsed = 0.0_C_DOUBLE

    nlines = size(ys_reduced_rows_send, 2)
    comm_send_elems = max(ys_hier_leaf_send_elems, ys_hier_global_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_column_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_column_a2a_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_column_return_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_group_return_send_elems)
    comm_send_elems = max(comm_send_elems, ys_hier_leaf_return_send_elems)
    comm_recv_elems = max(ys_hier_leaf_recv_elems, ys_hier_global_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_column_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_column_a2a_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_column_return_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_group_return_recv_elems)
    comm_recv_elems = max(comm_recv_elems, ys_hier_leaf_return_recv_elems)
    call ensure_ycomm_buffers(comm_send_elems, comm_recv_elems)

    call ys_hier_pack_leaf_rows(nlines)
    call roctxPush("MPI_Alltoallv hier_leaf_rows_to_group")
    if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
    !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
    call MPI_Alltoallv(ycomm_sendbuf(1:ys_hier_leaf_send_elems), &
                       ys_hier_leaf_send_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       ys_hier_leaf_send_displs(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       MPI_DOUBLE_COMPLEX, ycomm_recvbuf(1:ys_hier_leaf_recv_elems), &
                       ys_hier_leaf_recv_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       ys_hier_leaf_recv_displs(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       MPI_DOUBLE_COMPLEX, ys_hier_group_comm, ierr_local)
    !$omp end target data
    if (ys_hier_comm_stats_enabled) then
      comm_elapsed = MPI_Wtime() - comm_t0
      call ys_hier_report_comm_stats("hier_leaf_rows_to_group", ys_hier_group_comm, &
                                     ys_hier_leaf_send_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                                     ys_hier_leaf_recv_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                                     comm_elapsed)
    end if
    call roctxPop("MPI_Alltoallv hier_leaf_rows_to_group")
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv hier_leaf_rows_to_group failed"

    call ys_hier_compose_node_rows()
    if (ys_hier_use_column_global) then
      select case (ys_hier_global_exchange)
      case (YS_HIER_EXCHANGE_ALLGATHER)
        call ys_hier_pack_node_rows_column()
        call roctxPush("MPI_Allgather hier_node_rows_to_global")
        if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
        !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
        call MPI_Allgather(ycomm_sendbuf(1:ys_hier_column_send_elems), ys_hier_column_send_elems, MPI_DOUBLE_COMPLEX, &
                           ycomm_recvbuf(1:ys_hier_column_recv_elems), ys_hier_column_send_elems, MPI_DOUBLE_COMPLEX, &
                           ys_hier_column_comm, ierr_local)
        !$omp end target data
        if (ys_hier_comm_stats_enabled) then
          comm_elapsed = MPI_Wtime() - comm_t0
          call ys_hier_report_comm_stats_elems("hier_node_rows_to_global", ys_hier_column_comm, &
                                               ys_hier_column_send_elems, ys_hier_column_recv_elems, comm_elapsed)
        end if
        call roctxPop("MPI_Allgather hier_node_rows_to_global")
        if (ierr_local /= MPI_SUCCESS) error stop "MPI_Allgather hier_node_rows_to_global failed"

        call ys_hier_solve_global_groups_column()
        call ys_hier_recover_leaf_returns_column()
      case (YS_HIER_EXCHANGE_ALLTOALL)
        call ys_hier_pack_node_rows_to_column_global()
        call roctxPush("MPI_Alltoallv hier_node_rows_to_global")
        if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
        !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
        call MPI_Alltoallv(ycomm_sendbuf(1:ys_hier_column_a2a_send_elems), ys_hier_column_global_send_counts, &
                           ys_hier_column_global_send_displs, MPI_DOUBLE_COMPLEX, &
                           ycomm_recvbuf(1:ys_hier_column_a2a_recv_elems), ys_hier_column_global_recv_counts, &
                           ys_hier_column_global_recv_displs, MPI_DOUBLE_COMPLEX, ys_hier_column_comm, ierr_local)
        !$omp end target data
        if (ys_hier_comm_stats_enabled) then
          comm_elapsed = MPI_Wtime() - comm_t0
          call ys_hier_report_comm_stats("hier_node_rows_to_global", ys_hier_column_comm, &
                                         ys_hier_column_global_send_counts, ys_hier_column_global_recv_counts, &
                                         comm_elapsed)
        end if
        call roctxPop("MPI_Alltoallv hier_node_rows_to_global")
        if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv hier_node_rows_to_global failed"

        call ys_hier_solve_global_groups_column_alltoall()
        call roctxPush("MPI_Alltoallv hier_group_values_from_global")
        if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
        !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
        call MPI_Alltoallv(ycomm_sendbuf(1:ys_hier_column_return_send_elems), ys_hier_column_return_send_counts, &
                           ys_hier_column_return_send_displs, MPI_DOUBLE_COMPLEX, &
                           ycomm_recvbuf(1:ys_hier_column_return_recv_elems), ys_hier_column_return_recv_counts, &
                           ys_hier_column_return_recv_displs, MPI_DOUBLE_COMPLEX, ys_hier_column_comm, ierr_local)
        !$omp end target data
        if (ys_hier_comm_stats_enabled) then
          comm_elapsed = MPI_Wtime() - comm_t0
          call ys_hier_report_comm_stats("hier_group_values_from_global", ys_hier_column_comm, &
                                         ys_hier_column_return_send_counts, ys_hier_column_return_recv_counts, &
                                         comm_elapsed)
        end if
        call roctxPop("MPI_Alltoallv hier_group_values_from_global")
        if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv hier_group_values_from_global failed"

        call ys_hier_recover_leaf_returns_column_alltoall()
      end select
    else
      call ys_hier_pack_node_rows_to_global(nlines)
      call roctxPush("MPI_Alltoallv hier_node_rows_to_global")
      if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
      !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
      call MPI_Alltoallv(ycomm_sendbuf(1:ys_hier_global_send_elems), ys_hier_global_send_counts, ys_hier_global_send_displs, &
                         MPI_DOUBLE_COMPLEX, ycomm_recvbuf(1:ys_hier_global_recv_elems), ys_hier_global_recv_counts, &
                         ys_hier_global_recv_displs, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr_local)
      !$omp end target data
      if (ys_hier_comm_stats_enabled) then
        comm_elapsed = MPI_Wtime() - comm_t0
        call ys_hier_report_comm_stats("hier_node_rows_to_global", MPI_COMM_Y, ys_hier_global_send_counts, &
                                       ys_hier_global_recv_counts, comm_elapsed)
      end if
      call roctxPop("MPI_Alltoallv hier_node_rows_to_global")
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv hier_node_rows_to_global failed"

      call ys_hier_solve_global_groups(nlines)
      call roctxPush("MPI_Alltoallv hier_group_values_from_global")
      if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
      !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
      call MPI_Alltoallv(ycomm_sendbuf(1:ys_hier_group_return_send_elems), ys_hier_group_return_send_counts, &
                         ys_hier_group_return_send_displs, MPI_DOUBLE_COMPLEX, &
                         ycomm_recvbuf(1:ys_hier_group_return_recv_elems), ys_hier_group_return_recv_counts, &
                         ys_hier_group_return_recv_displs, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr_local)
      !$omp end target data
      if (ys_hier_comm_stats_enabled) then
        comm_elapsed = MPI_Wtime() - comm_t0
        call ys_hier_report_comm_stats("hier_group_values_from_global", MPI_COMM_Y, ys_hier_group_return_send_counts, &
                                       ys_hier_group_return_recv_counts, comm_elapsed)
      end if
      call roctxPop("MPI_Alltoallv hier_group_values_from_global")
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv hier_group_values_from_global failed"

      call ys_hier_recover_leaf_returns()
    end if
    call roctxPush("MPI_Alltoallv hier_leaf_values_from_group")
    if (ys_hier_comm_stats_enabled) comm_t0 = MPI_Wtime()
    !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
    call MPI_Alltoallv(ycomm_sendbuf(1:ys_hier_leaf_return_send_elems), &
                       ys_hier_leaf_return_send_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       ys_hier_leaf_return_send_displs(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       MPI_DOUBLE_COMPLEX, ycomm_recvbuf(1:ys_hier_leaf_return_recv_elems), &
                       ys_hier_leaf_return_recv_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       ys_hier_leaf_return_recv_displs(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       MPI_DOUBLE_COMPLEX, ys_hier_group_comm, ierr_local)
    !$omp end target data
    if (ys_hier_comm_stats_enabled) then
      comm_elapsed = MPI_Wtime() - comm_t0
      call ys_hier_report_comm_stats("hier_leaf_values_from_group", ys_hier_group_comm, &
                                     ys_hier_leaf_return_send_counts(ys_hier_group_first + 1:ys_hier_group_first + &
                                                                     ys_hier_group_count), &
                                     ys_hier_leaf_return_recv_counts(ys_hier_group_first + 1:ys_hier_group_first + &
                                                                     ys_hier_group_count), &
                                     comm_elapsed)
    end if
    call roctxPop("MPI_Alltoallv hier_leaf_values_from_group")
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv hier_leaf_values_from_group failed"

    call ys_hier_unpack_leaf_returns(nlines)
  end subroutine ys_solve_reduced_interfaces_hier

  subroutine ys_hier_report_comm_stats(label, comm, send_counts, recv_counts, elapsed)
    implicit none
    character(*), intent(in) :: label
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: send_counts(:), recv_counts(:)
    real(C_DOUBLE), intent(in) :: elapsed

    call ys_hier_report_comm_stats_elems(label, comm, sum(send_counts), sum(recv_counts), elapsed)
  end subroutine ys_hier_report_comm_stats

  subroutine ys_hier_report_comm_stats_elems(label, comm, send_elems, recv_elems, elapsed)
    implicit none
    character(*), intent(in) :: label
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: send_elems, recv_elems
    real(C_DOUBLE), intent(in) :: elapsed
    integer :: comm_rank, comm_size, ierr_local
    real(C_DOUBLE) :: local_send_bytes, local_recv_bytes
    real(C_DOUBLE) :: send_bandwidth_gbs, recv_bandwidth_gbs

    local_send_bytes = 16.0_C_DOUBLE*real(send_elems, C_DOUBLE)
    local_recv_bytes = 16.0_C_DOUBLE*real(recv_elems, C_DOUBLE)
    call MPI_Comm_rank(comm, comm_rank, ierr_local)
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_rank y-Schur comm stats failed"
    call MPI_Comm_size(comm, comm_size, ierr_local)
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_size y-Schur comm stats failed"

    if (elapsed > 0.0_C_DOUBLE) then
      send_bandwidth_gbs = local_send_bytes/elapsed/1.0e9_C_DOUBLE
      recv_bandwidth_gbs = local_recv_bytes/elapsed/1.0e9_C_DOUBLE
    else
      send_bandwidth_gbs = 0.0_C_DOUBLE
      recv_bandwidth_gbs = 0.0_C_DOUBLE
    end if

    write (*, '("Y_SCHUR_COMM_STATS label=",a," ipy=",i0," comm_rank=",i0," comm_size=",i0,'// &
           '" send_MB=",f12.6," recv_MB=",f12.6," time_ms=",f12.6,'// &
           '" send_GBps=",f12.6," recv_GBps=",f12.6)') &
      trim(label), ipy, comm_rank, comm_size, local_send_bytes/1.0e6_C_DOUBLE, local_recv_bytes/1.0e6_C_DOUBLE, &
      elapsed*1.0e3_C_DOUBLE, send_bandwidth_gbs, recv_bandwidth_gbs
  end subroutine ys_hier_report_comm_stats_elems

  subroutine ys_node_leaf_pack_rows_column(nlines)
    implicit none
    integer(C_INT), intent(in) :: nlines
    integer(C_INT) :: iline, irow, offset

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_reduced_rows_send, ycomm_sendbuf, nlines) &
    !$omp private(iline, irow, offset)
    do iline = 1, nlines
      do irow = 1, 20
        offset = (iline - 1)*20 + irow
        ycomm_sendbuf(offset) = ys_reduced_rows_send(irow, iline)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_node_leaf_pack_rows_column

  subroutine ys_node_leaf_pack_rows_to_column_global(nlines)
    implicit none
    integer(C_INT), intent(in) :: nlines
    integer(C_INT) :: iline, irow, owner_rank, local_row0, offset
    integer(C_INT) :: base_count, remainder, split_line

    base_count = nlines/ys_hier_column_count
    remainder = mod(nlines, ys_hier_column_count)
    split_line = (base_count + 1_C_INT)*remainder
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_reduced_rows_send, ycomm_sendbuf, ys_hier_column_global_send_displs, nlines, &
    !$omp& base_count, remainder, split_line) &
    !$omp private(iline, irow, owner_rank, local_row0, offset)
    do iline = 1, nlines
      do irow = 1, 20
        if (iline <= split_line) then
          owner_rank = (iline - 1)/(base_count + 1_C_INT)
          local_row0 = iline - owner_rank*(base_count + 1_C_INT)
        else
          owner_rank = remainder + (iline - split_line - 1_C_INT)/base_count
          local_row0 = iline - split_line - (owner_rank - remainder)*base_count
        end if
        offset = ys_hier_column_global_send_displs(owner_rank + 1) + (local_row0 - 1)*20 + irow
        ycomm_sendbuf(offset) = ys_reduced_rows_send(irow, iline)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_node_leaf_pack_rows_to_column_global

  subroutine ys_node_leaf_solve_global_groups_column(nlines)
    implicit none
    integer(C_INT), intent(in) :: nlines
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: iline, iblock, row0, offset

    call roctxPush("ys_node_leaf_global_group_solve")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_recvbuf, ys_reduced_matrix_lu, ys_reduced_rhs, ys_left_interface_values, &
    !$omp& ys_right_interface_values, nlines, ys_hier_ngroups, ys_hier_group_id, bw) &
    !$omp private(iline, iblock, row0, offset)
    do iline = 1, nlines
      ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline) = (0.0d0, 0.0d0)
      ys_reduced_rhs(1:4*ys_hier_ngroups, iline) = (0.0d0, 0.0d0)
      do iblock = 0, ys_hier_ngroups - 1
        offset = iblock*(20*nlines) + (iline - 1)*20
        row0 = 4*iblock
        ys_reduced_matrix_lu(row0 + 1, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 2, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 3, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 4, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_rhs(row0 + 1, iline) = ycomm_recvbuf(offset + 1)
        ys_reduced_rhs(row0 + 2, iline) = ycomm_recvbuf(offset + 2)
        ys_reduced_rhs(row0 + 3, iline) = ycomm_recvbuf(offset + 3)
        ys_reduced_rhs(row0 + 4, iline) = ycomm_recvbuf(offset + 4)
        if (iblock > 0) then
          ys_reduced_matrix_lu(row0 + 1, bw - 1, iline) = -ycomm_recvbuf(offset + 5)
          ys_reduced_matrix_lu(row0 + 2, bw - 2, iline) = -ycomm_recvbuf(offset + 6)
          ys_reduced_matrix_lu(row0 + 3, bw - 3, iline) = -ycomm_recvbuf(offset + 7)
          ys_reduced_matrix_lu(row0 + 4, bw - 4, iline) = -ycomm_recvbuf(offset + 8)
          ys_reduced_matrix_lu(row0 + 1, bw, iline) = -ycomm_recvbuf(offset + 9)
          ys_reduced_matrix_lu(row0 + 2, bw - 1, iline) = -ycomm_recvbuf(offset + 10)
          ys_reduced_matrix_lu(row0 + 3, bw - 2, iline) = -ycomm_recvbuf(offset + 11)
          ys_reduced_matrix_lu(row0 + 4, bw - 3, iline) = -ycomm_recvbuf(offset + 12)
        end if
        if (iblock < ys_hier_ngroups - 1) then
          ys_reduced_matrix_lu(row0 + 1, bw + 5, iline) = -ycomm_recvbuf(offset + 13)
          ys_reduced_matrix_lu(row0 + 2, bw + 4, iline) = -ycomm_recvbuf(offset + 14)
          ys_reduced_matrix_lu(row0 + 3, bw + 3, iline) = -ycomm_recvbuf(offset + 15)
          ys_reduced_matrix_lu(row0 + 4, bw + 2, iline) = -ycomm_recvbuf(offset + 16)
          ys_reduced_matrix_lu(row0 + 1, bw + 6, iline) = -ycomm_recvbuf(offset + 17)
          ys_reduced_matrix_lu(row0 + 2, bw + 5, iline) = -ycomm_recvbuf(offset + 18)
          ys_reduced_matrix_lu(row0 + 3, bw + 4, iline) = -ycomm_recvbuf(offset + 19)
          ys_reduced_matrix_lu(row0 + 4, bw + 3, iline) = -ycomm_recvbuf(offset + 20)
        end if
      end do

      call ys_factor_banded_complex(ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))
      call ys_solve_factored_banded_complex(ys_reduced_rhs(1:4*ys_hier_ngroups, iline), &
                                            ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))

      row0 = 4*ys_hier_group_id
      ys_left_interface_values(:, iline) = (0.0d0, 0.0d0)
      ys_right_interface_values(:, iline) = (0.0d0, 0.0d0)
      if (ys_hier_group_id > 0) then
        ys_left_interface_values(1, iline) = ys_reduced_rhs(row0 - 1, iline)
        ys_left_interface_values(2, iline) = ys_reduced_rhs(row0, iline)
      end if
      if (ys_hier_group_id < ys_hier_ngroups - 1) then
        ys_right_interface_values(1, iline) = ys_reduced_rhs(row0 + 5, iline)
        ys_right_interface_values(2, iline) = ys_reduced_rhs(row0 + 6, iline)
      end if
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_node_leaf_global_group_solve")
  end subroutine ys_node_leaf_solve_global_groups_column

  subroutine ys_node_leaf_solve_global_groups_column_alltoall()
    implicit none
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: iline, iblock, row0, offset

    call roctxPush("ys_node_leaf_global_group_solve")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_recvbuf, ycomm_sendbuf, ys_hier_column_global_recv_displs, &
    !$omp& ys_hier_column_return_send_displs, ys_hier_column_solve_nlines, ys_hier_ngroups, &
    !$omp& ys_reduced_matrix_lu, ys_reduced_rhs, bw) &
    !$omp private(iline, iblock, row0, offset)
    do iline = 1, ys_hier_column_solve_nlines
      ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline) = (0.0d0, 0.0d0)
      ys_reduced_rhs(1:4*ys_hier_ngroups, iline) = (0.0d0, 0.0d0)
      do iblock = 0, ys_hier_ngroups - 1
        offset = ys_hier_column_global_recv_displs(iblock + 1) + (iline - 1)*20
        row0 = 4*iblock
        ys_reduced_matrix_lu(row0 + 1, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 2, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 3, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 4, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_rhs(row0 + 1, iline) = ycomm_recvbuf(offset + 1)
        ys_reduced_rhs(row0 + 2, iline) = ycomm_recvbuf(offset + 2)
        ys_reduced_rhs(row0 + 3, iline) = ycomm_recvbuf(offset + 3)
        ys_reduced_rhs(row0 + 4, iline) = ycomm_recvbuf(offset + 4)
        if (iblock > 0) then
          ys_reduced_matrix_lu(row0 + 1, bw - 1, iline) = -ycomm_recvbuf(offset + 5)
          ys_reduced_matrix_lu(row0 + 2, bw - 2, iline) = -ycomm_recvbuf(offset + 6)
          ys_reduced_matrix_lu(row0 + 3, bw - 3, iline) = -ycomm_recvbuf(offset + 7)
          ys_reduced_matrix_lu(row0 + 4, bw - 4, iline) = -ycomm_recvbuf(offset + 8)
          ys_reduced_matrix_lu(row0 + 1, bw, iline) = -ycomm_recvbuf(offset + 9)
          ys_reduced_matrix_lu(row0 + 2, bw - 1, iline) = -ycomm_recvbuf(offset + 10)
          ys_reduced_matrix_lu(row0 + 3, bw - 2, iline) = -ycomm_recvbuf(offset + 11)
          ys_reduced_matrix_lu(row0 + 4, bw - 3, iline) = -ycomm_recvbuf(offset + 12)
        end if
        if (iblock < ys_hier_ngroups - 1) then
          ys_reduced_matrix_lu(row0 + 1, bw + 5, iline) = -ycomm_recvbuf(offset + 13)
          ys_reduced_matrix_lu(row0 + 2, bw + 4, iline) = -ycomm_recvbuf(offset + 14)
          ys_reduced_matrix_lu(row0 + 3, bw + 3, iline) = -ycomm_recvbuf(offset + 15)
          ys_reduced_matrix_lu(row0 + 4, bw + 2, iline) = -ycomm_recvbuf(offset + 16)
          ys_reduced_matrix_lu(row0 + 1, bw + 6, iline) = -ycomm_recvbuf(offset + 17)
          ys_reduced_matrix_lu(row0 + 2, bw + 5, iline) = -ycomm_recvbuf(offset + 18)
          ys_reduced_matrix_lu(row0 + 3, bw + 4, iline) = -ycomm_recvbuf(offset + 19)
          ys_reduced_matrix_lu(row0 + 4, bw + 3, iline) = -ycomm_recvbuf(offset + 20)
        end if
      end do

      call ys_factor_banded_complex(ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))
      call ys_solve_factored_banded_complex(ys_reduced_rhs(1:4*ys_hier_ngroups, iline), &
                                            ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))
      do iblock = 0, ys_hier_ngroups - 1
        offset = ys_hier_column_return_send_displs(iblock + 1) + (iline - 1)*YS_REDUCED_RETURN_WIDTH
        row0 = 4*iblock
        ycomm_sendbuf(offset + 1) = ys_reduced_rhs(row0 + 1, iline)
        ycomm_sendbuf(offset + 2) = ys_reduced_rhs(row0 + 2, iline)
        ycomm_sendbuf(offset + 3) = ys_reduced_rhs(row0 + 3, iline)
        ycomm_sendbuf(offset + 4) = ys_reduced_rhs(row0 + 4, iline)
        ycomm_sendbuf(offset + 5) = (0.0d0, 0.0d0)
        ycomm_sendbuf(offset + 6) = (0.0d0, 0.0d0)
        ycomm_sendbuf(offset + 7) = (0.0d0, 0.0d0)
        ycomm_sendbuf(offset + 8) = (0.0d0, 0.0d0)
        if (iblock > 0) then
          ycomm_sendbuf(offset + 5) = ys_reduced_rhs(row0 - 1, iline)
          ycomm_sendbuf(offset + 6) = ys_reduced_rhs(row0, iline)
        end if
        if (iblock < ys_hier_ngroups - 1) then
          ycomm_sendbuf(offset + 7) = ys_reduced_rhs(row0 + 5, iline)
          ycomm_sendbuf(offset + 8) = ys_reduced_rhs(row0 + 6, iline)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_node_leaf_global_group_solve")
  end subroutine ys_node_leaf_solve_global_groups_column_alltoall

  subroutine ys_node_leaf_unpack_column_returns(nlines)
    implicit none
    integer(C_INT), intent(in) :: nlines
    integer(C_INT) :: iline, k, src_rank, src_first, offset, row0
    integer(C_INT) :: base_count, remainder, split_line

    row0 = 4*ys_hier_group_id
    base_count = nlines/ys_hier_column_count
    remainder = mod(nlines, ys_hier_column_count)
    split_line = (base_count + 1_C_INT)*remainder
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_recvbuf, ys_reduced_rhs, ys_left_interface_values, ys_right_interface_values, &
    !$omp& ys_hier_column_return_recv_displs, nlines, base_count, remainder, split_line, row0) &
    !$omp private(iline, k, src_rank, src_first, offset)
    do iline = 1, nlines
      if (iline <= split_line) then
        src_rank = (iline - 1)/(base_count + 1_C_INT)
        src_first = src_rank*(base_count + 1_C_INT) + 1
      else
        src_rank = remainder + (iline - split_line - 1_C_INT)/base_count
        src_first = split_line + (src_rank - remainder)*base_count + 1
      end if
      offset = ys_hier_column_return_recv_displs(src_rank + 1) + &
               (iline - src_first)*YS_REDUCED_RETURN_WIDTH
      do k = 1, 4
        ys_reduced_rhs(row0 + k, iline) = ycomm_recvbuf(offset + k)
      end do
      ys_left_interface_values(1, iline) = ycomm_recvbuf(offset + 5)
      ys_left_interface_values(2, iline) = ycomm_recvbuf(offset + 6)
      ys_right_interface_values(1, iline) = ycomm_recvbuf(offset + 7)
      ys_right_interface_values(2, iline) = ycomm_recvbuf(offset + 8)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_node_leaf_unpack_column_returns

  subroutine ys_node_leaf_transpose_source_to_group(field, slab, ny, nz, total_lines, nlines_z)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, total_lines, nlines_z
    complex(C_DOUBLE_COMPLEX), intent(in) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: slab(:, :)
    integer :: ierr_local
    integer(C_INT) :: rank, local_rank, y_first, y_last, rows, first_line, line_count
    integer(C_INT) :: send_offset, recv_offset, send_elems, recv_elems
    integer(C_INT) :: node_y0, node_yN, send_y_first, send_y_last, send_rows
    integer(C_INT) :: nx_lines, base_count, remainder, split_line, max_group_lines
    integer(C_INT) :: dest_local, src_local, ilocal, iline, ix, iz, iy, p, slab_y

    node_y0 = 1_C_INT + ys_hier_group_first*(ny - 1_C_INT)/npy_grid
    node_yN = (ys_hier_group_first + ys_hier_group_count)*(ny - 1_C_INT)/npy_grid
    send_y_first = ny0
    send_y_last = nyN
    if (ys_hier_group_rank == 0) send_y_first = ny0 - 2_C_INT
    if (ys_hier_group_rank == ys_hier_group_count - 1) send_y_last = nyN + 2_C_INT
    send_rows = send_y_last - send_y_first + 1_C_INT

    send_offset = 0_C_INT
    recv_offset = 0_C_INT
    do local_rank = 0, ys_hier_group_count - 1
      rank = ys_hier_group_first + local_rank
      nx_lines = total_lines/nlines_z
      call ys_split_line_range(local_rank, nx_lines, ys_hier_group_count, first_line, line_count)
      first_line = (first_line - 1_C_INT)*nlines_z + 1_C_INT
      line_count = line_count*nlines_z
      ys_hier_leaf_send_counts(rank + 1) = line_count*send_rows
      ys_hier_leaf_send_displs(rank + 1) = send_offset
      send_offset = send_offset + ys_hier_leaf_send_counts(rank + 1)

      y_first = 1_C_INT + rank*(ny - 1_C_INT)/npy_grid
      y_last = (rank + 1_C_INT)*(ny - 1_C_INT)/npy_grid
      if (local_rank == 0) y_first = y_first - 2_C_INT
      if (local_rank == ys_hier_group_count - 1) y_last = y_last + 2_C_INT
      rows = y_last - y_first + 1_C_INT
      ys_hier_leaf_recv_counts(rank + 1) = ys_workspace_nlines*rows
      ys_hier_leaf_recv_displs(rank + 1) = recv_offset
      recv_offset = recv_offset + ys_hier_leaf_recv_counts(rank + 1)
    end do
    send_elems = send_offset
    recv_elems = recv_offset
    call ensure_ycomm_buffers(send_elems, recv_elems)
    !$omp target update to(ys_hier_leaf_send_counts, ys_hier_leaf_send_displs, &
    !$omp& ys_hier_leaf_recv_counts, ys_hier_leaf_recv_displs)

    nx_lines = total_lines/nlines_z
    base_count = nx_lines/ys_hier_group_count
    remainder = mod(nx_lines, ys_hier_group_count)
    split_line = (base_count + 1_C_INT)*remainder*nlines_z
    max_group_lines = (base_count + merge(1_C_INT, 0_C_INT, remainder > 0))*nlines_z
    call roctxPush("node_leaf_source_to_group pack")
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(field, ycomm_sendbuf, ys_hier_leaf_send_displs, total_lines, nlines_z, nz, nx0, &
    !$omp& ys_hier_group_first, ys_hier_group_count, base_count, remainder, split_line, &
    !$omp& max_group_lines, send_y_first, send_y_last, send_rows) &
    !$omp private(dest_local, ilocal, iy, first_line, line_count, iline, ix, iz, p)
    do dest_local = 0, ys_hier_group_count - 1
      do ilocal = 1, max_group_lines
        do iy = send_y_first, send_y_last
          if (dest_local < remainder) then
            line_count = (base_count + 1_C_INT)*nlines_z
            first_line = dest_local*(base_count + 1_C_INT)*nlines_z + 1_C_INT
          else
            line_count = base_count*nlines_z
            first_line = split_line + (dest_local - remainder)*base_count*nlines_z + 1_C_INT
          end if
          if (ilocal > line_count) cycle
          iline = first_line + ilocal - 1_C_INT
          ix = (iline - 1_C_INT)/nlines_z + nx0
          iz = mod(iline - 1_C_INT, nlines_z) - nz
          p = ys_hier_leaf_send_displs(ys_hier_group_first + dest_local + 1) + &
              (ilocal - 1_C_INT)*send_rows + (iy - send_y_first + 1_C_INT)
          ycomm_sendbuf(p) = field(iy, iz, ix)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("node_leaf_source_to_group pack")

    call roctxPush("MPI_Alltoallv node_leaf_source_to_group")
    !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
    call MPI_Alltoallv(ycomm_sendbuf(1:send_elems), &
                       ys_hier_leaf_send_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       ys_hier_leaf_send_displs(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       MPI_DOUBLE_COMPLEX, ycomm_recvbuf(1:recv_elems), &
                       ys_hier_leaf_recv_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       ys_hier_leaf_recv_displs(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       MPI_DOUBLE_COMPLEX, ys_hier_group_comm, ierr_local)
    !$omp end target data
    call roctxPop("MPI_Alltoallv node_leaf_source_to_group")
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv node_leaf_source_to_group failed"

    call roctxPush("node_leaf_source_to_group unpack")
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(slab, ycomm_recvbuf, ys_hier_leaf_recv_displs, ny, npy_grid, ys_hier_group_first, &
    !$omp& ys_hier_group_count, ys_workspace_nlines, node_y0, node_yN) &
    !$omp private(src_local, ilocal, iy, rank, y_first, y_last, rows, p, slab_y)
    do src_local = 0, ys_hier_group_count - 1
      do ilocal = 1, ys_workspace_nlines
        do iy = node_y0 - 2_C_INT, node_yN + 2_C_INT
          rank = ys_hier_group_first + src_local
          y_first = 1_C_INT + rank*(ny - 1_C_INT)/npy_grid
          y_last = (rank + 1_C_INT)*(ny - 1_C_INT)/npy_grid
          if (src_local == 0) y_first = y_first - 2_C_INT
          if (src_local == ys_hier_group_count - 1) y_last = y_last + 2_C_INT
          if (iy < y_first .or. iy > y_last) cycle
          rows = y_last - y_first + 1_C_INT
          p = ys_hier_leaf_recv_displs(rank + 1) + (ilocal - 1_C_INT)*rows + (iy - y_first + 1_C_INT)
          slab_y = iy - node_y0 + 3_C_INT
          slab(slab_y, ilocal) = ycomm_recvbuf(p)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("node_leaf_source_to_group unpack")
  end subroutine ys_node_leaf_transpose_source_to_group

  subroutine ys_node_leaf_transpose_solution_from_group(slab, field, ny, nz, total_lines, nlines_z)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, total_lines, nlines_z
    complex(C_DOUBLE_COMPLEX), intent(in) :: slab(:, :)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer :: ierr_local
    integer(C_INT) :: rank, local_rank, y_first, y_last, rows, first_line, line_count
    integer(C_INT) :: local_y_first, local_y_last, local_rows, send_offset, recv_offset, send_elems, recv_elems
    integer(C_INT) :: node_y0, node_yN, nx_lines, base_count, remainder, split_line, max_group_lines
    integer(C_INT) :: dest_local, src_local, ilocal, iline, ix, iz, iy, p, slab_y

    node_y0 = 1_C_INT + ys_hier_group_first*(ny - 1_C_INT)/npy_grid
    node_yN = (ys_hier_group_first + ys_hier_group_count)*(ny - 1_C_INT)/npy_grid
    local_y_first = ny0 - 2_C_INT
    local_y_last = nyN + 2_C_INT
    local_rows = local_y_last - local_y_first + 1_C_INT

    send_offset = 0_C_INT
    recv_offset = 0_C_INT
    do local_rank = 0, ys_hier_group_count - 1
      rank = ys_hier_group_first + local_rank
      y_first = 1_C_INT + rank*(ny - 1_C_INT)/npy_grid - 2_C_INT
      y_last = (rank + 1_C_INT)*(ny - 1_C_INT)/npy_grid + 2_C_INT
      rows = y_last - y_first + 1_C_INT
      ys_hier_leaf_send_counts(rank + 1) = ys_workspace_nlines*rows
      ys_hier_leaf_send_displs(rank + 1) = send_offset
      send_offset = send_offset + ys_hier_leaf_send_counts(rank + 1)

      nx_lines = total_lines/nlines_z
      call ys_split_line_range(local_rank, nx_lines, ys_hier_group_count, first_line, line_count)
      first_line = (first_line - 1_C_INT)*nlines_z + 1_C_INT
      line_count = line_count*nlines_z
      ys_hier_leaf_recv_counts(rank + 1) = line_count*local_rows
      ys_hier_leaf_recv_displs(rank + 1) = recv_offset
      recv_offset = recv_offset + ys_hier_leaf_recv_counts(rank + 1)
    end do
    send_elems = send_offset
    recv_elems = recv_offset
    call ensure_ycomm_buffers(send_elems, recv_elems)
    !$omp target update to(ys_hier_leaf_send_counts, ys_hier_leaf_send_displs, &
    !$omp& ys_hier_leaf_recv_counts, ys_hier_leaf_recv_displs)

    call roctxPush("node_leaf_solution_from_group pack")
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(slab, ycomm_sendbuf, ys_hier_leaf_send_displs, ny, npy_grid, ys_hier_group_first, &
    !$omp& ys_hier_group_count, ys_workspace_nlines, node_y0, node_yN) &
    !$omp private(dest_local, ilocal, iy, rank, y_first, y_last, rows, p, slab_y)
    do dest_local = 0, ys_hier_group_count - 1
      do ilocal = 1, ys_workspace_nlines
        do iy = node_y0 - 2_C_INT, node_yN + 2_C_INT
          rank = ys_hier_group_first + dest_local
          y_first = 1_C_INT + rank*(ny - 1_C_INT)/npy_grid - 2_C_INT
          y_last = (rank + 1_C_INT)*(ny - 1_C_INT)/npy_grid + 2_C_INT
          if (iy < y_first .or. iy > y_last) cycle
          rows = y_last - y_first + 1_C_INT
          slab_y = iy - node_y0 + 3_C_INT
          p = ys_hier_leaf_send_displs(rank + 1) + (ilocal - 1_C_INT)*rows + (iy - y_first + 1_C_INT)
          ycomm_sendbuf(p) = slab(slab_y, ilocal)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("node_leaf_solution_from_group pack")

    call roctxPush("MPI_Alltoallv node_leaf_solution_from_group")
    !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
    call MPI_Alltoallv(ycomm_sendbuf(1:send_elems), &
                       ys_hier_leaf_send_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       ys_hier_leaf_send_displs(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       MPI_DOUBLE_COMPLEX, ycomm_recvbuf(1:recv_elems), &
                       ys_hier_leaf_recv_counts(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       ys_hier_leaf_recv_displs(ys_hier_group_first + 1:ys_hier_group_first + ys_hier_group_count), &
                       MPI_DOUBLE_COMPLEX, ys_hier_group_comm, ierr_local)
    !$omp end target data
    call roctxPop("MPI_Alltoallv node_leaf_solution_from_group")
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv node_leaf_solution_from_group failed"

    nx_lines = total_lines/nlines_z
    base_count = nx_lines/ys_hier_group_count
    remainder = mod(nx_lines, ys_hier_group_count)
    split_line = (base_count + 1_C_INT)*remainder*nlines_z
    max_group_lines = (base_count + merge(1_C_INT, 0_C_INT, remainder > 0))*nlines_z
    call roctxPush("node_leaf_solution_from_group unpack")
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(field, ycomm_recvbuf, ys_hier_leaf_recv_displs, total_lines, nlines_z, nz, nx0, &
    !$omp& ys_hier_group_first, ys_hier_group_count, base_count, remainder, split_line, &
    !$omp& max_group_lines, local_y_first, local_y_last, local_rows) &
    !$omp private(src_local, ilocal, iy, first_line, line_count, iline, ix, iz, p)
    do src_local = 0, ys_hier_group_count - 1
      do ilocal = 1, max_group_lines
        do iy = local_y_first, local_y_last
          if (src_local < remainder) then
            line_count = (base_count + 1_C_INT)*nlines_z
            first_line = src_local*(base_count + 1_C_INT)*nlines_z + 1_C_INT
          else
            line_count = base_count*nlines_z
            first_line = split_line + (src_local - remainder)*base_count*nlines_z + 1_C_INT
          end if
          if (ilocal > line_count) cycle
          iline = first_line + ilocal - 1_C_INT
          ix = (iline - 1_C_INT)/nlines_z + nx0
          iz = mod(iline - 1_C_INT, nlines_z) - nz
          p = ys_hier_leaf_recv_displs(ys_hier_group_first + src_local + 1) + &
              (ilocal - 1_C_INT)*local_rows + (iy - local_y_first + 1_C_INT)
          field(iy, iz, ix) = ycomm_recvbuf(p)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("node_leaf_solution_from_group unpack")
  end subroutine ys_node_leaf_transpose_solution_from_group

  subroutine ys_hier_pack_leaf_rows(nlines)
    implicit none
    integer(C_INT), intent(in) :: nlines
    integer(C_INT) :: iline, irow, owner_rank, owner_local, local_row0, offset
    integer(C_INT) :: base_count, remainder, split_line

    base_count = nlines/ys_hier_group_count
    remainder = mod(nlines, ys_hier_group_count)
    split_line = (base_count + 1_C_INT)*remainder
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_reduced_rows_send, ycomm_sendbuf, ys_hier_leaf_send_displs, nlines, ys_hier_group_first, &
    !$omp& base_count, remainder, split_line) &
    !$omp private(iline, irow, owner_rank, owner_local, local_row0, offset)
    do iline = 1, nlines
      do irow = 1, 20
        if (iline <= split_line) then
          owner_local = (iline - 1)/(base_count + 1_C_INT)
          local_row0 = iline - owner_local*(base_count + 1_C_INT)
        else
          owner_local = remainder + (iline - split_line - 1_C_INT)/base_count
          local_row0 = iline - split_line - (owner_local - remainder)*base_count
        end if
        owner_rank = ys_hier_group_first + owner_local
        offset = ys_hier_leaf_send_displs(owner_rank + 1) + (local_row0 - 1)*20 + irow
        ycomm_sendbuf(offset) = ys_reduced_rows_send(irow, iline)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_hier_pack_leaf_rows

  subroutine ys_hier_compose_node_rows()
    implicit none
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: iline, child, child_rank, row0, row, k, ext_col, exposed_var, offset, col

    call roctxPush("ys_hier_compose_node_rows")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_recvbuf, ys_hier_leaf_recv_displs, ys_reduced_matrix_lu, ys_hier_recover_basis, &
    !$omp& ys_hier_node_rows, ys_hier_node_owned_nlines, ys_hier_group_first, ys_hier_group_count, bw) &
    !$omp private(iline, child, child_rank, row0, row, k, ext_col, exposed_var, offset, col)
    do iline = 1, ys_hier_node_owned_nlines
      ys_reduced_matrix_lu(1:4*ys_hier_group_count, :, iline) = (0.0d0, 0.0d0)
      ys_hier_recover_basis(1:4*ys_hier_group_count, :, iline) = (0.0d0, 0.0d0)
      do child = 0, ys_hier_group_count - 1
        child_rank = ys_hier_group_first + child
        offset = ys_hier_leaf_recv_displs(child_rank + 1) + (iline - 1)*20
        row0 = 4*child
        do k = 1, 4
          row = row0 + k
          ys_reduced_matrix_lu(row, bw + 1, iline) = (1.0d0, 0.0d0)
          ys_hier_recover_basis(row, 1, iline) = ycomm_recvbuf(offset + k)
          if (child > 0) then
            col = row0 - 1
            ys_reduced_matrix_lu(row, bw + 1 + col - row, iline) = -ycomm_recvbuf(offset + 4 + k)
            col = row0
            ys_reduced_matrix_lu(row, bw + 1 + col - row, iline) = -ycomm_recvbuf(offset + 8 + k)
          else
            ys_hier_recover_basis(row, 2, iline) = ycomm_recvbuf(offset + 4 + k)
            ys_hier_recover_basis(row, 3, iline) = ycomm_recvbuf(offset + 8 + k)
          end if
          if (child < ys_hier_group_count - 1) then
            col = row0 + 5
            ys_reduced_matrix_lu(row, bw + 1 + col - row, iline) = -ycomm_recvbuf(offset + 12 + k)
            col = row0 + 6
            ys_reduced_matrix_lu(row, bw + 1 + col - row, iline) = -ycomm_recvbuf(offset + 16 + k)
          else
            ys_hier_recover_basis(row, 4, iline) = ycomm_recvbuf(offset + 12 + k)
            ys_hier_recover_basis(row, 5, iline) = ycomm_recvbuf(offset + 16 + k)
          end if
        end do
      end do

      call ys_factor_banded_complex(ys_reduced_matrix_lu(1:4*ys_hier_group_count, :, iline))
      call ys_solve_factored_banded_complex_multi(ys_hier_recover_basis(1:4*ys_hier_group_count, 1:5, iline), &
                                                  ys_reduced_matrix_lu(1:4*ys_hier_group_count, :, iline))

      do k = 1, 4
        if (k <= 2) then
          exposed_var = k
        else
          exposed_var = 4*(ys_hier_group_count - 1) + k
        end if
        ys_hier_node_rows(k, iline) = ys_hier_recover_basis(exposed_var, 1, iline)
        do ext_col = 1, 4
          ys_hier_node_rows(4*ext_col + k, iline) = ys_hier_recover_basis(exposed_var, ext_col + 1, iline)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_hier_compose_node_rows")
  end subroutine ys_hier_compose_node_rows

  subroutine ys_hier_pack_node_rows_column()
    implicit none
    integer(C_INT) :: iline, irow, offset

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_hier_node_rows, ycomm_sendbuf, ys_hier_node_owned_nlines) &
    !$omp private(iline, irow, offset)
    do iline = 1, ys_hier_node_owned_nlines
      do irow = 1, 20
        offset = (iline - 1)*20 + irow
        ycomm_sendbuf(offset) = ys_hier_node_rows(irow, iline)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_hier_pack_node_rows_column

  subroutine ys_hier_pack_node_rows_to_column_global()
    implicit none
    integer(C_INT) :: iline, irow, owner_rank, local_row0, offset
    integer(C_INT) :: base_count, remainder, split_line

    base_count = ys_hier_node_owned_nlines/ys_hier_column_count
    remainder = mod(ys_hier_node_owned_nlines, ys_hier_column_count)
    split_line = (base_count + 1_C_INT)*remainder
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_hier_node_rows, ycomm_sendbuf, ys_hier_column_global_send_displs, ys_hier_node_owned_nlines, &
    !$omp& base_count, remainder, split_line) &
    !$omp private(iline, irow, owner_rank, local_row0, offset)
    do iline = 1, ys_hier_node_owned_nlines
      do irow = 1, 20
        if (iline <= split_line) then
          owner_rank = (iline - 1)/(base_count + 1_C_INT)
          local_row0 = iline - owner_rank*(base_count + 1_C_INT)
        else
          owner_rank = remainder + (iline - split_line - 1_C_INT)/base_count
          local_row0 = iline - split_line - (owner_rank - remainder)*base_count
        end if
        offset = ys_hier_column_global_send_displs(owner_rank + 1) + (local_row0 - 1)*20 + irow
        ycomm_sendbuf(offset) = ys_hier_node_rows(irow, iline)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_hier_pack_node_rows_to_column_global

  subroutine ys_hier_pack_node_rows_to_global(nlines)
    implicit none
    integer(C_INT), intent(in) :: nlines
    integer(C_INT) :: iline, irow, global_line, owner_rank, local_row0, offset
    integer(C_INT) :: base_count, remainder, split_line, owner_first, intersection_first

    base_count = nlines/npy_grid
    remainder = mod(nlines, npy_grid)
    split_line = (base_count + 1_C_INT)*remainder
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_hier_node_rows, ycomm_sendbuf, ys_hier_global_send_displs, ys_hier_global_line_first, &
    !$omp& ys_hier_node_first_line, ys_hier_node_owned_nlines, base_count, remainder, split_line) &
    !$omp private(iline, irow, global_line, owner_rank, local_row0, offset, owner_first, intersection_first)
    do iline = 1, ys_hier_node_owned_nlines
      do irow = 1, 20
        global_line = ys_hier_node_first_line + iline - 1
        if (global_line <= split_line) then
          owner_rank = (global_line - 1)/(base_count + 1_C_INT)
          local_row0 = global_line - owner_rank*(base_count + 1_C_INT)
        else
          owner_rank = remainder + (global_line - split_line - 1_C_INT)/base_count
          local_row0 = global_line - split_line - (owner_rank - remainder)*base_count
        end if
        owner_first = ys_hier_global_line_first(owner_rank + 1)
        intersection_first = max(ys_hier_node_first_line, owner_first)
        offset = ys_hier_global_send_displs(owner_rank + 1) + (global_line - intersection_first)*20 + irow
        ycomm_sendbuf(offset) = ys_hier_node_rows(irow, iline)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_hier_pack_node_rows_to_global

  subroutine ys_hier_solve_global_groups_column()
    implicit none
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: iline, iblock, row0, offset

    call roctxPush("ys_hier_global_group_solve")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_recvbuf, ys_hier_group_values, ys_hier_node_owned_nlines, ys_hier_ngroups, &
    !$omp& ys_hier_group_id, ys_reduced_matrix_lu, ys_reduced_rhs, bw) &
    !$omp private(iline, iblock, row0, offset)
    do iline = 1, ys_hier_node_owned_nlines
      ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline) = (0.0d0, 0.0d0)
      ys_reduced_rhs(1:4*ys_hier_ngroups, iline) = (0.0d0, 0.0d0)
      do iblock = 0, ys_hier_ngroups - 1
        offset = iblock*(20*ys_hier_node_owned_nlines) + (iline - 1)*20
        row0 = 4*iblock
        ys_reduced_matrix_lu(row0 + 1, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 2, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 3, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 4, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_rhs(row0 + 1, iline) = ycomm_recvbuf(offset + 1)
        ys_reduced_rhs(row0 + 2, iline) = ycomm_recvbuf(offset + 2)
        ys_reduced_rhs(row0 + 3, iline) = ycomm_recvbuf(offset + 3)
        ys_reduced_rhs(row0 + 4, iline) = ycomm_recvbuf(offset + 4)
        if (iblock > 0) then
          ys_reduced_matrix_lu(row0 + 1, bw - 1, iline) = -ycomm_recvbuf(offset + 5)
          ys_reduced_matrix_lu(row0 + 2, bw - 2, iline) = -ycomm_recvbuf(offset + 6)
          ys_reduced_matrix_lu(row0 + 3, bw - 3, iline) = -ycomm_recvbuf(offset + 7)
          ys_reduced_matrix_lu(row0 + 4, bw - 4, iline) = -ycomm_recvbuf(offset + 8)
          ys_reduced_matrix_lu(row0 + 1, bw, iline) = -ycomm_recvbuf(offset + 9)
          ys_reduced_matrix_lu(row0 + 2, bw - 1, iline) = -ycomm_recvbuf(offset + 10)
          ys_reduced_matrix_lu(row0 + 3, bw - 2, iline) = -ycomm_recvbuf(offset + 11)
          ys_reduced_matrix_lu(row0 + 4, bw - 3, iline) = -ycomm_recvbuf(offset + 12)
        end if
        if (iblock < ys_hier_ngroups - 1) then
          ys_reduced_matrix_lu(row0 + 1, bw + 5, iline) = -ycomm_recvbuf(offset + 13)
          ys_reduced_matrix_lu(row0 + 2, bw + 4, iline) = -ycomm_recvbuf(offset + 14)
          ys_reduced_matrix_lu(row0 + 3, bw + 3, iline) = -ycomm_recvbuf(offset + 15)
          ys_reduced_matrix_lu(row0 + 4, bw + 2, iline) = -ycomm_recvbuf(offset + 16)
          ys_reduced_matrix_lu(row0 + 1, bw + 6, iline) = -ycomm_recvbuf(offset + 17)
          ys_reduced_matrix_lu(row0 + 2, bw + 5, iline) = -ycomm_recvbuf(offset + 18)
          ys_reduced_matrix_lu(row0 + 3, bw + 4, iline) = -ycomm_recvbuf(offset + 19)
          ys_reduced_matrix_lu(row0 + 4, bw + 3, iline) = -ycomm_recvbuf(offset + 20)
        end if
      end do

      call ys_factor_banded_complex(ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))
      call ys_solve_factored_banded_complex(ys_reduced_rhs(1:4*ys_hier_ngroups, iline), &
                                            ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))

      row0 = 4*ys_hier_group_id
      ys_hier_group_values(1, iline) = ys_reduced_rhs(row0 + 1, iline)
      ys_hier_group_values(2, iline) = ys_reduced_rhs(row0 + 2, iline)
      ys_hier_group_values(3, iline) = ys_reduced_rhs(row0 + 3, iline)
      ys_hier_group_values(4, iline) = ys_reduced_rhs(row0 + 4, iline)
      ys_hier_group_values(5, iline) = (0.0d0, 0.0d0)
      ys_hier_group_values(6, iline) = (0.0d0, 0.0d0)
      ys_hier_group_values(7, iline) = (0.0d0, 0.0d0)
      ys_hier_group_values(8, iline) = (0.0d0, 0.0d0)
      if (ys_hier_group_id > 0) then
        ys_hier_group_values(5, iline) = ys_reduced_rhs(row0 - 1, iline)
        ys_hier_group_values(6, iline) = ys_reduced_rhs(row0, iline)
      end if
      if (ys_hier_group_id < ys_hier_ngroups - 1) then
        ys_hier_group_values(7, iline) = ys_reduced_rhs(row0 + 5, iline)
        ys_hier_group_values(8, iline) = ys_reduced_rhs(row0 + 6, iline)
      end if
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_hier_global_group_solve")
  end subroutine ys_hier_solve_global_groups_column

  subroutine ys_hier_solve_global_groups_column_alltoall()
    implicit none
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: iline, iblock, row0, offset

    call roctxPush("ys_hier_global_group_solve")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_recvbuf, ycomm_sendbuf, ys_hier_column_global_recv_displs, &
    !$omp& ys_hier_column_return_send_displs, ys_hier_column_solve_nlines, ys_hier_ngroups, &
    !$omp& ys_reduced_matrix_lu, ys_reduced_rhs, bw) &
    !$omp private(iline, iblock, row0, offset)
    do iline = 1, ys_hier_column_solve_nlines
      ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline) = (0.0d0, 0.0d0)
      ys_reduced_rhs(1:4*ys_hier_ngroups, iline) = (0.0d0, 0.0d0)
      do iblock = 0, ys_hier_ngroups - 1
        offset = ys_hier_column_global_recv_displs(iblock + 1) + (iline - 1)*20
        row0 = 4*iblock
        ys_reduced_matrix_lu(row0 + 1, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 2, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 3, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 4, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_rhs(row0 + 1, iline) = ycomm_recvbuf(offset + 1)
        ys_reduced_rhs(row0 + 2, iline) = ycomm_recvbuf(offset + 2)
        ys_reduced_rhs(row0 + 3, iline) = ycomm_recvbuf(offset + 3)
        ys_reduced_rhs(row0 + 4, iline) = ycomm_recvbuf(offset + 4)
        if (iblock > 0) then
          ys_reduced_matrix_lu(row0 + 1, bw - 1, iline) = -ycomm_recvbuf(offset + 5)
          ys_reduced_matrix_lu(row0 + 2, bw - 2, iline) = -ycomm_recvbuf(offset + 6)
          ys_reduced_matrix_lu(row0 + 3, bw - 3, iline) = -ycomm_recvbuf(offset + 7)
          ys_reduced_matrix_lu(row0 + 4, bw - 4, iline) = -ycomm_recvbuf(offset + 8)
          ys_reduced_matrix_lu(row0 + 1, bw, iline) = -ycomm_recvbuf(offset + 9)
          ys_reduced_matrix_lu(row0 + 2, bw - 1, iline) = -ycomm_recvbuf(offset + 10)
          ys_reduced_matrix_lu(row0 + 3, bw - 2, iline) = -ycomm_recvbuf(offset + 11)
          ys_reduced_matrix_lu(row0 + 4, bw - 3, iline) = -ycomm_recvbuf(offset + 12)
        end if
        if (iblock < ys_hier_ngroups - 1) then
          ys_reduced_matrix_lu(row0 + 1, bw + 5, iline) = -ycomm_recvbuf(offset + 13)
          ys_reduced_matrix_lu(row0 + 2, bw + 4, iline) = -ycomm_recvbuf(offset + 14)
          ys_reduced_matrix_lu(row0 + 3, bw + 3, iline) = -ycomm_recvbuf(offset + 15)
          ys_reduced_matrix_lu(row0 + 4, bw + 2, iline) = -ycomm_recvbuf(offset + 16)
          ys_reduced_matrix_lu(row0 + 1, bw + 6, iline) = -ycomm_recvbuf(offset + 17)
          ys_reduced_matrix_lu(row0 + 2, bw + 5, iline) = -ycomm_recvbuf(offset + 18)
          ys_reduced_matrix_lu(row0 + 3, bw + 4, iline) = -ycomm_recvbuf(offset + 19)
          ys_reduced_matrix_lu(row0 + 4, bw + 3, iline) = -ycomm_recvbuf(offset + 20)
        end if
      end do

      call ys_factor_banded_complex(ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))
      call ys_solve_factored_banded_complex(ys_reduced_rhs(1:4*ys_hier_ngroups, iline), &
                                            ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))

      do iblock = 0, ys_hier_ngroups - 1
        offset = ys_hier_column_return_send_displs(iblock + 1) + (iline - 1)*YS_REDUCED_RETURN_WIDTH
        row0 = 4*iblock
        ycomm_sendbuf(offset + 1) = ys_reduced_rhs(row0 + 1, iline)
        ycomm_sendbuf(offset + 2) = ys_reduced_rhs(row0 + 2, iline)
        ycomm_sendbuf(offset + 3) = ys_reduced_rhs(row0 + 3, iline)
        ycomm_sendbuf(offset + 4) = ys_reduced_rhs(row0 + 4, iline)
        ycomm_sendbuf(offset + 5) = (0.0d0, 0.0d0)
        ycomm_sendbuf(offset + 6) = (0.0d0, 0.0d0)
        ycomm_sendbuf(offset + 7) = (0.0d0, 0.0d0)
        ycomm_sendbuf(offset + 8) = (0.0d0, 0.0d0)
        if (iblock > 0) then
          ycomm_sendbuf(offset + 5) = ys_reduced_rhs(row0 - 1, iline)
          ycomm_sendbuf(offset + 6) = ys_reduced_rhs(row0, iline)
        end if
        if (iblock < ys_hier_ngroups - 1) then
          ycomm_sendbuf(offset + 7) = ys_reduced_rhs(row0 + 5, iline)
          ycomm_sendbuf(offset + 8) = ys_reduced_rhs(row0 + 6, iline)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_hier_global_group_solve")
  end subroutine ys_hier_solve_global_groups_column_alltoall

  subroutine ys_hier_solve_global_groups(nlines)
    implicit none
    integer(C_INT), intent(in) :: nlines
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: iline, global_line, iblock, row0, group_first, group_count, owner_local
    integer(C_INT) :: src_rank, src_first, intersection_first, offset, owner_rank, local_row0
    integer(C_INT) :: base_count, remainder, split_line

    call roctxPush("ys_hier_global_group_solve")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_recvbuf, ycomm_sendbuf, ys_hier_global_recv_displs, ys_hier_group_return_send_displs, &
    !$omp& ys_hier_group_first_by_group, ys_hier_group_count_by_group, ys_hier_rank_line_first, &
    !$omp& ys_hier_global_first_line, ys_hier_global_owned_nlines, ys_hier_ngroups, ys_reduced_matrix_lu, &
    !$omp& ys_reduced_rhs, nlines, bw) &
    !$omp private(iline, global_line, iblock, row0, group_first, group_count, owner_local, src_rank, src_first, &
    !$omp& intersection_first, offset, owner_rank, local_row0, base_count, remainder, split_line)
    do iline = 1, ys_hier_global_owned_nlines
      global_line = ys_hier_global_first_line + iline - 1
      ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline) = (0.0d0, 0.0d0)
      ys_reduced_rhs(1:4*ys_hier_ngroups, iline) = (0.0d0, 0.0d0)
      do iblock = 0, ys_hier_ngroups - 1
        group_first = ys_hier_group_first_by_group(iblock + 1)
        group_count = ys_hier_group_count_by_group(iblock + 1)
        base_count = nlines/group_count
        remainder = mod(nlines, group_count)
        split_line = (base_count + 1_C_INT)*remainder
        if (global_line <= split_line) then
          owner_local = (global_line - 1)/(base_count + 1_C_INT)
        else
          owner_local = remainder + (global_line - split_line - 1_C_INT)/base_count
        end if
        src_rank = group_first + owner_local
        src_first = ys_hier_rank_line_first(src_rank + 1)
        intersection_first = max(ys_hier_global_first_line, src_first)
        offset = ys_hier_global_recv_displs(src_rank + 1) + (global_line - intersection_first)*20
        row0 = 4*iblock
        ys_reduced_matrix_lu(row0 + 1, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 2, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 3, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 4, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_rhs(row0 + 1, iline) = ycomm_recvbuf(offset + 1)
        ys_reduced_rhs(row0 + 2, iline) = ycomm_recvbuf(offset + 2)
        ys_reduced_rhs(row0 + 3, iline) = ycomm_recvbuf(offset + 3)
        ys_reduced_rhs(row0 + 4, iline) = ycomm_recvbuf(offset + 4)
        if (iblock > 0) then
          ys_reduced_matrix_lu(row0 + 1, bw - 1, iline) = -ycomm_recvbuf(offset + 5)
          ys_reduced_matrix_lu(row0 + 2, bw - 2, iline) = -ycomm_recvbuf(offset + 6)
          ys_reduced_matrix_lu(row0 + 3, bw - 3, iline) = -ycomm_recvbuf(offset + 7)
          ys_reduced_matrix_lu(row0 + 4, bw - 4, iline) = -ycomm_recvbuf(offset + 8)
          ys_reduced_matrix_lu(row0 + 1, bw, iline) = -ycomm_recvbuf(offset + 9)
          ys_reduced_matrix_lu(row0 + 2, bw - 1, iline) = -ycomm_recvbuf(offset + 10)
          ys_reduced_matrix_lu(row0 + 3, bw - 2, iline) = -ycomm_recvbuf(offset + 11)
          ys_reduced_matrix_lu(row0 + 4, bw - 3, iline) = -ycomm_recvbuf(offset + 12)
        end if
        if (iblock < ys_hier_ngroups - 1) then
          ys_reduced_matrix_lu(row0 + 1, bw + 5, iline) = -ycomm_recvbuf(offset + 13)
          ys_reduced_matrix_lu(row0 + 2, bw + 4, iline) = -ycomm_recvbuf(offset + 14)
          ys_reduced_matrix_lu(row0 + 3, bw + 3, iline) = -ycomm_recvbuf(offset + 15)
          ys_reduced_matrix_lu(row0 + 4, bw + 2, iline) = -ycomm_recvbuf(offset + 16)
          ys_reduced_matrix_lu(row0 + 1, bw + 6, iline) = -ycomm_recvbuf(offset + 17)
          ys_reduced_matrix_lu(row0 + 2, bw + 5, iline) = -ycomm_recvbuf(offset + 18)
          ys_reduced_matrix_lu(row0 + 3, bw + 4, iline) = -ycomm_recvbuf(offset + 19)
          ys_reduced_matrix_lu(row0 + 4, bw + 3, iline) = -ycomm_recvbuf(offset + 20)
        end if
      end do

      call ys_factor_banded_complex(ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))
      call ys_solve_factored_banded_complex(ys_reduced_rhs(1:4*ys_hier_ngroups, iline), &
                                            ys_reduced_matrix_lu(1:4*ys_hier_ngroups, :, iline))
      do iblock = 0, ys_hier_ngroups - 1
        group_first = ys_hier_group_first_by_group(iblock + 1)
        group_count = ys_hier_group_count_by_group(iblock + 1)
        base_count = nlines/group_count
        remainder = mod(nlines, group_count)
        split_line = (base_count + 1_C_INT)*remainder
        if (global_line <= split_line) then
          owner_local = (global_line - 1)/(base_count + 1_C_INT)
          local_row0 = global_line - owner_local*(base_count + 1_C_INT)
        else
          owner_local = remainder + (global_line - split_line - 1_C_INT)/base_count
          local_row0 = global_line - split_line - (owner_local - remainder)*base_count
        end if
        owner_rank = group_first + owner_local
        src_first = ys_hier_rank_line_first(owner_rank + 1)
        intersection_first = max(ys_hier_global_first_line, src_first)
        offset = ys_hier_group_return_send_displs(owner_rank + 1) + &
                 (global_line - intersection_first)*YS_REDUCED_RETURN_WIDTH
        row0 = 4*iblock
        ycomm_sendbuf(offset + 1) = ys_reduced_rhs(row0 + 1, iline)
        ycomm_sendbuf(offset + 2) = ys_reduced_rhs(row0 + 2, iline)
        ycomm_sendbuf(offset + 3) = ys_reduced_rhs(row0 + 3, iline)
        ycomm_sendbuf(offset + 4) = ys_reduced_rhs(row0 + 4, iline)
        ycomm_sendbuf(offset + 5) = (0.0d0, 0.0d0)
        ycomm_sendbuf(offset + 6) = (0.0d0, 0.0d0)
        ycomm_sendbuf(offset + 7) = (0.0d0, 0.0d0)
        ycomm_sendbuf(offset + 8) = (0.0d0, 0.0d0)
        if (iblock > 0) then
          ycomm_sendbuf(offset + 5) = ys_reduced_rhs(row0 - 1, iline)
          ycomm_sendbuf(offset + 6) = ys_reduced_rhs(row0, iline)
        end if
        if (iblock < ys_hier_ngroups - 1) then
          ycomm_sendbuf(offset + 7) = ys_reduced_rhs(row0 + 5, iline)
          ycomm_sendbuf(offset + 8) = ys_reduced_rhs(row0 + 6, iline)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_hier_global_group_solve")
  end subroutine ys_hier_solve_global_groups

  subroutine ys_hier_recover_leaf_returns_column()
    implicit none
    integer(C_INT) :: iline, child, child_rank, row0, k, offset
    complex(C_DOUBLE_COMPLEX) :: prev1, prev2, next1, next2

    call roctxPush("ys_hier_recover_leaf_returns")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_sendbuf, ys_hier_leaf_return_send_displs, ys_hier_node_owned_nlines, ys_hier_group_first, &
    !$omp& ys_hier_group_count, ys_hier_recover_basis, ys_hier_group_values, ys_reduced_rhs) &
    !$omp private(iline, child, child_rank, row0, k, offset, prev1, prev2, next1, next2)
    do iline = 1, ys_hier_node_owned_nlines
      prev1 = ys_hier_group_values(5, iline)
      prev2 = ys_hier_group_values(6, iline)
      next1 = ys_hier_group_values(7, iline)
      next2 = ys_hier_group_values(8, iline)
      do child = 0, ys_hier_group_count - 1
        row0 = 4*child
        do k = 1, 4
          ys_reduced_rhs(row0 + k, iline) = ys_hier_recover_basis(row0 + k, 1, iline) + &
                                            ys_hier_recover_basis(row0 + k, 2, iline)*prev1 + &
                                            ys_hier_recover_basis(row0 + k, 3, iline)*prev2 + &
                                            ys_hier_recover_basis(row0 + k, 4, iline)*next1 + &
                                            ys_hier_recover_basis(row0 + k, 5, iline)*next2
        end do
      end do
      do child = 0, ys_hier_group_count - 1
        child_rank = ys_hier_group_first + child
        row0 = 4*child
        offset = ys_hier_leaf_return_send_displs(child_rank + 1) + (iline - 1)*YS_REDUCED_RETURN_WIDTH
        ycomm_sendbuf(offset + 1) = ys_reduced_rhs(row0 + 1, iline)
        ycomm_sendbuf(offset + 2) = ys_reduced_rhs(row0 + 2, iline)
        ycomm_sendbuf(offset + 3) = ys_reduced_rhs(row0 + 3, iline)
        ycomm_sendbuf(offset + 4) = ys_reduced_rhs(row0 + 4, iline)
        if (child > 0) then
          ycomm_sendbuf(offset + 5) = ys_reduced_rhs(row0 - 1, iline)
          ycomm_sendbuf(offset + 6) = ys_reduced_rhs(row0, iline)
        else
          ycomm_sendbuf(offset + 5) = ys_hier_group_values(5, iline)
          ycomm_sendbuf(offset + 6) = ys_hier_group_values(6, iline)
        end if
        if (child < ys_hier_group_count - 1) then
          ycomm_sendbuf(offset + 7) = ys_reduced_rhs(row0 + 5, iline)
          ycomm_sendbuf(offset + 8) = ys_reduced_rhs(row0 + 6, iline)
        else
          ycomm_sendbuf(offset + 7) = ys_hier_group_values(7, iline)
          ycomm_sendbuf(offset + 8) = ys_hier_group_values(8, iline)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_hier_recover_leaf_returns")
  end subroutine ys_hier_recover_leaf_returns_column

  subroutine ys_hier_recover_leaf_returns_column_alltoall()
    implicit none
    integer(C_INT) :: iline, child, child_rank, row0, k, offset
    integer(C_INT) :: src_rank, src_first, base_count, remainder, split_line
    complex(C_DOUBLE_COMPLEX) :: prev1, prev2, next1, next2

    call roctxPush("ys_hier_recover_leaf_returns")
    base_count = ys_hier_node_owned_nlines/ys_hier_column_count
    remainder = mod(ys_hier_node_owned_nlines, ys_hier_column_count)
    split_line = (base_count + 1_C_INT)*remainder
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_recvbuf, ycomm_sendbuf, ys_hier_column_return_recv_displs, &
    !$omp& ys_hier_leaf_return_send_displs, ys_hier_node_owned_nlines, ys_hier_group_first, &
    !$omp& ys_hier_group_count, ys_hier_recover_basis, ys_hier_group_values, ys_reduced_rhs, &
    !$omp& base_count, remainder, split_line) &
    !$omp private(iline, child, child_rank, row0, k, offset, src_rank, src_first, prev1, prev2, next1, next2)
    do iline = 1, ys_hier_node_owned_nlines
      if (iline <= split_line) then
        src_rank = (iline - 1)/(base_count + 1_C_INT)
        src_first = src_rank*(base_count + 1_C_INT) + 1
      else
        src_rank = remainder + (iline - split_line - 1_C_INT)/base_count
        src_first = split_line + (src_rank - remainder)*base_count + 1
      end if
      offset = ys_hier_column_return_recv_displs(src_rank + 1) + &
               (iline - src_first)*YS_REDUCED_RETURN_WIDTH
      do k = 1, YS_REDUCED_RETURN_WIDTH
        ys_hier_group_values(k, iline) = ycomm_recvbuf(offset + k)
      end do
      prev1 = ys_hier_group_values(5, iline)
      prev2 = ys_hier_group_values(6, iline)
      next1 = ys_hier_group_values(7, iline)
      next2 = ys_hier_group_values(8, iline)
      do child = 0, ys_hier_group_count - 1
        row0 = 4*child
        do k = 1, 4
          ys_reduced_rhs(row0 + k, iline) = ys_hier_recover_basis(row0 + k, 1, iline) + &
                                            ys_hier_recover_basis(row0 + k, 2, iline)*prev1 + &
                                            ys_hier_recover_basis(row0 + k, 3, iline)*prev2 + &
                                            ys_hier_recover_basis(row0 + k, 4, iline)*next1 + &
                                            ys_hier_recover_basis(row0 + k, 5, iline)*next2
        end do
      end do
      do child = 0, ys_hier_group_count - 1
        child_rank = ys_hier_group_first + child
        row0 = 4*child
        offset = ys_hier_leaf_return_send_displs(child_rank + 1) + (iline - 1)*YS_REDUCED_RETURN_WIDTH
        ycomm_sendbuf(offset + 1) = ys_reduced_rhs(row0 + 1, iline)
        ycomm_sendbuf(offset + 2) = ys_reduced_rhs(row0 + 2, iline)
        ycomm_sendbuf(offset + 3) = ys_reduced_rhs(row0 + 3, iline)
        ycomm_sendbuf(offset + 4) = ys_reduced_rhs(row0 + 4, iline)
        if (child > 0) then
          ycomm_sendbuf(offset + 5) = ys_reduced_rhs(row0 - 1, iline)
          ycomm_sendbuf(offset + 6) = ys_reduced_rhs(row0, iline)
        else
          ycomm_sendbuf(offset + 5) = ys_hier_group_values(5, iline)
          ycomm_sendbuf(offset + 6) = ys_hier_group_values(6, iline)
        end if
        if (child < ys_hier_group_count - 1) then
          ycomm_sendbuf(offset + 7) = ys_reduced_rhs(row0 + 5, iline)
          ycomm_sendbuf(offset + 8) = ys_reduced_rhs(row0 + 6, iline)
        else
          ycomm_sendbuf(offset + 7) = ys_hier_group_values(7, iline)
          ycomm_sendbuf(offset + 8) = ys_hier_group_values(8, iline)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_hier_recover_leaf_returns")
  end subroutine ys_hier_recover_leaf_returns_column_alltoall

  subroutine ys_hier_recover_leaf_returns()
    implicit none
    integer(C_INT) :: iline, global_line, src_rank, src_first, intersection_first, offset
    integer(C_INT) :: child, child_rank, row0, k
    integer(C_INT) :: base_count, remainder, split_line
    complex(C_DOUBLE_COMPLEX) :: prev1, prev2, next1, next2

    call roctxPush("ys_hier_recover_leaf_returns")
    base_count = size(ys_reduced_rows_send, 2)/npy_grid
    remainder = mod(size(ys_reduced_rows_send, 2), npy_grid)
    split_line = (base_count + 1_C_INT)*remainder
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ycomm_recvbuf, ycomm_sendbuf, ys_hier_group_return_recv_displs, ys_hier_leaf_return_send_displs, &
    !$omp& ys_hier_global_line_first, ys_hier_node_first_line, ys_hier_node_owned_nlines, ys_hier_group_first, &
    !$omp& ys_hier_group_count, ys_hier_recover_basis, ys_hier_group_values, ys_reduced_rhs, base_count, remainder, split_line) &
    !$omp private(iline, global_line, src_rank, src_first, intersection_first, offset, child, child_rank, &
    !$omp& row0, k, prev1, prev2, next1, next2)
    do iline = 1, ys_hier_node_owned_nlines
      global_line = ys_hier_node_first_line + iline - 1
      if (global_line <= split_line) then
        src_rank = (global_line - 1)/(base_count + 1_C_INT)
      else
        src_rank = remainder + (global_line - split_line - 1_C_INT)/base_count
      end if
      src_first = ys_hier_global_line_first(src_rank + 1)
      intersection_first = max(ys_hier_node_first_line, src_first)
      offset = ys_hier_group_return_recv_displs(src_rank + 1) + &
               (global_line - intersection_first)*YS_REDUCED_RETURN_WIDTH
      do k = 1, YS_REDUCED_RETURN_WIDTH
        ys_hier_group_values(k, iline) = ycomm_recvbuf(offset + k)
      end do
      prev1 = ys_hier_group_values(5, iline)
      prev2 = ys_hier_group_values(6, iline)
      next1 = ys_hier_group_values(7, iline)
      next2 = ys_hier_group_values(8, iline)
      do child = 0, ys_hier_group_count - 1
        row0 = 4*child
        do k = 1, 4
          ys_reduced_rhs(row0 + k, iline) = ys_hier_recover_basis(row0 + k, 1, iline) + &
                                            ys_hier_recover_basis(row0 + k, 2, iline)*prev1 + &
                                            ys_hier_recover_basis(row0 + k, 3, iline)*prev2 + &
                                            ys_hier_recover_basis(row0 + k, 4, iline)*next1 + &
                                            ys_hier_recover_basis(row0 + k, 5, iline)*next2
        end do
      end do
      do child = 0, ys_hier_group_count - 1
        child_rank = ys_hier_group_first + child
        row0 = 4*child
        offset = ys_hier_leaf_return_send_displs(child_rank + 1) + (iline - 1)*YS_REDUCED_RETURN_WIDTH
        ycomm_sendbuf(offset + 1) = ys_reduced_rhs(row0 + 1, iline)
        ycomm_sendbuf(offset + 2) = ys_reduced_rhs(row0 + 2, iline)
        ycomm_sendbuf(offset + 3) = ys_reduced_rhs(row0 + 3, iline)
        ycomm_sendbuf(offset + 4) = ys_reduced_rhs(row0 + 4, iline)
        if (child > 0) then
          ycomm_sendbuf(offset + 5) = ys_reduced_rhs(row0 - 1, iline)
          ycomm_sendbuf(offset + 6) = ys_reduced_rhs(row0, iline)
        else
          ycomm_sendbuf(offset + 5) = ys_hier_group_values(5, iline)
          ycomm_sendbuf(offset + 6) = ys_hier_group_values(6, iline)
        end if
        if (child < ys_hier_group_count - 1) then
          ycomm_sendbuf(offset + 7) = ys_reduced_rhs(row0 + 5, iline)
          ycomm_sendbuf(offset + 8) = ys_reduced_rhs(row0 + 6, iline)
        else
          ycomm_sendbuf(offset + 7) = ys_hier_group_values(7, iline)
          ycomm_sendbuf(offset + 8) = ys_hier_group_values(8, iline)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_hier_recover_leaf_returns")
  end subroutine ys_hier_recover_leaf_returns

  subroutine ys_hier_unpack_leaf_returns(nlines)
    implicit none
    integer(C_INT), intent(in) :: nlines
    integer(C_INT) :: iline, src_local, src_rank, src_first, offset, row0
    integer(C_INT) :: base_count, remainder, split_line

    row0 = 4*ipy
    base_count = nlines/ys_hier_group_count
    remainder = mod(nlines, ys_hier_group_count)
    split_line = (base_count + 1_C_INT)*remainder
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_reduced_rhs, ys_left_interface_values, ys_right_interface_values, ycomm_recvbuf, &
    !$omp& ys_hier_leaf_return_recv_displs, ys_hier_rank_line_first, ys_hier_group_first, base_count, remainder, &
    !$omp& split_line, row0, nlines) &
    !$omp private(iline, src_local, src_rank, src_first, offset)
    do iline = 1, nlines
      if (iline <= split_line) then
        src_local = (iline - 1)/(base_count + 1_C_INT)
      else
        src_local = remainder + (iline - split_line - 1_C_INT)/base_count
      end if
      src_rank = ys_hier_group_first + src_local
      src_first = ys_hier_rank_line_first(src_rank + 1)
      offset = ys_hier_leaf_return_recv_displs(src_rank + 1) + (iline - src_first)*YS_REDUCED_RETURN_WIDTH
      ys_reduced_rhs(row0 + 1, iline) = ycomm_recvbuf(offset + 1)
      ys_reduced_rhs(row0 + 2, iline) = ycomm_recvbuf(offset + 2)
      ys_reduced_rhs(row0 + 3, iline) = ycomm_recvbuf(offset + 3)
      ys_reduced_rhs(row0 + 4, iline) = ycomm_recvbuf(offset + 4)
      ys_left_interface_values(1, iline) = ycomm_recvbuf(offset + 5)
      ys_left_interface_values(2, iline) = ycomm_recvbuf(offset + 6)
      ys_right_interface_values(1, iline) = ycomm_recvbuf(offset + 7)
      ys_right_interface_values(2, iline) = ycomm_recvbuf(offset + 8)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_hier_unpack_leaf_returns

  subroutine ys_factor_banded_complex(a)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: a(:, :)
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: n, i, j, t, last
    complex(C_DOUBLE_COMPLEX) :: piv, factor

    n = size(a, 1)
    do i = 1, n
      piv = a(i, bw + 1)
      last = min(bw, n - i)
      do j = 1, last
        factor = a(i + j, bw + 1 - j)/piv
        a(i + j, bw + 1 - j) = factor
        do t = 1, last
          a(i + j, bw + 1 + t - j) = a(i + j, bw + 1 + t - j) - factor*a(i, bw + 1 + t)
        end do
      end do
    end do
  end subroutine ys_factor_banded_complex

  subroutine ys_solve_factored_banded_complex(rhs, a)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: a(:, :)
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
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

  subroutine ys_solve_factored_banded_complex_multi(rhs, a)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:, :)
    complex(C_DOUBLE_COMPLEX), intent(in) :: a(:, :)
    integer(C_INT), parameter :: bw = YS_REDUCED_BW
    integer(C_INT) :: n, nrhs, i, j, irhs
    complex(C_DOUBLE_COMPLEX) :: factor, piv

    n = size(a, 1)
    nrhs = size(rhs, 2)
    do i = 1, n
      do j = max(1_C_INT, i - bw), i - 1
        factor = a(i, bw + 1 + j - i)
        do irhs = 1, nrhs
          rhs(i, irhs) = rhs(i, irhs) - factor*rhs(j, irhs)
        end do
      end do
    end do

    do i = n, 1, -1
      do j = i + 1, min(n, i + bw)
        factor = a(i, bw + 1 + j - i)
        do irhs = 1, nrhs
          rhs(i, irhs) = rhs(i, irhs) - factor*rhs(j, irhs)
        end do
      end do
      piv = a(i, bw + 1)
      do irhs = 1, nrhs
        rhs(i, irhs) = rhs(i, irhs)/piv
      end do
    end do
  end subroutine ys_solve_factored_banded_complex_multi

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
