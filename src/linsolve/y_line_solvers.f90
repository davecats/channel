#include "header.h"

module y_line_solvers

  use, intrinsic :: iso_c_binding
  use channel_grid, only: ny0, nyN, nx0, nxN, nxB, npy_grid, ipy, iproc
  ! ycomm_* and ensure_ycomm_buffers exist in every build; only MPI_COMM_Y is
  ! MPI-specific.  This use must stay outside the #ifdef below: with it inside,
  ! nvfortran rejects ycomm_sendbuf/ycomm_recvbuf in the shared() clauses of the
  ! target regions further down ("must appear in a SHARED or PRIVATE clause"),
  ! while gfortran accepts the same source.
  use mpi_transpose, only: ycomm_sendbuf, ycomm_recvbuf, ensure_ycomm_buffers
  use y_schur_solver, only: ys_schur_config, ys_schur_workspace, ys_schur_configure, &
                            ys_schur_prepare, ys_schur_release, ys_schur_solve_from_packed, &
                            YS_SCHUR_EXCHANGE_AUTO
#ifdef HAVE_MPI
  use mpi_transpose, only: MPI_COMM_Y
  use y_pipeline_nccl, only: channel_comm_use_nccl, channel_comm_p2p_ensure, &
                             channel_comm_send, channel_comm_recv, channel_comm_sendrecv, &
                             channel_comm_cache, channel_comm_cache_reserve, channel_comm_cache_finalize
#endif
  use roctx, only: roctxPush, roctxPop
  use byte_workspace, only: workspace_request, workspace_release, workspace_slice, workspace_align_offset
  use env_options, only: env_flag, env_int, env_int64
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

  ! Which wall-normal solver to use, and how many pipeline batches.
  !
  ! mpi_autotune decides this at startup and writes it here; the compact solve
  ! driver reads it.  It lives with the solver rather than with the autotuner
  ! so that the DNS state module does not have to depend on the autotuner just
  ! to learn which solver was chosen.
  integer(C_INT), parameter, public :: Y_SOLVER_SCHUR = 1_C_INT
  integer(C_INT), parameter, public :: Y_SOLVER_PIPELINED_LU = 2_C_INT
  integer(C_INT), save, public :: ys_selected_solver = Y_SOLVER_SCHUR
  integer(C_INT), save, public :: ys_selected_batches = 0_C_INT

  integer(C_INT), parameter :: YS_ENDPOINT_RESPONSE_CONST = 1_C_INT
  integer(C_INT), parameter :: YS_ENDPOINT_RESPONSE_EVEN_Z = 2_C_INT

  public :: ys_prepare_assembled_workspace, ys_release_workspace, ys_get_workspace_bytes
  public :: ys_get_gpusparse_buffer_bytes
  public :: ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs
  public :: ys_eqm1, ys_eq0, ys_eqn, ys_eqnp1
  public :: ys_boundary_lower_rhs0, ys_boundary_upper_rhsn, ys_boundary_lower_eq, ys_boundary_upper_eq
  public :: ys_boundary_lower_rhs0_owner, ys_boundary_upper_rhsn_owner, ys_boundary_lower_eq_owner, ys_boundary_upper_eq_owner
  public :: ys_lower_ghost_owner, ys_lower_boundary_owner, ys_upper_boundary_owner, ys_upper_ghost_owner
  public :: ys_eqm1_owner, ys_eq0_owner, ys_eqn_owner, ys_eqnp1_owner
  public :: ys_gpsv_matrix, ys_gpsv_rhs, ys_gpsv_line_matrix, ys_gpsv_line_rhs
  public :: ys_gpsv_owner_matrix, ys_gpsv_owner_rhs, ys_owner_nz, ys_owner_nx, ys_owner_ix0, ys_owner_ixN
  public :: ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x
  public :: ys_solve_packed_pentadiagonal, ys_solve_endpoint_schur, ys_solve_pipelined_lu
  public :: ys_finalize_nccl_contexts

  integer(C_INT), save :: ys_workspace_ny = -1
  integer(C_INT), save :: ys_workspace_nz = -1
  integer(C_INT), save :: ys_workspace_nx = -1
  integer(C_INT), save :: ys_workspace_nlines = 0
  integer(C_INT), save :: ys_workspace_active_n = 0
  integer(C_INT), save :: ys_workspace_reduced_node_size = -1
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  logical, save :: ys_gpsv_buffer_reported = .false.
#endif
  integer(C_INT), save :: ys_owner_nz = 0
  integer(C_INT), save :: ys_owner_nx = 0
  integer(C_INT), save :: ys_owner_ix0 = 1
  integer(C_INT), save :: ys_owner_ixN = 0
  integer(C_SIZE_T), parameter :: YS_BATCH_DEFAULT_COMPLEX_CAP = 120000000_C_SIZE_T
  integer(C_SIZE_T), save :: ys_batch_complex_cap = YS_BATCH_DEFAULT_COMPLEX_CAP
  logical, save :: ys_batch_complex_cap_initialized = .false.
  integer(C_INT), save :: ys_pipeline_timing_solve_id = 0_C_INT
#ifdef HAVE_MPI
  type(channel_comm_cache), save :: ys_pipeline_nccl_contexts
#endif

  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_gpsv_matrix_store(:) => null(), ys_gpsv_rhs_store(:) => null()
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_gpsv_ds(:), ys_gpsv_dl(:), ys_gpsv_d(:), ys_gpsv_du(:), ys_gpsv_dw(:), ys_gpsv_x(:)
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_gpsv_matrix(:, :, :, :), ys_gpsv_rhs(:, :, :)
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_gpsv_line_matrix(:, :, :), ys_gpsv_line_rhs(:, :)
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_gpsv_owner_matrix(:, :, :, :), ys_gpsv_owner_rhs(:, :, :)

  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_lower_ghost_store(:) => null(), ys_lower_boundary_store(:) => null()
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_upper_boundary_store(:) => null(), ys_upper_ghost_store(:) => null()
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_lower_ghost_rhs(:), ys_lower_boundary_rhs(:)
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_upper_boundary_rhs(:), ys_upper_ghost_rhs(:)
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_lower_ghost_owner(:, :), ys_lower_boundary_owner(:, :)
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_upper_boundary_owner(:, :), ys_upper_ghost_owner(:, :)

  real(C_DOUBLE), pointer, contiguous, save :: ys_eqm1_store(:) => null(), ys_eq0_store(:) => null()
  real(C_DOUBLE), pointer, contiguous, save :: ys_eqn_store(:) => null(), ys_eqnp1_store(:) => null()
  real(C_DOUBLE), pointer, contiguous, save :: ys_eqm1(:, :), ys_eq0(:, :), ys_eqn(:, :), ys_eqnp1(:, :)
  real(C_DOUBLE), pointer, contiguous, save :: ys_eqm1_owner(:, :, :), ys_eq0_owner(:, :, :), ys_eqn_owner(:, :, :), ys_eqnp1_owner(:, :, :)

  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_boundary_lower_rhs0_store(:) => null(), ys_boundary_upper_rhsn_store(:) => null()
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_boundary_lower_rhs0(:), ys_boundary_upper_rhsn(:)
  real(C_DOUBLE), pointer, contiguous, save :: ys_boundary_lower_eq_store(:) => null(), ys_boundary_upper_eq_store(:) => null()
  real(C_DOUBLE), pointer, contiguous, save :: ys_boundary_lower_eq(:, :), ys_boundary_upper_eq(:, :)
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_boundary_lower_rhs0_owner(:, :), ys_boundary_upper_rhsn_owner(:, :)
  real(C_DOUBLE), pointer, contiguous, save :: ys_boundary_lower_eq_owner(:, :, :), ys_boundary_upper_eq_owner(:, :, :)

  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_reduced_rows_send(:, :) => null()
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_left_interface_values(:, :) => null(), ys_right_interface_values(:, :) => null()
  complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_reduced_rhs(:, :) => null()
  integer(C_INT), allocatable, save :: ys_reduced_pass_counts(:)
  integer(C_INT), save :: ys_reduced_exchange_mode = YS_SCHUR_EXCHANGE_AUTO
  type(ys_schur_workspace), save :: ys_reduced_schur_ws

#ifdef HAVE_CUDA
  type(cusparseHandle), save :: ys_gpsv_handle
  logical, save :: ys_gpsv_handle_created = .false.
#elif defined(HAVE_HIP)
  type(c_ptr), save :: ys_gpsv_handle = c_null_ptr
  type(c_ptr), save :: ys_gpsv_buffer = c_null_ptr
  logical, save :: ys_gpsv_handle_created = .false.
#endif
  integer(C_INT), save :: ys_gpsv_n = -1, ys_gpsv_batch = -1
  integer(C_INT64_T), save :: ys_batch_capacity = 0_C_INT64_T
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  integer(C_SIZE_T), save :: ys_gpsv_buffer_peak_bytes = 0_C_SIZE_T
#endif
#if defined(HAVE_CUDA)
  integer(8), save :: ys_gpsv_buffer_size = 0_8
  integer(8), save :: ys_gpsv_buffer_capacity = 0_8
#elif defined(HAVE_HIP)
  integer(C_SIZE_T), save :: ys_gpsv_buffer_size = 0_C_SIZE_T
  integer(C_SIZE_T), save :: ys_gpsv_buffer_capacity = 0_C_SIZE_T
#endif
 complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_batch_ds(:) => null(), ys_batch_dl(:) => null(), ys_batch_d(:) => null()
 complex(C_DOUBLE_COMPLEX), pointer, contiguous, save :: ys_batch_du(:) => null(), ys_batch_dw(:) => null(), ys_batch_x(:) => null()
#if defined(HAVE_CUDA)
  character(c_char), allocatable, target, save :: ys_gpsv_buffer(:)
  ! Keep these in this translation unit: nvlink does not resolve device-side
  ! calls to declare-target routines pulled from another archive member, so
  ! moving the kernels into their own module breaks the NVHPC GPU link.
  !$omp declare target(ys_factor_penta_interleaved)
  !$omp declare target(ys_factor_penta_interleaved_continue)
  !$omp declare target(ys_forward_substitute_penta_interleaved)
  !$omp declare target(ys_forward_substitute_penta_interleaved_continue)
  !$omp declare target(ys_backward_substitute_penta_interleaved)
  !$omp declare target(ys_backward_substitute_penta_interleaved_continue)
  !$omp declare target(ys_solve_factored_penta_interleaved)
#endif
  integer(C_INT), save :: ys_forced_chunk_nx = 0_C_INT
  logical, save :: ys_forced_chunk_nx_initialized = .false.
  logical, save :: ys_force_custom_gpsv = .false.
  logical, save :: ys_force_custom_gpsv_initialized = .false.

contains

  subroutine ys_reset_workspace_state()
    implicit none

    ys_workspace_ny = -1
    ys_workspace_nz = -1
    ys_workspace_nx = -1
    ys_workspace_nlines = 0
    ys_workspace_active_n = 0
    ys_workspace_reduced_node_size = -1
    ys_owner_nz = 0
    ys_owner_nx = 0
    ys_owner_ix0 = 1
    ys_owner_ixN = 0
  end subroutine ys_reset_workspace_state

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

  integer(C_SIZE_T) function ys_get_batch_complex_cap()
    implicit none
    integer(C_INT64_T) :: parsed

    if (.not. ys_batch_complex_cap_initialized) then
      if (env_int64("CHANNEL_YS_BATCH_MAX_COMPLEX", parsed)) then
        if (parsed > 0_C_INT64_T) ys_batch_complex_cap = min(int(parsed, C_SIZE_T), int(huge(0), C_SIZE_T))
      end if
      ys_batch_complex_cap_initialized = .true.
    end if
    ys_get_batch_complex_cap = ys_batch_complex_cap
  end function ys_get_batch_complex_cap

  integer(C_INT) function ys_get_forced_chunk_nx()
    implicit none
    integer(C_INT) :: parsed

    if (.not. ys_forced_chunk_nx_initialized) then
      if (env_int("CHANNEL_YS_CHUNK_NX", parsed)) then
        if (parsed > 0_C_INT) ys_forced_chunk_nx = parsed
      end if
      ys_forced_chunk_nx_initialized = .true.
    end if
    ys_get_forced_chunk_nx = ys_forced_chunk_nx
  end function ys_get_forced_chunk_nx

  integer(C_SIZE_T) function ys_endpoint_batch_capacity(active_n, nlines, nz)
    implicit none
    integer(C_INT), intent(in) :: active_n, nlines, nz
    integer(C_INT) :: forced_chunk_nx
    integer(C_SIZE_T) :: nx_count, chunk_nx, nlines_z, uncapped

    nlines_z = int(2*nz + 1, C_SIZE_T)
    nx_count = int(nlines/(2*nz + 1), C_SIZE_T)
    chunk_nx = nx_count
    forced_chunk_nx = ys_get_forced_chunk_nx()
    if (forced_chunk_nx > 0_C_INT) chunk_nx = min(chunk_nx, int(forced_chunk_nx, C_SIZE_T))
    uncapped = int(active_n, C_SIZE_T)*chunk_nx*(nlines_z + 4_C_SIZE_T*int(nz + 1, C_SIZE_T))
    ys_endpoint_batch_capacity = min(uncapped, ys_get_batch_complex_cap())
  end function ys_endpoint_batch_capacity

  logical function ys_use_custom_gpsv()
    implicit none

    if (.not. ys_force_custom_gpsv_initialized) then
      call env_flag("CHANNEL_YS_FORCE_CUSTOM_GPSV", ys_force_custom_gpsv)
      ys_force_custom_gpsv_initialized = .true.
    end if
    ys_use_custom_gpsv = ys_force_custom_gpsv
  end function ys_use_custom_gpsv

  ! Single description of the y-solver workspace layout.
  !
  ! With bind = .false. the walk only accumulates the byte count; with
  ! bind = .true. it additionally points every view at its slice of the arena
  ! that workspace_request has already handed back.  Both passes run the same
  ! list of slots, which is what keeps the reported size and the actual
  ! bindings from drifting apart.
  subroutine ys_workspace_layout(active_n, nlines, nz, bind, nbytes)
    implicit none
    integer(C_INT), intent(in) :: active_n, nlines, nz
    logical, intent(in) :: bind
    integer(C_SIZE_T), intent(out) :: nbytes
    type(C_PTR) :: ptr
    complex(C_DOUBLE_COMPLEX), pointer :: cbuf(:)
    real(C_DOUBLE), pointer :: rbuf(:)
    integer(C_SIZE_T) :: offset, nall, batch_capacity, n_complex, n_real
    integer(C_SIZE_T) :: complex_bytes, real_bytes

    complex_bytes = int(C_SIZEOF((0.0_C_DOUBLE, 0.0_C_DOUBLE)), C_SIZE_T)
    real_bytes = int(C_SIZEOF(0.0_C_DOUBLE), C_SIZE_T)
    nall = int(nlines, C_SIZE_T)*int(active_n, C_SIZE_T)
    batch_capacity = ys_endpoint_batch_capacity(active_n, nlines, nz)
    offset = 0_C_SIZE_T

    call slot_complex(5_C_SIZE_T*nall, ys_gpsv_matrix_store)
    call slot_complex(nall, ys_gpsv_rhs_store)

    n_complex = int(nlines, C_SIZE_T)
    call slot_complex(n_complex, ys_lower_ghost_store)
    call slot_complex(n_complex, ys_lower_boundary_store)
    call slot_complex(n_complex, ys_upper_boundary_store)
    call slot_complex(n_complex, ys_upper_ghost_store)

    n_real = 5_C_SIZE_T*int(nlines, C_SIZE_T)
    call slot_real(n_real, ys_eqm1_store)
    call slot_real(n_real, ys_eq0_store)
    call slot_real(n_real, ys_eqn_store)
    call slot_real(n_real, ys_eqnp1_store)

    call slot_complex(n_complex, ys_boundary_lower_rhs0_store)
    call slot_complex(n_complex, ys_boundary_upper_rhsn_store)
    n_real = 4_C_SIZE_T*int(nlines, C_SIZE_T)
    call slot_real(n_real, ys_boundary_lower_eq_store)
    call slot_real(n_real, ys_boundary_upper_eq_store)

    call slot_complex(batch_capacity, ys_batch_ds)
    call slot_complex(batch_capacity, ys_batch_dl)
    call slot_complex(batch_capacity, ys_batch_d)
    call slot_complex(batch_capacity, ys_batch_du)
    call slot_complex(batch_capacity, ys_batch_dw)
    call slot_complex(batch_capacity, ys_batch_x)

    if (npy_grid > 1) then
      call slot_complex_2d(20_C_INT, nlines, ys_reduced_rows_send)
      call slot_complex_2d(2_C_INT, nlines, ys_left_interface_values)
      call slot_complex_2d(2_C_INT, nlines, ys_right_interface_values)
      call slot_complex_2d(4_C_INT*npy_grid, nlines, ys_reduced_rhs)
    end if

    nbytes = offset

  contains
    subroutine slot_complex(count, view)
      integer(C_SIZE_T), intent(in) :: count
      complex(C_DOUBLE_COMPLEX), pointer, contiguous, intent(inout) :: view(:)

      if (bind) then
        call workspace_slice(offset, ptr)
        call c_f_pointer(ptr, cbuf, [int(count)])
        view(1:int(count)) => cbuf
      end if
      offset = workspace_align_offset(offset + count*complex_bytes)
    end subroutine slot_complex

    subroutine slot_real(count, view)
      integer(C_SIZE_T), intent(in) :: count
      real(C_DOUBLE), pointer, contiguous, intent(inout) :: view(:)

      if (bind) then
        call workspace_slice(offset, ptr)
        call c_f_pointer(ptr, rbuf, [int(count)])
        view(1:int(count)) => rbuf
      end if
      offset = workspace_align_offset(offset + count*real_bytes)
    end subroutine slot_real

    subroutine slot_complex_2d(nrows, ncols, view)
      integer(C_INT), intent(in) :: nrows, ncols
      complex(C_DOUBLE_COMPLEX), pointer, contiguous, intent(inout) :: view(:, :)
      integer(C_SIZE_T) :: count

      count = int(nrows, C_SIZE_T)*int(ncols, C_SIZE_T)
      if (bind) then
        call workspace_slice(offset, ptr)
        call c_f_pointer(ptr, cbuf, [int(count)])
        view(1:nrows, 1:ncols) => cbuf
      end if
      offset = workspace_align_offset(offset + count*complex_bytes)
    end subroutine slot_complex_2d
  end subroutine ys_workspace_layout

  subroutine ys_bind_workspace_storage(active_n, nlines, nz)
    implicit none
    integer(C_INT), intent(in) :: active_n, nlines, nz
    type(C_PTR) :: base
    integer(C_SIZE_T) :: total_bytes

    ys_batch_capacity = int(ys_endpoint_batch_capacity(active_n, nlines, nz), C_INT64_T)

    call ys_workspace_layout(active_n, nlines, nz, .false., total_bytes)
    call workspace_request(total_bytes, "y_line_solver", base)
    call ys_workspace_layout(active_n, nlines, nz, .true., total_bytes)
  end subroutine ys_bind_workspace_storage

  subroutine ys_get_workspace_bytes(ny, nz, row_start, row_end, line_start, nlines, nbytes)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, row_start, row_end, line_start, nlines
    integer(C_SIZE_T), intent(out) :: nbytes
    integer(C_INT) :: active_n, nlines_z

    active_n = row_end - row_start + 1
    nlines_z = 2*nz + 1
    if (ny < 1) error stop "ys_get_workspace_bytes requires ny >= 1"
    if (active_n < 1) error stop "ys_get_workspace_bytes requires at least one row"
    if (mod(line_start - 1, nlines_z) /= 0) error stop "ys_get_workspace_bytes requires ix-aligned line_start"
    if (mod(nlines, nlines_z) /= 0) error stop "ys_get_workspace_bytes requires full ix columns"

    call ys_workspace_layout(active_n, nlines, nz, .false., nbytes)
  end subroutine ys_get_workspace_bytes

  subroutine ys_get_gpusparse_buffer_bytes(nbytes)
    implicit none
    integer(C_SIZE_T), intent(out) :: nbytes

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    nbytes = ys_gpsv_buffer_peak_bytes
#else
    nbytes = 0_C_SIZE_T
#endif
  end subroutine ys_get_gpusparse_buffer_bytes

 subroutine ys_allocate_workspace(row_start, row_end, line_start, nlines, nz, schur_pass_counts, schur_exchange_mode, prepare_schur)
    implicit none
    integer(C_INT), intent(in) :: row_start, row_end, line_start, nlines, nz
    integer(C_INT), optional, intent(in) :: schur_pass_counts(:)
    integer(C_INT), intent(in) :: schur_exchange_mode
    logical, optional, intent(in) :: prepare_schur
    integer(C_INT) :: active_n
    type(ys_schur_config) :: cfg
    logical :: prepare_schur_value

    active_n = row_end - row_start + 1
    prepare_schur_value = .true.
    if (present(prepare_schur)) prepare_schur_value = prepare_schur

    if (associated(ys_gpsv_rhs_store)) error stop "ys_allocate_workspace called with core workspace already allocated"
    if (associated(ys_batch_ds)) error stop "ys_allocate_workspace called with batch workspace already allocated"
    if (associated(ys_reduced_rows_send)) error stop "ys_allocate_workspace called with reduced workspace already allocated"
    call ys_bind_workspace_storage(active_n, nlines, nz)

    call ys_bind_core_views(row_start, row_end, line_start, nlines)

    !$omp target enter data map(alloc: ys_gpsv_matrix_store, ys_gpsv_rhs_store, &
    !$omp& ys_lower_ghost_store, ys_lower_boundary_store, ys_upper_boundary_store, ys_upper_ghost_store, &
    !$omp& ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store, &
    !$omp& ys_boundary_lower_rhs0_store, ys_boundary_upper_rhsn_store, &
    !$omp& ys_boundary_lower_eq_store, ys_boundary_upper_eq_store)

    !$omp target enter data map(alloc: ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)

    if (npy_grid > 1 .and. prepare_schur_value) then
      if (allocated(ys_reduced_pass_counts)) deallocate (ys_reduced_pass_counts)
      allocate (ys_reduced_pass_counts(size(schur_pass_counts)))
      ys_reduced_pass_counts = schur_pass_counts
      ys_reduced_exchange_mode = schur_exchange_mode

      call ys_schur_configure(cfg, schur_pass_counts, exchange_mode=schur_exchange_mode)
      call ys_schur_prepare(ys_reduced_schur_ws, cfg, nlines, npy_grid)

      !$omp target enter data map(alloc: ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, &
      !$omp& ys_reduced_rhs)
      ys_workspace_reduced_node_size = 0_C_INT
    end if
  end subroutine ys_allocate_workspace

  subroutine ys_prepare_assembled_workspace(ny, nz, row_start, row_end, line_start, nlines, schur_pass_counts, schur_exchange_mode, prepare_schur)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, row_start, row_end, line_start, nlines
    integer(C_INT), optional, intent(in) :: schur_pass_counts(:)
    integer(C_INT), optional, intent(in) :: schur_exchange_mode
    logical, optional, intent(in) :: prepare_schur
    integer(C_INT) :: active_n, nlines_z, wanted_exchange_mode
    logical :: prepare_schur_value

    active_n = row_end - row_start + 1
    nlines_z = 2*nz + 1
    wanted_exchange_mode = YS_SCHUR_EXCHANGE_AUTO
    if (present(schur_exchange_mode)) wanted_exchange_mode = schur_exchange_mode
    prepare_schur_value = .true.
    if (present(prepare_schur)) prepare_schur_value = prepare_schur

    if (active_n < 1) error stop "ys_prepare_assembled_workspace requires at least one row"
    if (mod(line_start - 1, nlines_z) /= 0) error stop "ys_prepare_assembled_workspace requires ix-aligned line_start"
    if (mod(nlines, nlines_z) /= 0) error stop "ys_prepare_assembled_workspace requires full ix columns"
    if (npy_grid > 1 .and. prepare_schur_value) then
      if (npy_grid > 1 .and. active_n < 4) error stop "ys_solve_ghost_field requires at least four y rows per rank"
      if (active_n < 4) error stop "ys_solve_ghost_field requires at least four active y rows"
      if (.not. present(schur_pass_counts)) error stop "distributed y-Schur requires explicit pass config"
    end if

    ys_owner_nz = nlines_z
    ys_owner_nx = nlines/nlines_z
    ys_owner_ix0 = nx0 + (line_start - 1)/nlines_z
    ys_owner_ixN = ys_owner_ix0 + ys_owner_nx - 1

    ys_workspace_ny = ny
    ys_workspace_nz = nz
    ys_workspace_nx = nlines/nlines_z
    ys_workspace_nlines = nlines
    ys_workspace_active_n = active_n

    if (associated(ys_gpsv_rhs_store) .or. associated(ys_reduced_rows_send)) then
      error stop "ys_prepare_assembled_workspace should only be called once after ys_release_workspace"
    end if

call ys_allocate_workspace(row_start, row_end, line_start, nlines, nz, schur_pass_counts, wanted_exchange_mode, prepare_schur_value)
  end subroutine ys_prepare_assembled_workspace

  subroutine ys_release_workspace(finalize_external)
    implicit none
    logical, intent(in), optional :: finalize_external
    logical :: release_byte_workspace
    logical :: release_external

    release_byte_workspace = .false.
    release_external = .false.
    if (present(finalize_external)) release_external = finalize_external

    if (associated(ys_gpsv_rhs_store)) then
      release_byte_workspace = .true.
      !$omp target exit data map(release: ys_gpsv_matrix_store, ys_gpsv_rhs_store, &
      !$omp& ys_lower_ghost_store, ys_lower_boundary_store, ys_upper_boundary_store, ys_upper_ghost_store, &
      !$omp& ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store, &
      !$omp& ys_boundary_lower_rhs0_store, ys_boundary_upper_rhsn_store, &
      !$omp& ys_boundary_lower_eq_store, ys_boundary_upper_eq_store)

      call ys_nullify_core_views()

      nullify (ys_gpsv_matrix_store, ys_gpsv_rhs_store, &
               ys_lower_ghost_store, ys_lower_boundary_store, ys_upper_boundary_store, ys_upper_ghost_store, &
               ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store, &
               ys_boundary_lower_rhs0_store, ys_boundary_upper_rhsn_store, &
               ys_boundary_lower_eq_store, ys_boundary_upper_eq_store)
    end if

    ! The reduced-system views live in the shared byte workspace, so they must be
    ! dropped whenever it is handed back.  The Schur level hierarchy behind them
    ! does not: it is rebuilt only when the decomposition actually changes, which
    ! is why it is torn down for finalization alone.
    if (associated(ys_reduced_rows_send)) then
      !$omp target exit data map(release: ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, &
      !$omp& ys_reduced_rhs)
      nullify (ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values, ys_reduced_rhs)
    end if
    if (release_external) then
      call ys_schur_release(ys_reduced_schur_ws)
      if (allocated(ys_reduced_pass_counts)) deallocate (ys_reduced_pass_counts)
      ys_reduced_exchange_mode = YS_SCHUR_EXCHANGE_AUTO
    end if
    ys_workspace_reduced_node_size = -1

    if (associated(ys_batch_ds)) then
      !$omp target exit data map(release: ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
      nullify (ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x)
    end if
    if (release_byte_workspace) then
      call workspace_release("y_line_solver")
    end if
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    if (release_external) then
      call ys_release_gpusparse_external()
      call ys_destroy_gpusparse_handle()
    end if
#endif
    ys_batch_capacity = 0_C_INT64_T
    if (release_external) then
      ys_gpsv_n = -1
      ys_gpsv_batch = -1
#ifdef HAVE_CUDA
      ys_gpsv_buffer_size = 0_8
      ys_gpsv_buffer_capacity = 0_8
#elif defined(HAVE_HIP)
      ys_gpsv_buffer_size = 0_C_SIZE_T
      ys_gpsv_buffer_capacity = 0_C_SIZE_T
      ys_gpsv_buffer = c_null_ptr
#endif
    end if
    call ys_reset_workspace_state()
  end subroutine ys_release_workspace

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  subroutine ys_release_gpusparse_external()
    implicit none
#ifdef HAVE_CUDA
    if (allocated(ys_gpsv_buffer)) then
      !$omp target exit data map(delete: ys_gpsv_buffer)
      deallocate (ys_gpsv_buffer)
    end if
    ys_gpsv_buffer_capacity = 0_8
#elif defined(HAVE_HIP)
    if (c_associated(ys_gpsv_buffer)) then
      call omp_target_free(ys_gpsv_buffer, omp_get_default_device())
      ys_gpsv_buffer = c_null_ptr
    end if
    ys_gpsv_buffer_capacity = 0_C_SIZE_T
#endif
  end subroutine ys_release_gpusparse_external

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

    ! Small systems are solved inline in ys_solve_interleaved_pentadiagonal, and
    ! ROCm rejects gpsv workspace queries for n < 3.
    if (n < 3) then
      ys_gpsv_n = n
      ys_gpsv_batch = batch_count
      return
    end if

    if (ys_gpsv_n == n .and. ys_gpsv_batch == batch_count &
#ifdef HAVE_CUDA
        .and. allocated(ys_gpsv_buffer) &
#elif defined(HAVE_HIP)
        .and. c_associated(ys_gpsv_buffer) &
#endif
        ) return

    call ys_create_gpusparse_handle()
    call ys_query_gpsv_buffer_size(ds, dl, d, du, dw, x, n, batch_count, buffer_size)

    ys_gpsv_buffer_size = buffer_size
    ys_gpsv_buffer_peak_bytes = max(ys_gpsv_buffer_peak_bytes, int(buffer_size, C_SIZE_T))
    if (iproc == 0 .and. .not. ys_gpsv_buffer_reported) then
      print '(A,F12.3,A)', "Sparse gpsv external buffer: exact queried=", &
        real(buffer_size, C_DOUBLE)/(1024.0d0*1024.0d0), &
        " MiB (outside shared byte workspace)"
      ys_gpsv_buffer_reported = .true.
    end if
#ifdef HAVE_CUDA
    if (.not. allocated(ys_gpsv_buffer) .or. ys_gpsv_buffer_capacity < buffer_size) then
      call ys_release_gpusparse_external()
      allocate (ys_gpsv_buffer(max(1_8, buffer_size)))
      !$omp target enter data map(alloc: ys_gpsv_buffer)
      ys_gpsv_buffer_capacity = buffer_size
    end if
#elif defined(HAVE_HIP)
    ! On MI300A, hipSPARSE gpsv rejects a workspace buffer created by mapping a
    ! Fortran character array with OpenMP target data, even though the operand
    ! arrays from use_device_addr() are accepted. omp_target_alloc() and hipMalloc()
    ! both work for this buffer; use the OpenMP allocator here to match mpi_transpose.
    if (.not. c_associated(ys_gpsv_buffer) .or. ys_gpsv_buffer_capacity < buffer_size) then
      call ys_release_gpusparse_external()
      ys_gpsv_buffer = omp_target_alloc(max(1_C_SIZE_T, buffer_size), omp_get_default_device())
      if (.not. c_associated(ys_gpsv_buffer)) then
        print *, "OpenMP target allocation failed in ys_prepare_gpusparse_workspace"
        error stop
      end if
      ys_gpsv_buffer_capacity = buffer_size
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

  subroutine ys_solve_interleaved_pentadiagonal(ds, dl, d, du, dw, x, n, batch_count, label)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout), target :: ds(:), dl(:), d(:), du(:), dw(:), x(:)
    integer(C_INT), intent(in) :: n, batch_count
    character(*), intent(in) :: label
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    integer(C_INT) :: iline
    integer(C_INT64_T) :: p, p1
    complex(C_DOUBLE_COMPLEX) :: d0inv, d1inv, lfac, rhs0, rhs1
#endif

    call roctxPush(label)
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    ! ROCSparse fails if the system size is 1: https://github.com/ROCm/rocSPARSE/blob/develop_deprecated/library/src/precond/rocsparse_gtsv.cpp
    if (ys_use_custom_gpsv()) then
      call ys_solve_custom_interleaved_pentadiagonal(ds, dl, d, du, dw, x, n, batch_count)
    else if (n < 3) then
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(ds, dl, d, du, dw, x, n, batch_count) private(iline, p, p1, d0inv, d1inv, lfac, rhs0, rhs1)
      do iline = 1, batch_count
        p = int(iline, C_INT64_T)
        if (n == 1) then
          d0inv = 1.0d0/d(p)
          d(p) = d0inv
          x(p) = x(p)*d0inv
        else
          p1 = p + int(batch_count, C_INT64_T)
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
    call ys_solve_custom_interleaved_pentadiagonal(ds, dl, d, du, dw, x, n, batch_count)
#endif
    call roctxPop(label)
  end subroutine ys_solve_interleaved_pentadiagonal

  subroutine ys_solve_custom_interleaved_pentadiagonal(ds, dl, d, du, dw, x, n, batch_count)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:), dw(:), x(:)
    integer(C_INT), intent(in) :: n, batch_count
    integer(C_INT) :: iline

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ds, dl, d, du, dw, x, n, batch_count) private(iline)
    do iline = 1, batch_count
      call ys_factor_penta_interleaved(ds, dl, d, du, dw, batch_count, iline, n, .false.)
      call ys_solve_factored_penta_interleaved(x, ds, dl, d, du, dw, batch_count, iline, n)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_solve_custom_interleaved_pentadiagonal

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
    integer(C_INT) :: row_start, active_n, nlines, nlines_z, nz, nx_count, dst_row_base, response_mode
    integer(C_INT) :: reduced_row0
    integer(C_INT) :: nI, exposed_n, left_count, interior_base
    integer(C_INT) :: max_batch_count, max_chunk_nx, forced_chunk_nx, chunk_x0, chunk_nx
    logical :: has_left_interface, has_right_interface, has_padded_dst, single_chunk

    row_start = ny0
    active_n = ys_workspace_active_n
    nlines = ys_workspace_nlines
    nlines_z = size(dst, 2)
    nz = (nlines_z - 1)/2
    nx_count = nlines/nlines_z
    has_padded_dst = (size(dst, 1) == active_n + 4)
    dst_row_base = merge(3_C_INT, 1_C_INT, has_padded_dst)
    response_mode = merge(YS_ENDPOINT_RESPONSE_EVEN_Z, YS_ENDPOINT_RESPONSE_CONST, symmetric_operator)

    reduced_row0 = 4*ipy

    has_left_interface = (ipy > 0)
    has_right_interface = (ipy < npy_grid - 1)
    left_count = merge(2_C_INT, 0_C_INT, has_left_interface)
    exposed_n = left_count + merge(2_C_INT, 0_C_INT, has_right_interface)
    interior_base = left_count
    nI = active_n - exposed_n
    if (nI <= 0) error stop "endpoint Schur solve needs at least one interior row"

    max_batch_count = int(min(int(huge(0_C_INT), C_INT64_T), ys_batch_capacity/int(nI, C_INT64_T)), C_INT)
    select case (response_mode)
    case (YS_ENDPOINT_RESPONSE_CONST)
      max_chunk_nx = (max_batch_count - exposed_n)/nlines_z
    case default
      max_chunk_nx = max_batch_count/(nlines_z + exposed_n*(nz + 1_C_INT))
    end select
    if (max_chunk_nx < 1_C_INT) then
      error stop "endpoint Schur batch workspace cannot hold one x column"
    end if
    max_chunk_nx = min(max_chunk_nx, nx_count)
    forced_chunk_nx = ys_get_forced_chunk_nx()
    if (forced_chunk_nx > 0_C_INT) max_chunk_nx = min(max_chunk_nx, forced_chunk_nx)
    single_chunk = (max_chunk_nx >= nx_count)

    chunk_x0 = 0_C_INT
    do while (chunk_x0 < nx_count)
      chunk_nx = min(max_chunk_nx, nx_count - chunk_x0)
      call solve_endpoint_chunk(chunk_x0, chunk_nx)
      if (exposed_n > 0) then
        call pack_schur_chunk(chunk_x0, chunk_nx)
      else
        call reconstruct_chunk(chunk_x0, chunk_nx)
      end if
      chunk_x0 = chunk_x0 + chunk_nx
    end do

    if (exposed_n > 0) then
      call ys_solve_reduced_interfaces()

      if (single_chunk) then
        call reconstruct_chunk(0_C_INT, nx_count)
      else
        chunk_x0 = 0_C_INT
        do while (chunk_x0 < nx_count)
          chunk_nx = min(max_chunk_nx, nx_count - chunk_x0)
          call solve_endpoint_chunk(chunk_x0, chunk_nx)
          call reconstruct_chunk(chunk_x0, chunk_nx)
          chunk_x0 = chunk_x0 + chunk_nx
        end do
      end if
    end if

  contains

    subroutine endpoint_chunk_shape(ix0_chunk, nx_chunk, line_first, line_count, response_count, chunk_batch_count)
      implicit none
      integer(C_INT), intent(in) :: ix0_chunk, nx_chunk
      integer(C_INT), intent(out) :: line_first, line_count, response_count, chunk_batch_count

      line_first = ix0_chunk*nlines_z + 1_C_INT
      line_count = nx_chunk*nlines_z
      select case (response_mode)
      case (YS_ENDPOINT_RESPONSE_CONST)
        response_count = 1_C_INT
      case default
        response_count = nx_chunk*(nz + 1_C_INT)
      end select
      chunk_batch_count = line_count + exposed_n*response_count
      if (int(nI, C_INT64_T)*int(chunk_batch_count, C_INT64_T) > ys_batch_capacity) then
        error stop "endpoint Schur chunk exceeds batch workspace capacity"
      end if
    end subroutine endpoint_chunk_shape

    subroutine solve_endpoint_chunk(ix0_chunk, nx_chunk)
      implicit none
      integer(C_INT), intent(in) :: ix0_chunk, nx_chunk
      complex(C_DOUBLE_COMPLEX) :: rhs_value, coeff, row_coeffs(-2:2)
      integer(C_INT) :: line_first, line_count, response_count, chunk_batch_count
      integer(C_INT) :: sys, iline, ref_iline, resp, local_i, local_idx, col, coupled_row
      integer(C_INT) :: j, offset, response_slot, exposed_slot
      integer(C_INT) :: ix_local, abs_iz
      integer(C_INT64_T) :: p, actual_p
      logical :: is_actual

      call endpoint_chunk_shape(ix0_chunk, nx_chunk, line_first, line_count, response_count, chunk_batch_count)

      call roctxPush("ys_endpoint_pack_plus_response")
      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, &
      !$omp& ys_batch_ds, ys_batch_dl, ys_batch_d, ys_batch_du, ys_batch_dw, ys_batch_x, &
      !$omp& active_n, nI, nlines, nlines_z, nz, response_count, chunk_batch_count, response_mode, &
      !$omp& line_first, line_count, ix0_chunk, has_left_interface, has_right_interface, left_count, exposed_n, interior_base) &
      !$omp private(sys, iline, ref_iline, resp, local_i, local_idx, col, coupled_row, p, j, offset, &
      !$omp& actual_p, response_slot, exposed_slot, ix_local, abs_iz, rhs_value, coeff, row_coeffs, is_actual)
      do local_i = 0, nI - 1
        do sys = 1, chunk_batch_count
          is_actual = (sys <= line_count)
          if (is_actual) then
            iline = line_first + sys - 1_C_INT
            ref_iline = iline
            response_slot = 0_C_INT
          else
            resp = mod(sys - line_count - 1_C_INT, response_count) + 1_C_INT
            response_slot = (sys - line_count - 1_C_INT)/response_count + 1_C_INT
            select case (response_mode)
            case (YS_ENDPOINT_RESPONSE_CONST)
              ref_iline = 1_C_INT
            case default
              ix_local = ix0_chunk + (resp - 1_C_INT)/(nz + 1_C_INT)
              abs_iz = mod(resp - 1_C_INT, nz + 1_C_INT)
              ref_iline = ix_local*nlines_z + abs_iz + nz + 1_C_INT
            end select
            iline = ref_iline
          end if

          local_idx = interior_base + local_i
          actual_p = int(local_idx, C_INT64_T)*int(nlines, C_INT64_T) + int(ref_iline, C_INT64_T)
          row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), &
                         ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)
          rhs_value = (0.0d0, 0.0d0)
          if (is_actual) rhs_value = ys_gpsv_x(actual_p)

          p = int(local_i, C_INT64_T)*int(chunk_batch_count, C_INT64_T) + int(sys, C_INT64_T)
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
                                              nI, chunk_batch_count, "ys_endpoint_gpsv_plus_response")
    end subroutine solve_endpoint_chunk

    subroutine pack_schur_chunk(ix0_chunk, nx_chunk)
      implicit none
      integer(C_INT), intent(in) :: ix0_chunk, nx_chunk
      complex(C_DOUBLE_COMPLEX) :: coeff, row_coeffs(-2:2), s4(4, 4), rhs4(4, 5)
      complex(C_DOUBLE_COMPLEX) :: pivot4, factor4
      integer(C_INT) :: line_first, line_count, response_count, chunk_batch_count
      integer(C_INT) :: iline, local_line, ix_local, iz, abs_iz, resp_index, iface, slot, local_idx
      integer(C_INT) :: col, coupled_row, j, exposed_slot, ghost_col, rhs_col, k, m
      integer(C_INT64_T) :: actual_p, batch_p

      if (exposed_n <= 0) return
      call endpoint_chunk_shape(ix0_chunk, nx_chunk, line_first, line_count, response_count, chunk_batch_count)

      call roctxPush("ys_endpoint_pack_schur")
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, ys_batch_x, &
      !$omp& ys_reduced_rows_send, active_n, nI, nlines, nlines_z, nz, response_count, chunk_batch_count, response_mode, &
      !$omp& line_first, line_count, ix0_chunk, has_left_interface, has_right_interface, left_count, exposed_n, interior_base) &
      !$omp private(iline, local_line, ix_local, iz, abs_iz, resp_index, iface, slot, local_idx, actual_p, batch_p, &
      !$omp& row_coeffs, rhs4, s4, col, coeff, coupled_row, j, exposed_slot, ghost_col, &
      !$omp& k, m, rhs_col, pivot4, factor4)
      do local_line = 1, line_count
        iline = line_first + local_line - 1_C_INT
        ix_local = (iline - 1)/nlines_z
        iz = mod(iline - 1, nlines_z) - nz
        select case (response_mode)
        case (YS_ENDPOINT_RESPONSE_CONST)
          resp_index = 1_C_INT
        case default
          abs_iz = abs(iz)
          resp_index = (ix_local - ix0_chunk)*(nz + 1_C_INT) + abs_iz + 1_C_INT
        end select

        ys_reduced_rows_send(:, iline) = (0.0d0, 0.0d0)
        s4(:, :) = (0.0d0, 0.0d0)
        rhs4(:, :) = (0.0d0, 0.0d0)

        do iface = 1, exposed_n
          slot = iface
          if (.not. has_left_interface) slot = iface + 2
          local_idx = merge(slot - 1, active_n + slot - 5, slot <= 2)
          actual_p = int(local_idx, C_INT64_T)*int(nlines, C_INT64_T) + int(iline, C_INT64_T)
          row_coeffs = (/ys_gpsv_ds(actual_p), ys_gpsv_dl(actual_p), ys_gpsv_d(actual_p), &
                         ys_gpsv_du(actual_p), ys_gpsv_dw(actual_p)/)
          rhs4(iface, 1) = ys_gpsv_x(actual_p)

          do col = -2, 2
            coeff = row_coeffs(col)
            if (coeff == (0.0d0, 0.0d0)) cycle
            coupled_row = local_idx + col
            if (coupled_row >= interior_base .and. coupled_row <= interior_base + nI - 1) then
              j = coupled_row - interior_base
              batch_p = int(j, C_INT64_T)*int(chunk_batch_count, C_INT64_T) + int(local_line, C_INT64_T)
              rhs4(iface, 1) = rhs4(iface, 1) - coeff*ys_batch_x(batch_p)
              do exposed_slot = 1, exposed_n
                batch_p = int(j, C_INT64_T)*int(chunk_batch_count, C_INT64_T) + int(line_count, C_INT64_T) + &
                          int(exposed_slot - 1_C_INT, C_INT64_T)*int(response_count, C_INT64_T) + &
                          int(resp_index, C_INT64_T)
                s4(iface, exposed_slot) = s4(iface, exposed_slot) - &
                                          coeff*ys_batch_x(batch_p)
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
    end subroutine pack_schur_chunk

    subroutine reconstruct_chunk(ix0_chunk, nx_chunk)
      implicit none
      integer(C_INT), intent(in) :: ix0_chunk, nx_chunk
      complex(C_DOUBLE_COMPLEX) :: iface_value(4)
      integer(C_INT) :: line_first, line_count, response_count, chunk_batch_count
      integer(C_INT) :: iline, local_line, ix_local, iz_index, iz, abs_iz, resp_index
      integer(C_INT) :: local_i, local_idx, exposed_slot
      integer(C_INT64_T) :: p

      call endpoint_chunk_shape(ix0_chunk, nx_chunk, line_first, line_count, response_count, chunk_batch_count)

      call roctxPush("ys_endpoint_reconstruct")
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(dst, ys_batch_x, ys_reduced_rhs, ys_left_interface_values, ys_right_interface_values, active_n, nI, nlines, &
      !$omp& nlines_z, nz, response_count, chunk_batch_count, response_mode, dst_row_base, has_padded_dst, &
      !$omp& line_first, line_count, ix0_chunk, has_left_interface, has_right_interface, left_count, exposed_n, interior_base, reduced_row0) &
      !$omp private(iline, local_line, ix_local, iz_index, iz, abs_iz, resp_index, iface_value, local_i, local_idx, p, exposed_slot)
      do local_line = 1, line_count
        iline = line_first + local_line - 1_C_INT
        ix_local = (iline - 1)/nlines_z
        iz_index = mod(iline - 1, nlines_z) + 1
        iz = iz_index - nz - 1
        select case (response_mode)
        case (YS_ENDPOINT_RESPONSE_CONST)
          resp_index = 1_C_INT
        case default
          abs_iz = abs(iz)
          resp_index = (ix_local - ix0_chunk)*(nz + 1_C_INT) + abs_iz + 1_C_INT
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
          p = int(local_i, C_INT64_T)*int(chunk_batch_count, C_INT64_T) + int(local_line, C_INT64_T)
          dst(local_idx + dst_row_base, iz_index, ix_local + 1) = ys_batch_x(p)
          do exposed_slot = 1, exposed_n
            p = int(local_i, C_INT64_T)*int(chunk_batch_count, C_INT64_T) + int(line_count, C_INT64_T) + &
                int(exposed_slot - 1_C_INT, C_INT64_T)*int(response_count, C_INT64_T) + int(resp_index, C_INT64_T)
            dst(local_idx + dst_row_base, iz_index, ix_local + 1) = &
              dst(local_idx + dst_row_base, iz_index, ix_local + 1) - &
              ys_batch_x(p)*iface_value(exposed_slot)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
      call roctxPop("ys_endpoint_reconstruct")
    end subroutine reconstruct_chunk
  end subroutine ys_solve_endpoint_schur

  subroutine ys_solve_pipelined_lu(dst, requested_batches)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)
    integer(C_INT), optional, intent(in) :: requested_batches
    integer(C_INT) :: active_n, nlines, nlines_z, nx_count, dst_row_base
    integer(C_INT) :: nbatches, batch_max_lines

    active_n = ys_workspace_active_n
    nlines = ys_workspace_nlines
    nlines_z = size(dst, 2, kind=C_INT)
    nx_count = nlines/nlines_z
    dst_row_base = merge(3_C_INT, 1_C_INT, size(dst, 1, kind=C_INT) == active_n + 4_C_INT)

    if (active_n < 1_C_INT) error stop "pipelined LU solve requires active rows"
    if (nlines < 1_C_INT) error stop "pipelined LU solve requires active lines"

    if (present(requested_batches)) then
      if (requested_batches < 1_C_INT) error stop "pipelined LU requested_batches must be >= 1"
      nbatches = min(requested_batches, nlines)
    else
      nbatches = ys_pipeline_batch_count(nlines)
    end if
    batch_max_lines = (nlines + nbatches - 1_C_INT)/nbatches

    if (npy_grid == 1_C_INT) then
      call ys_solve_packed_pentadiagonal(active_n, nlines, "ys_pipelined_lu_local")
    else
#ifdef HAVE_MPI
      call ys_solve_pipelined_lu_distributed(active_n, nlines, nbatches, batch_max_lines)
#else
      error stop "pipelined LU distributed solve requires MPI"
#endif
    end if

    call ys_unpack_pipelined_solution(dst, active_n, nlines, nlines_z, nx_count, dst_row_base)
  end subroutine ys_solve_pipelined_lu

  integer(C_INT) function ys_pipeline_batch_count(nlines)
    implicit none
    integer(C_INT), intent(in) :: nlines
    integer(C_INT), parameter :: CANDIDATES(5) = [1_C_INT, 2_C_INT, 4_C_INT, 8_C_INT, 16_C_INT]
    integer(C_INT) :: requested
    integer :: i

    if (env_int("CHANNEL_Y_PIPELINE_BATCHES", requested)) then
      if (requested < 1_C_INT) error stop "CHANNEL_Y_PIPELINE_BATCHES must be >= 1"
      ys_pipeline_batch_count = min(requested, nlines)
      return
    end if

    ! Otherwise take the largest candidate that still gives every batch a line.
    ys_pipeline_batch_count = 1_C_INT
    do i = 1, size(CANDIDATES)
      if (CANDIDATES(i) <= nlines) ys_pipeline_batch_count = CANDIDATES(i)
    end do
  end function ys_pipeline_batch_count

  subroutine ys_pipeline_batch_range(batch, nbatches, nlines, first_line, line_count)
    implicit none
    integer(C_INT), intent(in) :: batch, nbatches, nlines
    integer(C_INT), intent(out) :: first_line, line_count
    integer(C_INT) :: base_count, remainder

    base_count = nlines/nbatches
    remainder = mod(nlines, nbatches)
    line_count = base_count
    if (batch <= remainder) line_count = line_count + 1_C_INT
    first_line = (batch - 1_C_INT)*base_count + min(batch - 1_C_INT, remainder) + 1_C_INT
  end subroutine ys_pipeline_batch_range

#ifdef HAVE_MPI
  subroutine ys_solve_pipelined_lu_distributed(active_n, nlines, nbatches, batch_max_lines)
    implicit none
    integer(C_INT), intent(in) :: active_n, nlines, nbatches, batch_max_lines
    ! Buckets of the timing array, in the order of the CSV header below.
    integer, parameter :: T_TOTAL = 1, T_FACTOR_WAIT = 2, T_FORWARD_WAIT = 3, T_BACKWARD_WAIT = 4
    integer, parameter :: T_FACTOR_KERNEL = 5, T_FORWARD_KERNEL = 6, T_BACKWARD_KERNEL = 7
    integer, parameter :: T_PACK = 8, T_COMM_POST = 9, T_HALO = 10
    integer(C_INT), parameter :: TAG_FACTOR = 8410_C_INT
    integer(C_INT), parameter :: TAG_FORWARD = 8510_C_INT
    integer(C_INT), parameter :: TAG_BACKWARD = 8610_C_INT
    integer(C_INT) :: factor_stride, solve_stride, factor_total, solve_total, send_elems, recv_elems
    integer(C_INT) :: batch, first_line, line_count, factor_count, solve_count
    integer(C_INT) :: factor_offset, forward_offset, backward_offset
    integer(C_INT) :: next_backward
    integer :: ierr
    type(MPI_Request), allocatable :: factor_send_req(:), forward_send_req(:), backward_recv_req(:), backward_send_req(:)
    type(MPI_Request) :: factor_recv_req, forward_recv_req
    type(MPI_Status) :: status
    logical :: timing_enabled
    real(C_DOUBLE) :: timing(10), time_start, time_t0

    if (active_n < 2_C_INT) error stop "pipelined LU distributed solve requires at least two local rows"
    if (channel_comm_use_nccl()) then
      call ys_solve_pipelined_lu_nccl_distributed(active_n, nlines, nbatches, batch_max_lines)
      return
    end if

    timing_enabled = .false.
    call env_flag("CHANNEL_Y_PIPELINE_TIMING", timing_enabled)
    timing = 0.0_C_DOUBLE
    if (timing_enabled) then
      ys_pipeline_timing_solve_id = ys_pipeline_timing_solve_id + 1_C_INT
      time_start = MPI_Wtime()
    end if

    factor_stride = 6_C_INT*batch_max_lines
    solve_stride = 2_C_INT*batch_max_lines
    factor_total = factor_stride*nbatches
    solve_total = solve_stride*nbatches
    send_elems = factor_total + solve_total
    recv_elems = factor_total + 2_C_INT*solve_total
    call ensure_ycomm_buffers(send_elems, recv_elems)
    allocate (factor_send_req(nbatches), forward_send_req(nbatches), backward_recv_req(nbatches), backward_send_req(nbatches))

    factor_send_req = MPI_REQUEST_NULL
    forward_send_req = MPI_REQUEST_NULL
    backward_recv_req = MPI_REQUEST_NULL
    backward_send_req = MPI_REQUEST_NULL

    if (ipy < npy_grid - 1_C_INT) then
      call roctxPush("ys_pipeline post_backward_recvs")
      do batch = 1_C_INT, nbatches
        call ys_pipeline_batch_range(batch, nbatches, nlines, first_line, line_count)
        solve_count = 2_C_INT*line_count
        backward_offset = factor_total + solve_total + (batch - 1_C_INT)*solve_stride
        !$omp target data use_device_addr(ycomm_recvbuf)
        call tic()
        call MPI_Irecv(ycomm_recvbuf(backward_offset + 1), solve_count, &
                       MPI_DOUBLE_COMPLEX, ipy + 1_C_INT, TAG_BACKWARD + batch, MPI_COMM_Y, backward_recv_req(batch), ierr)
        call toc(T_COMM_POST)
        !$omp end target data
      end do
      call roctxPop("ys_pipeline post_backward_recvs")
    end if

    next_backward = 1_C_INT
    call roctxPush("ys_pipeline forward_sweep")
    do batch = 1_C_INT, nbatches
      call roctxPush("ys_pipeline forward_batch")
      call ys_pipeline_batch_range(batch, nbatches, nlines, first_line, line_count)
      factor_count = 6_C_INT*line_count
      solve_count = 2_C_INT*line_count
      factor_offset = (batch - 1_C_INT)*factor_stride
      forward_offset = factor_total + (batch - 1_C_INT)*solve_stride
      if (ipy > 0_C_INT) then
        !$omp target data use_device_addr(ycomm_recvbuf)
        call tic()
        call MPI_Irecv(ycomm_recvbuf(factor_offset + 1), factor_count, &
                       MPI_DOUBLE_COMPLEX, ipy - 1_C_INT, TAG_FACTOR + batch, MPI_COMM_Y, factor_recv_req, ierr)
        call MPI_Irecv(ycomm_recvbuf(forward_offset + 1), solve_count, &
                       MPI_DOUBLE_COMPLEX, ipy - 1_C_INT, TAG_FORWARD + batch, MPI_COMM_Y, forward_recv_req, ierr)
        call toc(T_COMM_POST)
        !$omp end target data
        call roctxPush("ys_pipeline MPI_Wait factor_recv")
        call tic()
        call MPI_Wait(factor_recv_req, status, ierr)
        call toc(T_FACTOR_WAIT)
        call roctxPop("ys_pipeline MPI_Wait factor_recv")
        call roctxPush("ys_pipeline apply_factor_continuation")
        call tic()
        call ys_apply_factor_continuation(first_line, line_count, active_n, factor_offset)
        call toc(T_FACTOR_KERNEL)
        call roctxPop("ys_pipeline apply_factor_continuation")
      end if
      call roctxPush("ys_pipeline factor_batch")
      call tic()
      call ys_factor_pipelined_batch(first_line, line_count, active_n, ipy < npy_grid - 1_C_INT)
      call toc(T_FACTOR_KERNEL)
      call roctxPop("ys_pipeline factor_batch")
      if (ipy < npy_grid - 1_C_INT) then
        call roctxPush("ys_pipeline pack_factor_state")
        call tic()
        call ys_pack_factor_state(first_line, line_count, active_n, factor_offset)
        call toc(T_PACK)
        call roctxPop("ys_pipeline pack_factor_state")
        call roctxPush("ys_pipeline MPI_Isend factor_state")
        !$omp target data use_device_addr(ycomm_sendbuf)
        call tic()
        call MPI_Isend(ycomm_sendbuf(factor_offset + 1), factor_count, &
                       MPI_DOUBLE_COMPLEX, ipy + 1_C_INT, TAG_FACTOR + batch, MPI_COMM_Y, factor_send_req(batch), ierr)
        call toc(T_COMM_POST)
        !$omp end target data
        call roctxPop("ys_pipeline MPI_Isend factor_state")
      end if
      if (ipy > 0_C_INT) then
        call roctxPush("ys_pipeline MPI_Wait forward_recv")
        call tic()
        call MPI_Wait(forward_recv_req, status, ierr)
        call toc(T_FORWARD_WAIT)
        call roctxPop("ys_pipeline MPI_Wait forward_recv")
        call roctxPush("ys_pipeline forward_continue")
        call tic()
        call ys_forward_pipelined_batch_continue(first_line, line_count, active_n, forward_offset)
        call toc(T_FORWARD_KERNEL)
        call roctxPop("ys_pipeline forward_continue")
      else
        call roctxPush("ys_pipeline forward_first_rank")
        call tic()
        call ys_forward_pipelined_batch(first_line, line_count, active_n)
        call toc(T_FORWARD_KERNEL)
        call roctxPop("ys_pipeline forward_first_rank")
      end if
      if (ipy < npy_grid - 1_C_INT) then
        call roctxPush("ys_pipeline pack_forward_state")
        call tic()
        call ys_pack_forward_state(first_line, line_count, active_n, forward_offset)
        call toc(T_PACK)
        call roctxPop("ys_pipeline pack_forward_state")
        call roctxPush("ys_pipeline MPI_Isend forward_state")
        !$omp target data use_device_addr(ycomm_sendbuf)
        call tic()
        call MPI_Isend(ycomm_sendbuf(forward_offset + 1), solve_count, &
                       MPI_DOUBLE_COMPLEX, ipy + 1_C_INT, TAG_FORWARD + batch, MPI_COMM_Y, forward_send_req(batch), ierr)
        call toc(T_COMM_POST)
        !$omp end target data
        call roctxPop("ys_pipeline MPI_Isend forward_state")
      end if
      call roctxPush("ys_pipeline drain_backward_nonblocking")
      call drain_backward_batches(batch, .false.)
      call roctxPop("ys_pipeline drain_backward_nonblocking")
      call roctxPop("ys_pipeline forward_batch")
    end do
    call roctxPop("ys_pipeline forward_sweep")
    call roctxPush("ys_pipeline drain_backward_blocking")
    call drain_backward_batches(nbatches, .true.)
    call roctxPop("ys_pipeline drain_backward_blocking")

    if (ipy < npy_grid - 1_C_INT) then
      call roctxPush("ys_pipeline MPI_Waitall factor_sends")
      call tic()
      call MPI_Waitall(nbatches, factor_send_req, MPI_STATUSES_IGNORE, ierr)
      call toc(T_COMM_POST)
      call roctxPop("ys_pipeline MPI_Waitall factor_sends")
      call roctxPush("ys_pipeline MPI_Waitall forward_sends")
      call tic()
      call MPI_Waitall(nbatches, forward_send_req, MPI_STATUSES_IGNORE, ierr)
      call toc(T_COMM_POST)
      call roctxPop("ys_pipeline MPI_Waitall forward_sends")
    end if
    if (ipy > 0_C_INT) then
      call roctxPush("ys_pipeline MPI_Waitall backward_sends")
      call tic()
      call MPI_Waitall(nbatches, backward_send_req, MPI_STATUSES_IGNORE, ierr)
      call toc(T_COMM_POST)
      call roctxPop("ys_pipeline MPI_Waitall backward_sends")
    end if

    call roctxPush("ys_pipeline exchange_solution_halos")
    call tic()
    call ys_exchange_pipelined_solution_halos(active_n, nlines)
    call toc(T_HALO)
    if (timing_enabled) then
      timing(T_TOTAL) = MPI_Wtime() - time_start
      call ys_print_pipeline_timing(ys_pipeline_timing_solve_id, active_n, nlines, nbatches, timing)
    end if
    call roctxPop("ys_pipeline exchange_solution_halos")

    deallocate (factor_send_req, forward_send_req, backward_recv_req, backward_send_req)

  contains
    ! Timed regions never nest here, so one start stamp is enough.
    subroutine tic()
      if (timing_enabled) time_t0 = MPI_Wtime()
    end subroutine tic

    subroutine toc(bucket)
      integer, intent(in) :: bucket

      if (timing_enabled) timing(bucket) = timing(bucket) + (MPI_Wtime() - time_t0)
    end subroutine toc

    subroutine ys_print_pipeline_timing(solve_id, active_n_value, nlines_value, nbatches_value, timing_value)
      integer(C_INT), intent(in) :: solve_id, active_n_value, nlines_value, nbatches_value
      real(C_DOUBLE), intent(in) :: timing_value(10)
      integer :: print_rank

      do print_rank = 0, int(npy_grid) - 1
        call MPI_Barrier(MPI_COMM_Y, ierr)
        if (int(ipy) == print_rank) then
          if (ipy == 0_C_INT) then
            write (*, '(a)') &
              "Y_PIPELINE_TIMING solve, ipy, active_n, nlines, batches, total_us, factor_wait_us, forward_wait_us, backward_wait_us, factor_kernel_us, forward_kernel_us, backward_kernel_us, pack_us, send_post_wait_us, halo_us"
          end if
          write (*, '(i0,", ",i0,", ",i0,", ",i0,", ",i0,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3)') &
            solve_id, ipy, active_n_value, nlines_value, nbatches_value, &
            timing_value(1)*1.0d6, timing_value(2)*1.0d6, timing_value(3)*1.0d6, &
            timing_value(4)*1.0d6, timing_value(5)*1.0d6, timing_value(6)*1.0d6, &
            timing_value(7)*1.0d6, timing_value(8)*1.0d6, timing_value(9)*1.0d6, &
            timing_value(10)*1.0d6
        end if
      end do
      call MPI_Barrier(MPI_COMM_Y, ierr)
    end subroutine ys_print_pipeline_timing

    subroutine drain_backward_batches(completed_batch, blocking)
      integer(C_INT), intent(in) :: completed_batch
      logical, intent(in) :: blocking
      logical :: ready

      do while (next_backward <= completed_batch)
        call ys_pipeline_batch_range(next_backward, nbatches, nlines, first_line, line_count)
        solve_count = 2_C_INT*line_count
        forward_offset = factor_total + (next_backward - 1_C_INT)*solve_stride
        backward_offset = factor_total + solve_total + (next_backward - 1_C_INT)*solve_stride

        if (ipy < npy_grid - 1_C_INT) then
          if (blocking) then
            call roctxPush("ys_pipeline MPI_Wait backward_recv")
            call tic()
            call MPI_Wait(backward_recv_req(next_backward), status, ierr)
            call toc(T_BACKWARD_WAIT)
            call roctxPop("ys_pipeline MPI_Wait backward_recv")
            ready = .true.
          else
            call roctxPush("ys_pipeline MPI_Test backward_recv")
            call MPI_Test(backward_recv_req(next_backward), ready, status, ierr)
            call roctxPop("ys_pipeline MPI_Test backward_recv")
          end if
          if (.not. ready) exit
          call roctxPush("ys_pipeline backward_continue")
          call tic()
          call ys_backward_pipelined_batch_continue(first_line, line_count, active_n, backward_offset)
          call toc(T_BACKWARD_KERNEL)
          call roctxPop("ys_pipeline backward_continue")
        else
          call roctxPush("ys_pipeline backward_top_rank")
          call tic()
          call ys_backward_pipelined_batch(first_line, line_count, active_n)
          call toc(T_BACKWARD_KERNEL)
          call roctxPop("ys_pipeline backward_top_rank")
        end if

        if (ipy > 0_C_INT) then
          if (ipy < npy_grid - 1_C_INT) then
            call roctxPush("ys_pipeline MPI_Wait forward_send_reuse")
            call tic()
            call MPI_Wait(forward_send_req(next_backward), status, ierr)
            call toc(T_COMM_POST)
            call roctxPop("ys_pipeline MPI_Wait forward_send_reuse")
          end if
          call roctxPush("ys_pipeline pack_backward_state")
          call tic()
          call ys_pack_backward_state(first_line, line_count, forward_offset)
          call toc(T_PACK)
          call roctxPop("ys_pipeline pack_backward_state")
          call roctxPush("ys_pipeline MPI_Isend backward_state")
          !$omp target data use_device_addr(ycomm_sendbuf)
          call tic()
          call MPI_Isend(ycomm_sendbuf(forward_offset + 1), solve_count, &
                         MPI_DOUBLE_COMPLEX, ipy - 1_C_INT, TAG_BACKWARD + next_backward, &
                         MPI_COMM_Y, backward_send_req(next_backward), ierr)
          call toc(T_COMM_POST)
          !$omp end target data
          call roctxPop("ys_pipeline MPI_Isend backward_state")
        end if
        next_backward = next_backward + 1_C_INT
      end do
    end subroutine drain_backward_batches
  end subroutine ys_solve_pipelined_lu_distributed

  subroutine ys_solve_pipelined_lu_nccl_distributed(active_n, nlines, nbatches, batch_max_lines)
    implicit none
    integer(C_INT), intent(in) :: active_n, nlines, nbatches, batch_max_lines
    integer(C_INT) :: factor_stride, solve_stride, factor_total, solve_total, send_elems, recv_elems
    integer(C_INT) :: batch, first_line, line_count, factor_count, solve_count
    integer(C_INT) :: factor_offset, forward_offset, backward_offset
    type(c_ptr) :: nccl_ctx

    factor_stride = 6_C_INT*batch_max_lines
    solve_stride = 2_C_INT*batch_max_lines
    factor_total = factor_stride*nbatches
    solve_total = solve_stride*nbatches
    send_elems = factor_total + solve_total
    recv_elems = factor_total + 2_C_INT*solve_total
    call ensure_ycomm_buffers(send_elems, recv_elems)
    call ys_pipeline_nccl_context(nccl_ctx)

    call roctxPush("ys_pipeline_nccl forward_sweep")
    do batch = 1_C_INT, nbatches
      call ys_pipeline_batch_range(batch, nbatches, nlines, first_line, line_count)
      factor_count = 6_C_INT*line_count
      solve_count = 2_C_INT*line_count
      factor_offset = (batch - 1_C_INT)*factor_stride
      forward_offset = factor_total + (batch - 1_C_INT)*solve_stride

      if (ipy > 0_C_INT) then
        call roctxPush("ys_pipeline_nccl recv_factor")
        !$omp target data use_device_addr(ycomm_recvbuf)
        call channel_comm_recv(nccl_ctx, c_loc(ycomm_recvbuf(factor_offset + 1)), factor_count, ipy - 1_C_INT)
        !$omp end target data
        call roctxPop("ys_pipeline_nccl recv_factor")
        call roctxPush("ys_pipeline_nccl apply_factor_continuation")
        call ys_apply_factor_continuation(first_line, line_count, active_n, factor_offset)
        call roctxPop("ys_pipeline_nccl apply_factor_continuation")
      end if

      call roctxPush("ys_pipeline_nccl factor_batch")
      call ys_factor_pipelined_batch(first_line, line_count, active_n, ipy < npy_grid - 1_C_INT)
      call roctxPop("ys_pipeline_nccl factor_batch")

      if (ipy < npy_grid - 1_C_INT) then
        call roctxPush("ys_pipeline_nccl pack_factor_state")
        call ys_pack_factor_state(first_line, line_count, active_n, factor_offset)
        call roctxPop("ys_pipeline_nccl pack_factor_state")
        call roctxPush("ys_pipeline_nccl send_factor")
        !$omp target data use_device_addr(ycomm_sendbuf)
        call channel_comm_send(nccl_ctx, c_loc(ycomm_sendbuf(factor_offset + 1)), factor_count, ipy + 1_C_INT)
        !$omp end target data
        call roctxPop("ys_pipeline_nccl send_factor")
      end if

      if (ipy > 0_C_INT) then
        call roctxPush("ys_pipeline_nccl recv_forward")
        !$omp target data use_device_addr(ycomm_recvbuf)
        call channel_comm_recv(nccl_ctx, c_loc(ycomm_recvbuf(forward_offset + 1)), solve_count, ipy - 1_C_INT)
        !$omp end target data
        call roctxPop("ys_pipeline_nccl recv_forward")
        call roctxPush("ys_pipeline_nccl forward_continue")
        call ys_forward_pipelined_batch_continue(first_line, line_count, active_n, forward_offset)
        call roctxPop("ys_pipeline_nccl forward_continue")
      else
        call roctxPush("ys_pipeline_nccl forward_first_rank")
        call ys_forward_pipelined_batch(first_line, line_count, active_n)
        call roctxPop("ys_pipeline_nccl forward_first_rank")
      end if

      if (ipy < npy_grid - 1_C_INT) then
        call roctxPush("ys_pipeline_nccl pack_forward_state")
        call ys_pack_forward_state(first_line, line_count, active_n, forward_offset)
        call roctxPop("ys_pipeline_nccl pack_forward_state")
        call roctxPush("ys_pipeline_nccl send_forward")
        !$omp target data use_device_addr(ycomm_sendbuf)
        call channel_comm_send(nccl_ctx, c_loc(ycomm_sendbuf(forward_offset + 1)), solve_count, ipy + 1_C_INT)
        !$omp end target data
        call roctxPop("ys_pipeline_nccl send_forward")
      end if
    end do
    call roctxPop("ys_pipeline_nccl forward_sweep")

    call roctxPush("ys_pipeline_nccl backward_sweep")
    do batch = nbatches, 1_C_INT, -1_C_INT
      call ys_pipeline_batch_range(batch, nbatches, nlines, first_line, line_count)
      solve_count = 2_C_INT*line_count
      forward_offset = factor_total + (batch - 1_C_INT)*solve_stride
      backward_offset = factor_total + solve_total + (batch - 1_C_INT)*solve_stride

      if (ipy < npy_grid - 1_C_INT) then
        call roctxPush("ys_pipeline_nccl recv_backward")
        !$omp target data use_device_addr(ycomm_recvbuf)
        call channel_comm_recv(nccl_ctx, c_loc(ycomm_recvbuf(backward_offset + 1)), solve_count, ipy + 1_C_INT)
        !$omp end target data
        call roctxPop("ys_pipeline_nccl recv_backward")
        call roctxPush("ys_pipeline_nccl backward_continue")
        call ys_backward_pipelined_batch_continue(first_line, line_count, active_n, backward_offset)
        call roctxPop("ys_pipeline_nccl backward_continue")
      else
        call roctxPush("ys_pipeline_nccl backward_top_rank")
        call ys_backward_pipelined_batch(first_line, line_count, active_n)
        call roctxPop("ys_pipeline_nccl backward_top_rank")
      end if

      if (ipy > 0_C_INT) then
        call roctxPush("ys_pipeline_nccl pack_backward_state")
        call ys_pack_backward_state(first_line, line_count, forward_offset)
        call roctxPop("ys_pipeline_nccl pack_backward_state")
        call roctxPush("ys_pipeline_nccl send_backward")
        !$omp target data use_device_addr(ycomm_sendbuf)
        call channel_comm_send(nccl_ctx, c_loc(ycomm_sendbuf(forward_offset + 1)), solve_count, ipy - 1_C_INT)
        !$omp end target data
        call roctxPop("ys_pipeline_nccl send_backward")
      end if
    end do
    call roctxPop("ys_pipeline_nccl backward_sweep")

    call roctxPush("ys_pipeline_nccl exchange_solution_halos")
    call ys_exchange_pipelined_solution_halos_nccl(active_n, nlines, nccl_ctx)
    call roctxPop("ys_pipeline_nccl exchange_solution_halos")
  end subroutine ys_solve_pipelined_lu_nccl_distributed

  subroutine ys_pipeline_nccl_context(nccl_ctx)
    implicit none
    type(c_ptr), intent(out) :: nccl_ctx

    call channel_comm_cache_reserve(ys_pipeline_nccl_contexts, npy_grid)
    call channel_comm_p2p_ensure(MPI_COMM_Y, ys_pipeline_nccl_contexts%by_size(npy_grid))
    nccl_ctx = ys_pipeline_nccl_contexts%by_size(npy_grid)
  end subroutine ys_pipeline_nccl_context

  subroutine ys_finalize_nccl_contexts()
    implicit none
#ifdef HAVE_MPI

    call channel_comm_cache_finalize(ys_pipeline_nccl_contexts)
#endif
  end subroutine ys_finalize_nccl_contexts

  subroutine ys_exchange_pipelined_solution_halos(active_n, nlines)
    implicit none
    integer(C_INT), intent(in) :: active_n, nlines
    integer(C_INT), parameter :: TAG_HALO_LOW = 8710_C_INT
    integer(C_INT), parameter :: TAG_HALO_HIGH = 8711_C_INT
    integer(C_INT) :: ierr, count, lower_offset, upper_offset, nreq
    type(MPI_Request) :: req(4)
    type(MPI_Status) :: statuses(4)

    count = 2_C_INT*nlines
    lower_offset = 0_C_INT
    upper_offset = count
    req = MPI_REQUEST_NULL
    nreq = 0

    if (ipy > 0_C_INT) then
      nreq = nreq + 1
      !$omp target data use_device_addr(ycomm_recvbuf)
      call MPI_Irecv(ycomm_recvbuf(lower_offset + 1), count, MPI_DOUBLE_COMPLEX, &
                     ipy - 1_C_INT, TAG_HALO_HIGH, MPI_COMM_Y, req(nreq), ierr)
      !$omp end target data
      call ys_pack_solution_rows(nlines, first_row=1_C_INT, send_offset=lower_offset)
      nreq = nreq + 1
      !$omp target data use_device_addr(ycomm_sendbuf)
      call MPI_Isend(ycomm_sendbuf(lower_offset + 1), count, MPI_DOUBLE_COMPLEX, &
                     ipy - 1_C_INT, TAG_HALO_LOW, MPI_COMM_Y, req(nreq), ierr)
      !$omp end target data
    end if

    if (ipy < npy_grid - 1_C_INT) then
      nreq = nreq + 1
      !$omp target data use_device_addr(ycomm_recvbuf)
      call MPI_Irecv(ycomm_recvbuf(upper_offset + 1), count, MPI_DOUBLE_COMPLEX, &
                     ipy + 1_C_INT, TAG_HALO_LOW, MPI_COMM_Y, req(nreq), ierr)
      !$omp end target data
      call ys_pack_solution_rows(nlines, first_row=active_n - 1_C_INT, send_offset=upper_offset)
      nreq = nreq + 1
      !$omp target data use_device_addr(ycomm_sendbuf)
      call MPI_Isend(ycomm_sendbuf(upper_offset + 1), count, MPI_DOUBLE_COMPLEX, &
                     ipy + 1_C_INT, TAG_HALO_HIGH, MPI_COMM_Y, req(nreq), ierr)
      !$omp end target data
    end if

    if (nreq > 0) call MPI_Waitall(nreq, req(1:nreq), statuses(1:nreq), ierr)
    if (ipy > 0_C_INT) call ys_unpack_lower_solution_halo(nlines, lower_offset)
    if (ipy < npy_grid - 1_C_INT) call ys_unpack_upper_solution_halo(nlines, upper_offset)
  end subroutine ys_exchange_pipelined_solution_halos

  subroutine ys_exchange_pipelined_solution_halos_nccl(active_n, nlines, nccl_ctx)
    implicit none
    integer(C_INT), intent(in) :: active_n, nlines
    type(c_ptr), intent(in), value :: nccl_ctx
    integer(C_INT) :: count, lower_offset, upper_offset

    count = 2_C_INT*nlines
    lower_offset = 0_C_INT
    upper_offset = count

    if (ipy > 0_C_INT) call ys_pack_solution_rows(nlines, first_row=1_C_INT, send_offset=lower_offset)
    if (ipy < npy_grid - 1_C_INT) call ys_pack_solution_rows(nlines, first_row=active_n - 1_C_INT, send_offset=upper_offset)

    if (ipy > 0_C_INT) then
      call roctxPush("ys_pipeline_nccl halo_lower")
      !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
      call channel_comm_sendrecv(nccl_ctx, c_loc(ycomm_sendbuf(lower_offset + 1)), count, ipy - 1_C_INT, &
                                 c_loc(ycomm_recvbuf(lower_offset + 1)), count, ipy - 1_C_INT)
      !$omp end target data
      call roctxPop("ys_pipeline_nccl halo_lower")
      call ys_unpack_lower_solution_halo(nlines, lower_offset)
    end if

    if (ipy < npy_grid - 1_C_INT) then
      call roctxPush("ys_pipeline_nccl halo_upper")
      !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
      call channel_comm_sendrecv(nccl_ctx, c_loc(ycomm_sendbuf(upper_offset + 1)), count, ipy + 1_C_INT, &
                                 c_loc(ycomm_recvbuf(upper_offset + 1)), count, ipy + 1_C_INT)
      !$omp end target data
      call roctxPop("ys_pipeline_nccl halo_upper")
      call ys_unpack_upper_solution_halo(nlines, upper_offset)
    end if
  end subroutine ys_exchange_pipelined_solution_halos_nccl

  ! Copy the pair of solution rows starting at first_row into the halo send
  ! buffer.  The lower halo sends rows 1 and 2, the upper rows active_n-1 and
  ! active_n; the row the pair starts at is the only difference between them.
  subroutine ys_pack_solution_rows(nlines, first_row, send_offset)
    implicit none
    integer(C_INT), intent(in) :: nlines, first_row, send_offset
    integer(C_INT) :: iline
    integer(C_INT64_T) :: p0, p1, q

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_x, ycomm_sendbuf, nlines, first_row, send_offset) private(iline, p0, p1, q)
    do iline = 1_C_INT, nlines
      p0 = int(iline, C_INT64_T) + int(first_row - 1_C_INT, C_INT64_T)*int(nlines, C_INT64_T)
      p1 = p0 + int(nlines, C_INT64_T)
      q = int(send_offset + 2_C_INT*(iline - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = ys_gpsv_x(p0)
      ycomm_sendbuf(q + 2_C_INT64_T) = ys_gpsv_x(p1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_pack_solution_rows

  ! The unpack pair stays written twice on purpose: the two halves differ by
  ! their destination array, not by an index, and ys_left/right_interface_values
  ! are named separately everywhere else they are read.  Parameterising them
  ! would mean either an array dummy inside a target region or a branch in the
  ! kernel, and neither reads better than these two.
  subroutine ys_unpack_lower_solution_halo(nlines, recv_offset)
    implicit none
    integer(C_INT), intent(in) :: nlines, recv_offset
    integer(C_INT) :: iline
    integer(C_INT64_T) :: q

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_left_interface_values, ycomm_recvbuf, nlines, recv_offset) private(iline, q)
    do iline = 1_C_INT, nlines
      q = int(recv_offset + 2_C_INT*(iline - 1_C_INT), C_INT64_T)
      ys_left_interface_values(1, iline) = ycomm_recvbuf(q + 1_C_INT64_T)
      ys_left_interface_values(2, iline) = ycomm_recvbuf(q + 2_C_INT64_T)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_unpack_lower_solution_halo

  subroutine ys_unpack_upper_solution_halo(nlines, recv_offset)
    implicit none
    integer(C_INT), intent(in) :: nlines, recv_offset
    integer(C_INT) :: iline
    integer(C_INT64_T) :: q

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_right_interface_values, ycomm_recvbuf, nlines, recv_offset) private(iline, q)
    do iline = 1_C_INT, nlines
      q = int(recv_offset + 2_C_INT*(iline - 1_C_INT), C_INT64_T)
      ys_right_interface_values(1, iline) = ycomm_recvbuf(q + 1_C_INT64_T)
      ys_right_interface_values(2, iline) = ycomm_recvbuf(q + 2_C_INT64_T)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_unpack_upper_solution_halo

  subroutine ys_apply_factor_continuation(first_line, line_count, active_n, recv_offset)
    implicit none
    integer(C_INT), intent(in) :: first_line, line_count, active_n, recv_offset
    integer(C_INT) :: local_line, iline, stride
    integer(C_INT64_T) :: q

    stride = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ycomm_recvbuf, &
    !$omp& first_line, line_count, active_n, stride, recv_offset) private(local_line, iline, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      q = int(recv_offset + 6_C_INT*(local_line - 1_C_INT), C_INT64_T)
      call ys_factor_penta_interleaved_continue(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, &
                                                stride, iline, active_n, &
                                                ycomm_recvbuf(q + 1_C_INT64_T), ycomm_recvbuf(q + 2_C_INT64_T), &
                                                ycomm_recvbuf(q + 3_C_INT64_T), ycomm_recvbuf(q + 4_C_INT64_T), &
                                                ycomm_recvbuf(q + 5_C_INT64_T), ycomm_recvbuf(q + 6_C_INT64_T))
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_apply_factor_continuation

  subroutine ys_factor_pipelined_batch(first_line, line_count, active_n, open_right)
    implicit none
    integer(C_INT), intent(in) :: first_line, line_count, active_n
    logical, intent(in) :: open_right
    integer(C_INT) :: local_line, iline, stride

    stride = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, stride, first_line, line_count, active_n, open_right) &
    !$omp private(local_line, iline)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      call ys_factor_penta_interleaved(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, &
                                       stride, iline, active_n, open_right)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_factor_pipelined_batch

  subroutine ys_forward_pipelined_batch(first_line, line_count, active_n)
    implicit none
    integer(C_INT), intent(in) :: first_line, line_count, active_n
    integer(C_INT) :: local_line, iline, stride

    stride = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_x, ys_gpsv_ds, ys_gpsv_dl, stride, first_line, line_count, active_n) private(local_line, iline)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      call ys_forward_substitute_penta_interleaved(ys_gpsv_x, ys_gpsv_ds, ys_gpsv_dl, stride, iline, active_n)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_forward_pipelined_batch

  subroutine ys_forward_pipelined_batch_continue(first_line, line_count, active_n, recv_offset)
    implicit none
    integer(C_INT), intent(in) :: first_line, line_count, active_n, recv_offset
    integer(C_INT) :: local_line, iline, stride
    integer(C_INT64_T) :: q

    stride = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_x, ys_gpsv_ds, ys_gpsv_dl, ycomm_recvbuf, stride, first_line, line_count, active_n, recv_offset) &
    !$omp private(local_line, iline, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      q = int(recv_offset + 2_C_INT*(local_line - 1_C_INT), C_INT64_T)
      call ys_forward_substitute_penta_interleaved_continue(ys_gpsv_x, ys_gpsv_ds, ys_gpsv_dl, stride, iline, active_n, &
                                                            ycomm_recvbuf(q + 1_C_INT64_T), ycomm_recvbuf(q + 2_C_INT64_T))
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_forward_pipelined_batch_continue

  subroutine ys_backward_pipelined_batch(first_line, line_count, active_n)
    implicit none
    integer(C_INT), intent(in) :: first_line, line_count, active_n
    integer(C_INT) :: local_line, iline, stride

    stride = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_x, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, stride, first_line, line_count, active_n) private(local_line, iline)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      call ys_backward_substitute_penta_interleaved(ys_gpsv_x, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, stride, iline, active_n)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_backward_pipelined_batch

  subroutine ys_backward_pipelined_batch_continue(first_line, line_count, active_n, recv_offset)
    implicit none
    integer(C_INT), intent(in) :: first_line, line_count, active_n, recv_offset
    integer(C_INT) :: local_line, iline, stride
    integer(C_INT64_T) :: q

    stride = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_x, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ycomm_recvbuf, &
    !$omp& stride, first_line, line_count, active_n, recv_offset) private(local_line, iline, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      q = int(recv_offset + 2_C_INT*(local_line - 1_C_INT), C_INT64_T)
      call ys_backward_substitute_penta_interleaved_continue(ys_gpsv_x, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, &
                                                             stride, iline, active_n, &
                                                             ycomm_recvbuf(q + 1_C_INT64_T), ycomm_recvbuf(q + 2_C_INT64_T))
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_backward_pipelined_batch_continue

  subroutine ys_pack_factor_state(first_line, line_count, active_n, send_offset)
    implicit none
    integer(C_INT), intent(in) :: first_line, line_count, active_n, send_offset
    integer(C_INT) :: local_line, iline, stride
    integer(C_INT64_T) :: p0, p1, q

    stride = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ycomm_sendbuf, stride, first_line, line_count, active_n, send_offset) &
    !$omp private(local_line, iline, p0, p1, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      p0 = int(iline, C_INT64_T) + int(active_n - 2_C_INT, C_INT64_T)*int(stride, C_INT64_T)
      p1 = p0 + int(stride, C_INT64_T)
      q = int(send_offset + 6_C_INT*(local_line - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = ys_gpsv_d(p0)
      ycomm_sendbuf(q + 2_C_INT64_T) = ys_gpsv_du(p0)
      ycomm_sendbuf(q + 3_C_INT64_T) = ys_gpsv_dw(p0)
      ycomm_sendbuf(q + 4_C_INT64_T) = ys_gpsv_d(p1)
      ycomm_sendbuf(q + 5_C_INT64_T) = ys_gpsv_du(p1)
      ycomm_sendbuf(q + 6_C_INT64_T) = ys_gpsv_dw(p1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_pack_factor_state

  subroutine ys_pack_forward_state(first_line, line_count, active_n, send_offset)
    implicit none
    integer(C_INT), intent(in) :: first_line, line_count, active_n, send_offset
    integer(C_INT) :: local_line, iline, stride
    integer(C_INT64_T) :: p0, p1, q

    stride = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_x, ycomm_sendbuf, stride, first_line, line_count, active_n, send_offset) private(local_line, iline, p0, p1, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      p0 = int(iline, C_INT64_T) + int(active_n - 2_C_INT, C_INT64_T)*int(stride, C_INT64_T)
      p1 = p0 + int(stride, C_INT64_T)
      q = int(send_offset + 2_C_INT*(local_line - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = ys_gpsv_x(p0)
      ycomm_sendbuf(q + 2_C_INT64_T) = ys_gpsv_x(p1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_pack_forward_state

  subroutine ys_pack_backward_state(first_line, line_count, send_offset)
    implicit none
    integer(C_INT), intent(in) :: first_line, line_count, send_offset
    integer(C_INT) :: local_line, iline, stride
    integer(C_INT64_T) :: p0, p1, q

    stride = ys_workspace_nlines

    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_x, ycomm_sendbuf, stride, first_line, line_count, send_offset) private(local_line, iline, p0, p1, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      p0 = int(iline, C_INT64_T)
      p1 = p0 + int(stride, C_INT64_T)
      q = int(send_offset + 2_C_INT*(local_line - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = ys_gpsv_x(p0)
      ycomm_sendbuf(q + 2_C_INT64_T) = ys_gpsv_x(p1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine ys_pack_backward_state
#endif

  subroutine ys_unpack_pipelined_solution(dst, active_n, nlines, nlines_z, nx_count, dst_row_base)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)
    integer(C_INT), intent(in) :: active_n, nlines, nlines_z, nx_count, dst_row_base
    integer(C_INT) :: ix, iz_index, irow, iline
    integer(C_INT64_T) :: p
    logical :: has_padded_dst, has_left_interface, has_right_interface

    has_padded_dst = (dst_row_base == 3_C_INT)
    has_left_interface = (ipy > 0_C_INT)
    has_right_interface = (ipy < npy_grid - 1_C_INT)
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(dst, ys_gpsv_x, active_n, nlines, nlines_z, nx_count, dst_row_base) private(ix, iz_index, irow, iline, p)
    do ix = 1_C_INT, nx_count
      do iz_index = 1_C_INT, nlines_z
        do irow = 1_C_INT, active_n
          iline = (ix - 1_C_INT)*nlines_z + iz_index
          p = int(irow - 1_C_INT, C_INT64_T)*int(nlines, C_INT64_T) + int(iline, C_INT64_T)
          dst(irow - 1_C_INT + dst_row_base, iz_index, ix) = ys_gpsv_x(p)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do

    if (has_padded_dst .and. npy_grid > 1_C_INT) then
      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp shared(dst, ys_left_interface_values, ys_right_interface_values, active_n, nlines_z, nx_count, has_left_interface, has_right_interface) &
      !$omp private(ix, iz_index, iline)
      do ix = 1_C_INT, nx_count
        do iz_index = 1_C_INT, nlines_z
          iline = (ix - 1_C_INT)*nlines_z + iz_index
          if (has_left_interface) then
            dst(1, iz_index, ix) = ys_left_interface_values(1, iline)
            dst(2, iz_index, ix) = ys_left_interface_values(2, iline)
          end if
          if (has_right_interface) then
            dst(active_n + 3_C_INT, iz_index, ix) = ys_right_interface_values(1, iline)
            dst(active_n + 4_C_INT, iz_index, ix) = ys_right_interface_values(2, iline)
          end if
        end do
      end do
      !$omp end target teams distribute parallel do
    end if
  end subroutine ys_unpack_pipelined_solution

  subroutine ys_solve_reduced_interfaces()
    implicit none

    call ys_pack_reduced_leaf_rows()
    call ys_schur_solve_from_packed(ys_reduced_schur_ws)
    call ys_unpack_reduced_leaf_values()
  end subroutine ys_solve_reduced_interfaces

  subroutine ys_pack_reduced_leaf_rows()
    implicit none
#ifdef HAVE_MPI
    integer(C_INT) :: arity, nlines, base_count, remainder
    integer(C_INT) :: dest, local_line, irow, first_line, line_count, global_line, offset

    arity = ys_reduced_pass_counts(1)
    nlines = size(ys_reduced_rows_send, 2, kind=C_INT)
    base_count = nlines/arity
    remainder = mod(nlines, arity)

    call roctxPush("ys_schur_pack_leaf_rows")
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(arity, base_count, remainder, ys_reduced_rows_send, ycomm_sendbuf) &
    !$omp private(dest, local_line, irow, first_line, line_count, global_line, offset)
    do dest = 0, arity - 1
      do local_line = 1, base_count + 1
        do irow = 1, 20_C_INT
          line_count = base_count
          if (dest < remainder) line_count = line_count + 1_C_INT
          if (local_line <= line_count) then
            first_line = dest*base_count + min(dest, remainder) + 1_C_INT
            global_line = first_line + local_line - 1_C_INT
            offset = 20_C_INT*(first_line - 1_C_INT) + (local_line - 1_C_INT)*20_C_INT + irow
            ycomm_sendbuf(offset) = ys_reduced_rows_send(irow, global_line)
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_schur_pack_leaf_rows")
#else
    error stop "ys_pack_reduced_leaf_rows requires MPI"
#endif
  end subroutine ys_pack_reduced_leaf_rows

  subroutine ys_unpack_reduced_leaf_values()
    implicit none
#ifdef HAVE_MPI
    integer(C_INT) :: arity, nlines, base_count, remainder, row0
    integer(C_INT) :: src, local_line, k, first_line, line_count, global_line, offset

    arity = ys_reduced_pass_counts(1)
    nlines = size(ys_reduced_rows_send, 2, kind=C_INT)
    base_count = nlines/arity
    remainder = mod(nlines, arity)
    row0 = 4_C_INT*ipy

    call roctxPush("ys_schur_unpack_leaf_values")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(arity, base_count, remainder, row0, ycomm_recvbuf, ys_reduced_rhs, &
    !$omp& ys_left_interface_values, ys_right_interface_values) &
    !$omp private(src, local_line, k, first_line, line_count, global_line, offset)
    do src = 0, arity - 1
      do local_line = 1, base_count + 1
        line_count = base_count
        if (src < remainder) line_count = line_count + 1_C_INT
        if (local_line <= line_count) then
          first_line = src*base_count + min(src, remainder) + 1_C_INT
          global_line = first_line + local_line - 1_C_INT
          offset = 8_C_INT*(first_line - 1_C_INT) + (local_line - 1_C_INT)*8_C_INT
          do k = 1, 4_C_INT
            ys_reduced_rhs(row0 + k, global_line) = ycomm_recvbuf(offset + k)
          end do
          ys_left_interface_values(1, global_line) = ycomm_recvbuf(offset + 5_C_INT)
          ys_left_interface_values(2, global_line) = ycomm_recvbuf(offset + 6_C_INT)
          ys_right_interface_values(1, global_line) = ycomm_recvbuf(offset + 7_C_INT)
          ys_right_interface_values(2, global_line) = ycomm_recvbuf(offset + 8_C_INT)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_schur_unpack_leaf_values")
#else
    error stop "ys_unpack_reduced_leaf_values requires MPI"
#endif
  end subroutine ys_unpack_reduced_leaf_values

  ! Pentadiagonal LU factorization of one interleaved line.
  !
  ! open_right = .true. leaves the trailing row coupled to the rows owned by the
  ! next y rank, so its upper diagonal is carried forward as well; that is the
  ! only difference between the local and the pipelined-continuation variants.

  subroutine ys_factor_penta_interleaved(ds, dl, d, du, dw, stride, first, n, open_right)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    logical, intent(in) :: open_right
    integer(C_INT) :: i
    integer(C_INT64_T) :: p, p1, p2
    complex(C_DOUBLE_COMPLEX) :: factor

    if (n <= 0) return
    if (n == 1) then
      d(first) = 1.0d0/d(first)
      return
    end if

    do i = 0, n - 3
      p = int(first, C_INT64_T) + int(i, C_INT64_T)*int(stride, C_INT64_T)
      d(p) = 1.0d0/d(p)

      p1 = p + int(stride, C_INT64_T)
      factor = dl(p1)*d(p)
      dl(p1) = factor
      d(p1) = d(p1) - factor*du(p)
      du(p1) = du(p1) - factor*dw(p)

      p2 = p + 2_C_INT64_T*int(stride, C_INT64_T)
      factor = ds(p2)*d(p)
      ds(p2) = factor
      dl(p2) = dl(p2) - factor*du(p)
      d(p2) = d(p2) - factor*dw(p)
    end do

    p = int(first, C_INT64_T) + int(n - 2, C_INT64_T)*int(stride, C_INT64_T)
    d(p) = 1.0d0/d(p)
    p1 = p + int(stride, C_INT64_T)
    factor = dl(p1)*d(p)
    dl(p1) = factor
    d(p1) = d(p1) - factor*du(p)
    if (open_right) du(p1) = du(p1) - factor*dw(p)

    d(p1) = 1.0d0/d(p1)
  end subroutine ys_factor_penta_interleaved

  subroutine ys_factor_penta_interleaved_continue(ds, dl, d, du, stride, first, n, &
                                                  prev0_d, prev0_du, prev0_dw, prev1_d, prev1_du, prev1_dw)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:)
    integer(C_INT), intent(in) :: stride, first, n
    complex(C_DOUBLE_COMPLEX), intent(in) :: prev0_d, prev0_du, prev0_dw, prev1_d, prev1_du, prev1_dw
    integer(C_INT64_T) :: p0, p1
    complex(C_DOUBLE_COMPLEX) :: factor

    if (n <= 0) return

    p0 = int(first, C_INT64_T)
    factor = ds(p0)*prev0_d
    ds(p0) = factor
    dl(p0) = dl(p0) - factor*prev0_du
    d(p0) = d(p0) - factor*prev0_dw

    factor = dl(p0)*prev1_d
    dl(p0) = factor
    d(p0) = d(p0) - factor*prev1_du
    du(p0) = du(p0) - factor*prev1_dw

    if (n >= 2_C_INT) then
      p1 = p0 + int(stride, C_INT64_T)
      factor = ds(p1)*prev1_d
      ds(p1) = factor
      dl(p1) = dl(p1) - factor*prev1_du
      d(p1) = d(p1) - factor*prev1_dw
    end if
  end subroutine ys_factor_penta_interleaved_continue

  subroutine ys_forward_substitute_penta_interleaved(rhs, ds, dl, stride, first, n)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: ds(:), dl(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i
    integer(C_INT64_T) :: p

    if (n <= 0) return

    if (n >= 2) then
      p = int(first, C_INT64_T) + int(stride, C_INT64_T)
      rhs(p) = rhs(p) - dl(p)*rhs(first)
    end if
    do i = 2, n - 1
      p = int(first, C_INT64_T) + int(i, C_INT64_T)*int(stride, C_INT64_T)
      rhs(p) = rhs(p) - ds(p)*rhs(int(first, C_INT64_T) + int(i - 2, C_INT64_T)*int(stride, C_INT64_T)) - &
               dl(p)*rhs(int(first, C_INT64_T) + int(i - 1, C_INT64_T)*int(stride, C_INT64_T))
    end do
  end subroutine ys_forward_substitute_penta_interleaved

  subroutine ys_forward_substitute_penta_interleaved_continue(rhs, ds, dl, stride, first, n, prev0_rhs, prev1_rhs)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: ds(:), dl(:)
    integer(C_INT), intent(in) :: stride, first, n
    complex(C_DOUBLE_COMPLEX), intent(in) :: prev0_rhs, prev1_rhs
    integer(C_INT) :: i
    integer(C_INT64_T) :: p0, p1, p

    if (n <= 0) return

    p0 = int(first, C_INT64_T)
    rhs(p0) = rhs(p0) - ds(p0)*prev0_rhs - dl(p0)*prev1_rhs
    if (n >= 2_C_INT) then
      p1 = p0 + int(stride, C_INT64_T)
      rhs(p1) = rhs(p1) - ds(p1)*prev1_rhs - dl(p1)*rhs(p0)
    end if
    do i = 2, n - 1
      p = int(first, C_INT64_T) + int(i, C_INT64_T)*int(stride, C_INT64_T)
      rhs(p) = rhs(p) - ds(p)*rhs(int(first, C_INT64_T) + int(i - 2, C_INT64_T)*int(stride, C_INT64_T)) - &
               dl(p)*rhs(int(first, C_INT64_T) + int(i - 1, C_INT64_T)*int(stride, C_INT64_T))
    end do
  end subroutine ys_forward_substitute_penta_interleaved_continue

  subroutine ys_backward_substitute_penta_interleaved(rhs, d, du, dw, stride, first, n)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i
    integer(C_INT64_T) :: p

    if (n <= 0) return

    p = int(first, C_INT64_T) + int(n - 1, C_INT64_T)*int(stride, C_INT64_T)
    rhs(p) = rhs(p)*d(p)
    if (n >= 2) then
      p = int(first, C_INT64_T) + int(n - 2, C_INT64_T)*int(stride, C_INT64_T)
      rhs(p) = (rhs(p) - du(p)*rhs(p + int(stride, C_INT64_T)))*d(p)
    end if
    do i = n - 3, 0, -1
      p = int(first, C_INT64_T) + int(i, C_INT64_T)*int(stride, C_INT64_T)
      rhs(p) = (rhs(p) - du(p)*rhs(p + int(stride, C_INT64_T)) - &
                dw(p)*rhs(p + 2_C_INT64_T*int(stride, C_INT64_T)))*d(p)
    end do
  end subroutine ys_backward_substitute_penta_interleaved

  subroutine ys_backward_substitute_penta_interleaved_continue(rhs, d, du, dw, stride, first, n, next0_rhs, next1_rhs)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    complex(C_DOUBLE_COMPLEX), intent(in) :: next0_rhs, next1_rhs
    integer(C_INT) :: i
    integer(C_INT64_T) :: p, p_last, p_prev

    if (n <= 0) return

    p_last = int(first, C_INT64_T) + int(n - 1, C_INT64_T)*int(stride, C_INT64_T)
    rhs(p_last) = (rhs(p_last) - du(p_last)*next0_rhs - dw(p_last)*next1_rhs)*d(p_last)
    if (n >= 2_C_INT) then
      p_prev = p_last - int(stride, C_INT64_T)
      rhs(p_prev) = (rhs(p_prev) - du(p_prev)*rhs(p_last) - dw(p_prev)*next0_rhs)*d(p_prev)
    end if
    do i = n - 3, 0, -1
      p = int(first, C_INT64_T) + int(i, C_INT64_T)*int(stride, C_INT64_T)
      rhs(p) = (rhs(p) - du(p)*rhs(p + int(stride, C_INT64_T)) - &
                dw(p)*rhs(p + 2_C_INT64_T*int(stride, C_INT64_T)))*d(p)
    end do
  end subroutine ys_backward_substitute_penta_interleaved_continue

  subroutine ys_solve_factored_penta_interleaved(rhs, ds, dl, d, du, dw, stride, first, n)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n

    if (n <= 0) return

    call ys_forward_substitute_penta_interleaved(rhs, ds, dl, stride, first, n)
    call ys_backward_substitute_penta_interleaved(rhs, d, du, dw, stride, first, n)
  end subroutine ys_solve_factored_penta_interleaved

end module y_line_solvers
