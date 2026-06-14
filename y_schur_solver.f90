#include "header.h"

module y_schur_solver

  use, intrinsic :: iso_c_binding
  use roctx, only: roctxPush, roctxPop
#ifdef HAVE_MPI
  use mpi_transpose, only: MPI_COMM_Y, ensure_ycomm_buffers, ycomm_sendbuf, ycomm_recvbuf, ipy
  use mpi_f08
#endif

  implicit none
  private

  integer(C_INT), parameter :: YS_SCHUR_ROW_WIDTH = 20_C_INT
  integer(C_INT), parameter :: YS_SCHUR_VALUE_WIDTH = 8_C_INT
  integer(C_INT), parameter :: YS_SCHUR_BW = 5_C_INT
  integer(C_INT), parameter :: YS_SCHUR_MAX_ARITY = 16_C_INT
  integer(C_INT), parameter :: YS_SCHUR_MAX_ROWS = 4_C_INT*YS_SCHUR_MAX_ARITY
  integer(C_INT), parameter, public :: YS_SCHUR_EXCHANGE_AUTO = 0_C_INT
  integer(C_INT), parameter, public :: YS_SCHUR_EXCHANGE_ALLTOALL = 1_C_INT
  integer(C_INT), parameter, public :: YS_SCHUR_EXCHANGE_ALLGATHER = 2_C_INT
  integer, parameter :: YS_SCHUR_MAX_PASSES = 16

  type, public :: ys_schur_config
    integer(C_INT), allocatable :: pass_node_counts(:)
    logical :: comm_stats_enabled = .false.
    integer(C_INT) :: exchange_mode = YS_SCHUR_EXCHANGE_AUTO
  end type ys_schur_config

  type, public :: ys_schur_workspace
    integer(C_INT) :: npy_count = 1_C_INT
    integer(C_INT) :: nlines = 0_C_INT
    integer(C_INT) :: pass_count = 0_C_INT
    logical :: prepared = .false.
  end type ys_schur_workspace

  public :: ys_schur_configure
  public :: ys_schur_default_pass_counts
  public :: ys_schur_prepare
  public :: ys_schur_release
  public :: ys_schur_solve_from_packed

  integer(C_INT), allocatable, save :: s_pass_counts(:)
  integer(C_INT), allocatable, save :: s_level_arity(:)
  integer(C_INT), allocatable, save :: s_level_exchange_mode(:), s_level_child_id(:)
  integer(C_INT), allocatable, save :: s_level_prev_first(:), s_level_prev_count(:)
  integer(C_INT), allocatable, save :: s_level_owned_first(:), s_level_owned_count(:)
  integer(C_INT), allocatable, save :: s_level_line_first(:, :), s_level_line_count(:, :)
  integer, allocatable, save :: s_row_send_counts(:, :), s_row_send_displs(:, :)
  integer, allocatable, save :: s_row_recv_counts(:, :), s_row_recv_displs(:, :)
  integer, allocatable, save :: s_value_send_counts(:, :), s_value_send_displs(:, :)
  integer, allocatable, save :: s_value_recv_counts(:, :), s_value_recv_displs(:, :)
  integer(C_INT), allocatable, save :: s_row_send_elems(:), s_row_recv_elems(:)
  integer(C_INT), allocatable, save :: s_value_send_elems(:), s_value_recv_elems(:)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: s_rows(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: s_recover_basis(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: s_values(:, :, :)
  logical, save :: s_comm_stats_enabled = .false.

#ifdef HAVE_MPI
  type(MPI_Comm), allocatable, save :: s_level_comm(:)
  logical, allocatable, save :: s_level_comm_active(:)
#else
  integer(C_INT), parameter :: ipy = 0_C_INT
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ycomm_sendbuf(:), ycomm_recvbuf(:)
#endif

#if defined(HAVE_CUDA)
  !$omp declare target(ys_schur_factor_banded_complex_fixed)
  !$omp declare target(ys_schur_solve_factored_banded_complex_multi_fixed)
#endif

contains

  subroutine ys_schur_configure(cfg, pass_node_counts, comm_stats_enabled, exchange_mode)
    implicit none
    type(ys_schur_config), intent(out) :: cfg
    integer(C_INT), intent(in) :: pass_node_counts(:)
    logical, optional, intent(in) :: comm_stats_enabled
    integer(C_INT), optional, intent(in) :: exchange_mode

    allocate (cfg%pass_node_counts(size(pass_node_counts)))
    cfg%pass_node_counts = pass_node_counts
    cfg%comm_stats_enabled = .false.
    if (present(comm_stats_enabled)) cfg%comm_stats_enabled = comm_stats_enabled
    cfg%exchange_mode = YS_SCHUR_EXCHANGE_AUTO
    if (present(exchange_mode)) cfg%exchange_mode = exchange_mode
    select case (cfg%exchange_mode)
    case (YS_SCHUR_EXCHANGE_AUTO, YS_SCHUR_EXCHANGE_ALLTOALL, YS_SCHUR_EXCHANGE_ALLGATHER)
    case default
      error stop "unknown y-Schur exchange mode"
    end select
  end subroutine ys_schur_configure

  subroutine ys_schur_default_pass_counts(npy_count, pass_node_counts)
    implicit none
    integer(C_INT), intent(in) :: npy_count
    integer(C_INT), allocatable, intent(out) :: pass_node_counts(:)
    integer(C_INT) :: remaining, factor
    integer(C_INT) :: tmp(YS_SCHUR_MAX_PASSES)
    integer :: n

    if (npy_count < 1_C_INT) error stop "ys_schur_default_pass_counts requires npy_count >= 1"
    if (npy_count == 1_C_INT) then
      allocate (pass_node_counts(0))
      return
    end if

    remaining = npy_count
    n = 0
    do while (remaining > 1_C_INT)
      if (mod(remaining, 4_C_INT) == 0_C_INT) then
        factor = 4_C_INT
      else if (mod(remaining, 3_C_INT) == 0_C_INT) then
        factor = 3_C_INT
      else if (mod(remaining, 2_C_INT) == 0_C_INT) then
        factor = 2_C_INT
      else
        factor = remaining
      end if
      n = n + 1
      if (n > YS_SCHUR_MAX_PASSES) error stop "too many default y-Schur passes"
      tmp(n) = factor
      remaining = remaining/factor
    end do

    allocate (pass_node_counts(n))
    pass_node_counts = tmp(1:n)
  end subroutine ys_schur_default_pass_counts

  subroutine ys_schur_release(ws)
    implicit none
    type(ys_schur_workspace), intent(inout) :: ws
#ifdef HAVE_MPI
    integer :: ilevel, ierr_local
#endif

#ifdef HAVE_MPI
    if (allocated(s_level_comm)) then
      do ilevel = 1, size(s_level_comm)
        if (s_level_comm_active(ilevel)) then
          call MPI_Comm_free(s_level_comm(ilevel), ierr_local)
          if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_free y-Schur level communicator failed"
          s_level_comm(ilevel) = MPI_COMM_NULL
          s_level_comm_active(ilevel) = .false.
        end if
      end do
      deallocate (s_level_comm, s_level_comm_active)
    end if
#endif

    if (allocated(s_rows)) then
      !$omp target exit data map(delete: s_rows, s_recover_basis, s_values)
      deallocate (s_rows, s_recover_basis, s_values)
    end if

    if (allocated(s_pass_counts)) then
      !$omp target exit data map(delete: s_pass_counts, s_level_arity, s_level_exchange_mode, s_level_child_id, &
      !$omp& s_level_prev_first, s_level_prev_count, &
      !$omp& s_level_owned_first, s_level_owned_count, s_level_line_first, s_level_line_count, &
      !$omp& s_row_send_counts, s_row_send_displs, s_row_recv_counts, s_row_recv_displs, &
      !$omp& s_value_send_counts, s_value_send_displs, s_value_recv_counts, s_value_recv_displs, &
      !$omp& s_row_send_elems, s_row_recv_elems, s_value_send_elems, s_value_recv_elems)
      deallocate (s_pass_counts, s_level_arity, s_level_exchange_mode, s_level_child_id, &
                  s_level_prev_first, s_level_prev_count, &
                  s_level_owned_first, s_level_owned_count, s_level_line_first, s_level_line_count)
      deallocate (s_row_send_counts, s_row_send_displs, s_row_recv_counts, s_row_recv_displs)
      deallocate (s_value_send_counts, s_value_send_displs, s_value_recv_counts, s_value_recv_displs)
      deallocate (s_row_send_elems, s_row_recv_elems, s_value_send_elems, s_value_recv_elems)
    end if

    s_comm_stats_enabled = .false.
    ws%npy_count = 1_C_INT
    ws%nlines = 0_C_INT
    ws%pass_count = 0_C_INT
    ws%prepared = .false.
  end subroutine ys_schur_release

  subroutine ys_schur_prepare(ws, cfg, nlines, npy_count)
    implicit none
    type(ys_schur_workspace), intent(inout) :: ws
    type(ys_schur_config), intent(in) :: cfg
    integer(C_INT), intent(in) :: nlines, npy_count
    integer(C_INT) :: pass_product, prev_first, prev_count, owned_rel_first, owned_count
    integer(C_INT) :: span, child_span, rank_in_parent, child_id, child_pos, parent_group
    integer(C_INT) :: max_owned, max_arity, max_send_elems, max_recv_elems
    integer(C_INT) :: rel_first, line_count, offset_send, offset_recv
    integer :: ilevel, irank, ierr_local, comm_size

    call ys_schur_release(ws)
    if (npy_count < 1_C_INT) error stop "y-Schur requires at least one y rank"
    if (nlines < 1_C_INT) error stop "y-Schur requires at least one line"

    ws%npy_count = npy_count
    ws%nlines = nlines
    ws%pass_count = int(size(cfg%pass_node_counts), C_INT)
    ws%prepared = .true.
    s_comm_stats_enabled = cfg%comm_stats_enabled

    if (npy_count == 1_C_INT) then
      if (size(cfg%pass_node_counts) /= 0) error stop "single-rank y-Schur config must not contain passes"
      return
    end if
    if (size(cfg%pass_node_counts) < 1) error stop "distributed y-Schur config must contain at least one pass"
    if (size(cfg%pass_node_counts) > YS_SCHUR_MAX_PASSES) error stop "too many y-Schur passes"

    allocate (s_pass_counts(size(cfg%pass_node_counts)))
    s_pass_counts = cfg%pass_node_counts

    pass_product = 1_C_INT
    max_arity = 1_C_INT
    do ilevel = 1, size(s_pass_counts)
      if (s_pass_counts(ilevel) < 2_C_INT) error stop "distributed y-Schur pass arity must be at least two"
      if (s_pass_counts(ilevel) > YS_SCHUR_MAX_ARITY) error stop "y-Schur pass arity exceeds compiled maximum"
      if (mod(npy_count, pass_product*s_pass_counts(ilevel)) /= 0_C_INT) &
        error stop "y-Schur pass config must divide the y communicator exactly at every level"
      pass_product = pass_product*s_pass_counts(ilevel)
      max_arity = max(max_arity, s_pass_counts(ilevel))
    end do
    if (pass_product /= npy_count) error stop "y-Schur pass config product must equal the y communicator size"

    select case (cfg%exchange_mode)
    case (YS_SCHUR_EXCHANGE_AUTO, YS_SCHUR_EXCHANGE_ALLTOALL, YS_SCHUR_EXCHANGE_ALLGATHER)
    case default
      error stop "unknown y-Schur exchange mode"
    end select

    allocate (s_level_arity(ws%pass_count), s_level_exchange_mode(ws%pass_count), s_level_child_id(ws%pass_count))
    allocate (s_level_prev_first(ws%pass_count), s_level_prev_count(ws%pass_count))
    allocate (s_level_owned_first(ws%pass_count), s_level_owned_count(ws%pass_count))
    allocate (s_level_line_first(max_arity, ws%pass_count), s_level_line_count(max_arity, ws%pass_count))
    allocate (s_row_send_counts(max_arity, ws%pass_count), s_row_send_displs(max_arity, ws%pass_count))
    allocate (s_row_recv_counts(max_arity, ws%pass_count), s_row_recv_displs(max_arity, ws%pass_count))
    allocate (s_value_send_counts(max_arity, ws%pass_count), s_value_send_displs(max_arity, ws%pass_count))
    allocate (s_value_recv_counts(max_arity, ws%pass_count), s_value_recv_displs(max_arity, ws%pass_count))
    allocate (s_row_send_elems(ws%pass_count), s_row_recv_elems(ws%pass_count))
    allocate (s_value_send_elems(ws%pass_count), s_value_recv_elems(ws%pass_count))

    s_level_line_first = 1_C_INT
    s_level_line_count = 0_C_INT
    s_row_send_counts = 0
    s_row_send_displs = 0
    s_row_recv_counts = 0
    s_row_recv_displs = 0
    s_value_send_counts = 0
    s_value_send_displs = 0
    s_value_recv_counts = 0
    s_value_recv_displs = 0
    s_row_send_elems = 0_C_INT
    s_row_recv_elems = 0_C_INT
    s_value_send_elems = 0_C_INT
    s_value_recv_elems = 0_C_INT

#ifdef HAVE_MPI
    allocate (s_level_comm(ws%pass_count), s_level_comm_active(ws%pass_count))
    s_level_comm = MPI_COMM_NULL
    s_level_comm_active = .false.
#else
    error stop "distributed y-Schur requires MPI"
#endif

    max_owned = 1_C_INT
    max_send_elems = 0_C_INT
    max_recv_elems = 0_C_INT
    span = 1_C_INT
    prev_first = 1_C_INT
    prev_count = nlines

    do ilevel = 1, ws%pass_count
      s_level_arity(ilevel) = s_pass_counts(ilevel)
      child_span = span
      span = span*s_level_arity(ilevel)
      rank_in_parent = mod(ipy, span)
      child_id = rank_in_parent/child_span
      child_pos = mod(rank_in_parent, child_span)
      parent_group = ipy/span
      s_level_child_id(ilevel) = child_id
      select case (cfg%exchange_mode)
      case (YS_SCHUR_EXCHANGE_AUTO)
        if (ilevel == ws%pass_count .and. s_level_arity(ilevel) == 2_C_INT) then
          s_level_exchange_mode(ilevel) = YS_SCHUR_EXCHANGE_ALLGATHER
        else
          s_level_exchange_mode(ilevel) = YS_SCHUR_EXCHANGE_ALLTOALL
        end if
      case (YS_SCHUR_EXCHANGE_ALLTOALL)
        s_level_exchange_mode(ilevel) = cfg%exchange_mode
      case (YS_SCHUR_EXCHANGE_ALLGATHER)
        if (ilevel == ws%pass_count) then
          s_level_exchange_mode(ilevel) = YS_SCHUR_EXCHANGE_ALLGATHER
        else
          s_level_exchange_mode(ilevel) = YS_SCHUR_EXCHANGE_ALLTOALL
        end if
      end select

      if (prev_count < s_level_arity(ilevel)) &
        error stop "dense y-Schur level requires at least one line per participating rank"
      call ys_schur_split_range(child_id, prev_count, s_level_arity(ilevel), owned_rel_first, owned_count)

      s_level_prev_first(ilevel) = prev_first
      s_level_prev_count(ilevel) = prev_count
      s_level_owned_first(ilevel) = prev_first + owned_rel_first - 1_C_INT
      s_level_owned_count(ilevel) = owned_count
      max_owned = max(max_owned, prev_count)
      max_owned = max(max_owned, owned_count)

      offset_send = 0_C_INT
      offset_recv = 0_C_INT
      do irank = 0, int(s_level_arity(ilevel)) - 1
        call ys_schur_split_range(int(irank, C_INT), prev_count, s_level_arity(ilevel), rel_first, line_count)
        s_level_line_first(irank + 1, ilevel) = prev_first + rel_first - 1_C_INT
        s_level_line_count(irank + 1, ilevel) = line_count
        s_row_send_counts(irank + 1, ilevel) = int(YS_SCHUR_ROW_WIDTH*line_count)
        if (s_level_exchange_mode(ilevel) == YS_SCHUR_EXCHANGE_ALLGATHER) then
          s_row_recv_counts(irank + 1, ilevel) = int(YS_SCHUR_ROW_WIDTH*prev_count)
        else
          s_row_recv_counts(irank + 1, ilevel) = int(YS_SCHUR_ROW_WIDTH*owned_count)
        end if
        s_value_send_counts(irank + 1, ilevel) = int(YS_SCHUR_VALUE_WIDTH*owned_count)
        s_row_send_displs(irank + 1, ilevel) = int(offset_send)
        s_row_recv_displs(irank + 1, ilevel) = int(offset_recv)
        s_value_send_displs(irank + 1, ilevel) = int(YS_SCHUR_VALUE_WIDTH*owned_count*irank)
        offset_send = offset_send + YS_SCHUR_ROW_WIDTH*line_count
        offset_recv = offset_recv + s_row_recv_counts(irank + 1, ilevel)
      end do
      s_row_send_elems(ilevel) = offset_send
      s_row_recv_elems(ilevel) = offset_recv

      offset_recv = 0_C_INT
      do irank = 0, int(s_level_arity(ilevel)) - 1
        line_count = s_level_line_count(irank + 1, ilevel)
        if (s_level_exchange_mode(ilevel) == YS_SCHUR_EXCHANGE_ALLGATHER) then
          s_value_recv_counts(irank + 1, ilevel) = int(YS_SCHUR_VALUE_WIDTH*s_level_arity(ilevel)*line_count)
        else
          s_value_recv_counts(irank + 1, ilevel) = int(YS_SCHUR_VALUE_WIDTH*line_count)
        end if
        s_value_recv_displs(irank + 1, ilevel) = int(offset_recv)
        offset_recv = offset_recv + s_value_recv_counts(irank + 1, ilevel)
      end do
      s_value_send_elems(ilevel) = YS_SCHUR_VALUE_WIDTH*s_level_arity(ilevel)*owned_count
      s_value_recv_elems(ilevel) = offset_recv
      max_send_elems = max(max_send_elems, s_row_send_elems(ilevel), s_value_send_elems(ilevel))
      max_recv_elems = max(max_recv_elems, s_row_recv_elems(ilevel), s_value_recv_elems(ilevel))

#ifdef HAVE_MPI
      call MPI_Comm_split(MPI_COMM_Y, int(parent_group*child_span + child_pos), int(child_id), &
                          s_level_comm(ilevel), ierr_local)
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_split y-Schur level communicator failed"
      s_level_comm_active(ilevel) = .true.
      call MPI_Comm_size(s_level_comm(ilevel), comm_size, ierr_local)
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_size y-Schur level communicator failed"
      if (comm_size /= int(s_level_arity(ilevel))) error stop "unexpected y-Schur dense communicator size"
#endif

      prev_first = s_level_owned_first(ilevel)
      prev_count = s_level_owned_count(ilevel)
    end do

#ifdef HAVE_MPI
    call ensure_ycomm_buffers(max(1_C_INT, max_send_elems), max(1_C_INT, max_recv_elems))
#endif

    allocate (s_rows(YS_SCHUR_ROW_WIDTH, max_owned, ws%pass_count))
    allocate (s_recover_basis(4*max_arity, 5, max_owned, ws%pass_count))
    allocate (s_values(YS_SCHUR_VALUE_WIDTH, max_owned, ws%pass_count))
    !$omp target enter data map(alloc: s_rows, s_recover_basis, s_values)
    !$omp target enter data map(to: s_pass_counts, s_level_arity, s_level_exchange_mode, s_level_child_id, &
    !$omp& s_level_prev_first, s_level_prev_count, &
    !$omp& s_level_owned_first, s_level_owned_count, s_level_line_first, s_level_line_count, &
    !$omp& s_row_send_counts, s_row_send_displs, s_row_recv_counts, s_row_recv_displs, &
    !$omp& s_value_send_counts, s_value_send_displs, s_value_recv_counts, s_value_recv_displs, &
    !$omp& s_row_send_elems, s_row_recv_elems, s_value_send_elems, s_value_recv_elems)
  end subroutine ys_schur_prepare

  subroutine ys_schur_solve_from_packed(ws)
    implicit none
    type(ys_schur_workspace), intent(in) :: ws
    integer :: ilevel
    logical :: redundant_root

    if (.not. ws%prepared) error stop "ys_schur_solve_from_packed called before ys_schur_prepare"
    if (ws%npy_count == 1_C_INT) return

    redundant_root = (s_level_exchange_mode(ws%pass_count) == YS_SCHUR_EXCHANGE_ALLGATHER)

    call ys_schur_exchange_rows(1)
    call ys_schur_compose_level(1, redundant_root .and. ws%pass_count == 1_C_INT)

    do ilevel = 2, int(ws%pass_count)
      call ys_schur_pack_parent_rows(ilevel)
      call ys_schur_exchange_rows(ilevel)
      call ys_schur_compose_level(ilevel, redundant_root .and. ilevel == int(ws%pass_count))
    end do

    call ys_schur_seed_root_values(ws%pass_count, redundant_root)
    do ilevel = int(ws%pass_count), 1, -1
      if (redundant_root .and. ilevel == int(ws%pass_count)) then
        call ys_schur_recover_redundant_root_values(ilevel)
      else
        call ys_schur_pack_recovered_values(ilevel)
        call ys_schur_exchange_values(ilevel)
        if (ilevel > 1) then
          call ys_schur_unpack_parent_values(ilevel)
        end if
      end if
    end do
  end subroutine ys_schur_solve_from_packed

  subroutine ys_schur_pack_parent_rows(ilevel)
    implicit none
    integer, intent(in) :: ilevel
    integer(C_INT) :: dest, local_line, irow, first_line, line_count, global_line, src_line, offset

    call roctxPush("ys_schur_pack_rows")
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(ilevel, ycomm_sendbuf, s_rows, s_level_arity, s_level_line_first, &
    !$omp& s_level_line_count, s_level_prev_first, s_row_send_displs) &
    !$omp private(dest, local_line, irow, first_line, line_count, global_line, src_line, offset)
    do dest = 0, s_level_arity(ilevel) - 1
      do local_line = 1, s_level_line_count(dest + 1, ilevel)
        do irow = 1, YS_SCHUR_ROW_WIDTH
          first_line = s_level_line_first(dest + 1, ilevel)
          line_count = s_level_line_count(dest + 1, ilevel)
          global_line = first_line + local_line - 1
          offset = s_row_send_displs(dest + 1, ilevel) + (local_line - 1)*YS_SCHUR_ROW_WIDTH + irow
          src_line = global_line - s_level_prev_first(ilevel) + 1
          ycomm_sendbuf(offset) = s_rows(irow, src_line, ilevel - 1)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_schur_pack_rows")
  end subroutine ys_schur_pack_parent_rows

  subroutine ys_schur_exchange_rows(ilevel)
    implicit none
    integer, intent(in) :: ilevel
#ifdef HAVE_MPI
    integer :: ierr_local
    real(C_DOUBLE) :: comm_t0, elapsed

    comm_t0 = 0.0_C_DOUBLE
    elapsed = 0.0_C_DOUBLE
    if (s_level_exchange_mode(ilevel) == YS_SCHUR_EXCHANGE_ALLGATHER) then
      call roctxPush("MPI_Allgather ys_schur_rows")
      if (s_comm_stats_enabled) comm_t0 = MPI_Wtime()
      !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
      call MPI_Allgather(ycomm_sendbuf(1:s_row_send_elems(ilevel)), s_row_send_elems(ilevel), MPI_DOUBLE_COMPLEX, &
                         ycomm_recvbuf(1:s_row_recv_elems(ilevel)), s_row_send_elems(ilevel), MPI_DOUBLE_COMPLEX, &
                         s_level_comm(ilevel), ierr_local)
      !$omp end target data
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Allgather y-Schur rows failed"
      if (s_comm_stats_enabled) then
        elapsed = MPI_Wtime() - comm_t0
        call ys_schur_report_comm_stats("ys_schur_rows_allgather", s_level_comm(ilevel), &
                                        s_row_send_elems(ilevel), s_row_recv_elems(ilevel), elapsed)
      end if
      call roctxPop("MPI_Allgather ys_schur_rows")
    else
      call roctxPush("MPI_Alltoallv ys_schur_rows")
      if (s_comm_stats_enabled) comm_t0 = MPI_Wtime()
      !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
      call MPI_Alltoallv(ycomm_sendbuf(1:s_row_send_elems(ilevel)), &
                         s_row_send_counts(1:s_level_arity(ilevel), ilevel), &
                         s_row_send_displs(1:s_level_arity(ilevel), ilevel), MPI_DOUBLE_COMPLEX, &
                         ycomm_recvbuf(1:s_row_recv_elems(ilevel)), &
                         s_row_recv_counts(1:s_level_arity(ilevel), ilevel), &
                         s_row_recv_displs(1:s_level_arity(ilevel), ilevel), MPI_DOUBLE_COMPLEX, &
                         s_level_comm(ilevel), ierr_local)
      !$omp end target data
      if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv y-Schur rows failed"
      if (s_comm_stats_enabled) then
        elapsed = MPI_Wtime() - comm_t0
        call ys_schur_report_comm_stats("ys_schur_rows_alltoall", s_level_comm(ilevel), &
                                        s_row_send_elems(ilevel), s_row_recv_elems(ilevel), elapsed)
      end if
      call roctxPop("MPI_Alltoallv ys_schur_rows")
    end if
#else
    error stop "ys_schur_exchange_rows requires MPI"
#endif
  end subroutine ys_schur_exchange_rows

  subroutine ys_schur_exchange_values(ilevel)
    implicit none
    integer, intent(in) :: ilevel
#ifdef HAVE_MPI
    integer :: ierr_local
    real(C_DOUBLE) :: comm_t0, elapsed

    comm_t0 = 0.0_C_DOUBLE
    elapsed = 0.0_C_DOUBLE
    if (s_level_exchange_mode(ilevel) == YS_SCHUR_EXCHANGE_ALLGATHER) then
      error stop "root allgather Schur levels should recover values locally"
    end if

    call roctxPush("MPI_Alltoallv ys_schur_values")
    if (s_comm_stats_enabled) comm_t0 = MPI_Wtime()
    !$omp target data use_device_addr(ycomm_sendbuf, ycomm_recvbuf)
    call MPI_Alltoallv(ycomm_sendbuf(1:s_value_send_elems(ilevel)), &
                       s_value_send_counts(1:s_level_arity(ilevel), ilevel), &
                       s_value_send_displs(1:s_level_arity(ilevel), ilevel), MPI_DOUBLE_COMPLEX, &
                       ycomm_recvbuf(1:s_value_recv_elems(ilevel)), &
                       s_value_recv_counts(1:s_level_arity(ilevel), ilevel), &
                       s_value_recv_displs(1:s_level_arity(ilevel), ilevel), MPI_DOUBLE_COMPLEX, &
                       s_level_comm(ilevel), ierr_local)
    !$omp end target data
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Alltoallv y-Schur values failed"
    if (s_comm_stats_enabled) then
      elapsed = MPI_Wtime() - comm_t0
      call ys_schur_report_comm_stats("ys_schur_values_alltoall", s_level_comm(ilevel), &
                                      s_value_send_elems(ilevel), s_value_recv_elems(ilevel), elapsed)
    end if
    call roctxPop("MPI_Alltoallv ys_schur_values")
#else
    error stop "ys_schur_exchange_values requires MPI"
#endif
  end subroutine ys_schur_exchange_values

  subroutine ys_schur_compose_level(ilevel, solve_redundant)
    implicit none
    integer, intent(in) :: ilevel
    logical, intent(in) :: solve_redundant
    integer(C_INT) :: iline, child, row0, row, k, ext_col, exposed_var, offset, col, rel_line
    integer(C_INT) :: n, i, j, irhs, solve_count
    complex(C_DOUBLE_COMPLEX) :: a(YS_SCHUR_MAX_ROWS, 2*YS_SCHUR_BW + 1), rhs(YS_SCHUR_MAX_ROWS, 5)

    solve_count = s_level_owned_count(ilevel)
    if (solve_redundant) solve_count = s_level_prev_count(ilevel)

    call roctxPush("ys_schur_compose_level")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ilevel, ycomm_recvbuf, s_row_recv_displs, s_recover_basis, s_rows, &
    !$omp& s_level_owned_count, s_level_owned_first, s_level_prev_first, s_level_prev_count, &
    !$omp& s_level_arity, s_level_exchange_mode, solve_redundant, solve_count) &
    !$omp private(iline, child, row0, row, k, ext_col, exposed_var, offset, col, &
    !$omp& rel_line, n, i, j, irhs, a, rhs)
    do iline = 1, solve_count
      n = 4*s_level_arity(ilevel)
      do j = 1, 2*YS_SCHUR_BW + 1
        do i = 1, YS_SCHUR_MAX_ROWS
          a(i, j) = (0.0d0, 0.0d0)
        end do
      end do
      do irhs = 1, 5
        do i = 1, YS_SCHUR_MAX_ROWS
          rhs(i, irhs) = (0.0d0, 0.0d0)
        end do
      end do

      do child = 0, s_level_arity(ilevel) - 1
        if (s_level_exchange_mode(ilevel) == YS_SCHUR_EXCHANGE_ALLGATHER) then
          if (solve_redundant) then
            rel_line = iline
          else
            rel_line = s_level_owned_first(ilevel) - s_level_prev_first(ilevel) + iline
          end if
          offset = s_row_recv_displs(child + 1, ilevel) + (rel_line - 1)*YS_SCHUR_ROW_WIDTH
        else
          offset = s_row_recv_displs(child + 1, ilevel) + (iline - 1)*YS_SCHUR_ROW_WIDTH
        end if
        row0 = 4*child
        do k = 1, 4
          row = row0 + k
          a(row, YS_SCHUR_BW + 1) = (1.0d0, 0.0d0)
          rhs(row, 1) = ycomm_recvbuf(offset + k)
          if (child > 0) then
            col = row0 - 1
            a(row, YS_SCHUR_BW + 1 + col - row) = -ycomm_recvbuf(offset + 4 + k)
            col = row0
            a(row, YS_SCHUR_BW + 1 + col - row) = -ycomm_recvbuf(offset + 8 + k)
          else
            rhs(row, 2) = ycomm_recvbuf(offset + 4 + k)
            rhs(row, 3) = ycomm_recvbuf(offset + 8 + k)
          end if
          if (child < s_level_arity(ilevel) - 1) then
            col = row0 + 5
            a(row, YS_SCHUR_BW + 1 + col - row) = -ycomm_recvbuf(offset + 12 + k)
            col = row0 + 6
            a(row, YS_SCHUR_BW + 1 + col - row) = -ycomm_recvbuf(offset + 16 + k)
          else
            rhs(row, 4) = ycomm_recvbuf(offset + 12 + k)
            rhs(row, 5) = ycomm_recvbuf(offset + 16 + k)
          end if
        end do
      end do

      call ys_schur_factor_banded_complex_fixed(a, n)
      call ys_schur_solve_factored_banded_complex_multi_fixed(rhs, a, n, 5_C_INT)

      do irhs = 1, 5
        do row = 1, n
          s_recover_basis(row, irhs, iline, ilevel) = rhs(row, irhs)
        end do
      end do

      do k = 1, 4
        if (k <= 2) then
          exposed_var = k
        else
          exposed_var = 4*(s_level_arity(ilevel) - 1) + k
        end if
        s_rows(k, iline, ilevel) = rhs(exposed_var, 1)
        do ext_col = 1, 4
          s_rows(4*ext_col + k, iline, ilevel) = rhs(exposed_var, ext_col + 1)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_schur_compose_level")
  end subroutine ys_schur_compose_level

  subroutine ys_schur_seed_root_values(root_level, solve_redundant)
    implicit none
    integer(C_INT), intent(in) :: root_level
    logical, intent(in) :: solve_redundant
    integer(C_INT) :: iline, k, solve_count

    solve_count = s_level_owned_count(root_level)
    if (solve_redundant) solve_count = s_level_prev_count(root_level)

    call roctxPush("ys_schur_seed_root_values")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(root_level, solve_count, s_values, s_rows) private(iline, k)
    do iline = 1, solve_count
      do k = 1, YS_SCHUR_VALUE_WIDTH
        if (k <= 4) then
          s_values(k, iline, root_level) = s_rows(k, iline, root_level)
        else
          s_values(k, iline, root_level) = (0.0d0, 0.0d0)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_schur_seed_root_values")
  end subroutine ys_schur_seed_root_values

  subroutine ys_schur_recover_redundant_root_values(ilevel)
    implicit none
    integer, intent(in) :: ilevel
    integer(C_INT) :: iline, k, row0, child_id, offset
    complex(C_DOUBLE_COMPLEX) :: prev1, prev2, next1, next2
    complex(C_DOUBLE_COMPLEX) :: recovered(YS_SCHUR_MAX_ROWS)

    call roctxPush("ys_schur_recover_redundant_root_values")
    if (ilevel > 1) then
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(ilevel, s_values, s_recover_basis, s_level_prev_count, &
      !$omp& s_level_arity, s_level_child_id) &
      !$omp private(iline, k, row0, child_id, prev1, prev2, next1, next2, recovered)
      do iline = 1, s_level_prev_count(ilevel)
        child_id = s_level_child_id(ilevel)
        row0 = 4_C_INT*child_id
        prev1 = s_values(5, iline, ilevel)
        prev2 = s_values(6, iline, ilevel)
        next1 = s_values(7, iline, ilevel)
        next2 = s_values(8, iline, ilevel)
        do k = 1, 4_C_INT*s_level_arity(ilevel)
          recovered(k) = s_recover_basis(k, 1, iline, ilevel) + &
                         s_recover_basis(k, 2, iline, ilevel)*prev1 + &
                         s_recover_basis(k, 3, iline, ilevel)*prev2 + &
                         s_recover_basis(k, 4, iline, ilevel)*next1 + &
                         s_recover_basis(k, 5, iline, ilevel)*next2
        end do

        s_values(1, iline, ilevel - 1) = recovered(row0 + 1)
        s_values(2, iline, ilevel - 1) = recovered(row0 + 2)
        s_values(3, iline, ilevel - 1) = recovered(row0 + 3)
        s_values(4, iline, ilevel - 1) = recovered(row0 + 4)
        if (child_id > 0_C_INT) then
          s_values(5, iline, ilevel - 1) = recovered(row0 - 1)
          s_values(6, iline, ilevel - 1) = recovered(row0)
        else
          s_values(5, iline, ilevel - 1) = prev1
          s_values(6, iline, ilevel - 1) = prev2
        end if
        if (child_id < s_level_arity(ilevel) - 1_C_INT) then
          s_values(7, iline, ilevel - 1) = recovered(row0 + 5)
          s_values(8, iline, ilevel - 1) = recovered(row0 + 6)
        else
          s_values(7, iline, ilevel - 1) = next1
          s_values(8, iline, ilevel - 1) = next2
        end if
      end do
      !$omp end target teams distribute parallel do
    else
      !$omp target teams distribute parallel do default(none) &
      !$omp shared(ilevel, ycomm_recvbuf, s_values, s_recover_basis, s_level_prev_count, &
      !$omp& s_level_arity, s_level_child_id) &
      !$omp private(iline, k, row0, child_id, offset, prev1, prev2, next1, next2, recovered)
      do iline = 1, s_level_prev_count(ilevel)
        child_id = s_level_child_id(ilevel)
        row0 = 4_C_INT*child_id
        prev1 = s_values(5, iline, ilevel)
        prev2 = s_values(6, iline, ilevel)
        next1 = s_values(7, iline, ilevel)
        next2 = s_values(8, iline, ilevel)
        do k = 1, 4_C_INT*s_level_arity(ilevel)
          recovered(k) = s_recover_basis(k, 1, iline, ilevel) + &
                         s_recover_basis(k, 2, iline, ilevel)*prev1 + &
                         s_recover_basis(k, 3, iline, ilevel)*prev2 + &
                         s_recover_basis(k, 4, iline, ilevel)*next1 + &
                         s_recover_basis(k, 5, iline, ilevel)*next2
        end do

        offset = (iline - 1_C_INT)*YS_SCHUR_VALUE_WIDTH
        ycomm_recvbuf(offset + 1) = recovered(row0 + 1)
        ycomm_recvbuf(offset + 2) = recovered(row0 + 2)
        ycomm_recvbuf(offset + 3) = recovered(row0 + 3)
        ycomm_recvbuf(offset + 4) = recovered(row0 + 4)
        if (child_id > 0_C_INT) then
          ycomm_recvbuf(offset + 5) = recovered(row0 - 1)
          ycomm_recvbuf(offset + 6) = recovered(row0)
        else
          ycomm_recvbuf(offset + 5) = prev1
          ycomm_recvbuf(offset + 6) = prev2
        end if
        if (child_id < s_level_arity(ilevel) - 1_C_INT) then
          ycomm_recvbuf(offset + 7) = recovered(row0 + 5)
          ycomm_recvbuf(offset + 8) = recovered(row0 + 6)
        else
          ycomm_recvbuf(offset + 7) = next1
          ycomm_recvbuf(offset + 8) = next2
        end if
      end do
      !$omp end target teams distribute parallel do
    end if
    call roctxPop("ys_schur_recover_redundant_root_values")
  end subroutine ys_schur_recover_redundant_root_values

  subroutine ys_schur_pack_recovered_values(ilevel)
    implicit none
    integer, intent(in) :: ilevel
    integer(C_INT) :: iline, child, row0, k, offset
    complex(C_DOUBLE_COMPLEX) :: prev1, prev2, next1, next2
    complex(C_DOUBLE_COMPLEX) :: recovered(YS_SCHUR_MAX_ROWS)

    call roctxPush("ys_schur_pack_recovered_values")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ilevel, ycomm_sendbuf, s_value_send_displs, s_level_owned_count, s_level_arity, &
    !$omp& s_recover_basis, s_values) &
    !$omp private(iline, child, row0, k, offset, prev1, prev2, next1, next2, recovered)
    do iline = 1, s_level_owned_count(ilevel)
      prev1 = s_values(5, iline, ilevel)
      prev2 = s_values(6, iline, ilevel)
      next1 = s_values(7, iline, ilevel)
      next2 = s_values(8, iline, ilevel)
      do k = 1, 4*s_level_arity(ilevel)
        recovered(k) = s_recover_basis(k, 1, iline, ilevel) + &
                       s_recover_basis(k, 2, iline, ilevel)*prev1 + &
                       s_recover_basis(k, 3, iline, ilevel)*prev2 + &
                       s_recover_basis(k, 4, iline, ilevel)*next1 + &
                       s_recover_basis(k, 5, iline, ilevel)*next2
      end do

      do child = 0, s_level_arity(ilevel) - 1
        row0 = 4*child
        offset = s_value_send_displs(child + 1, ilevel) + (iline - 1)*YS_SCHUR_VALUE_WIDTH
        ycomm_sendbuf(offset + 1) = recovered(row0 + 1)
        ycomm_sendbuf(offset + 2) = recovered(row0 + 2)
        ycomm_sendbuf(offset + 3) = recovered(row0 + 3)
        ycomm_sendbuf(offset + 4) = recovered(row0 + 4)
        if (child > 0) then
          ycomm_sendbuf(offset + 5) = recovered(row0 - 1)
          ycomm_sendbuf(offset + 6) = recovered(row0)
        else
          ycomm_sendbuf(offset + 5) = prev1
          ycomm_sendbuf(offset + 6) = prev2
        end if
        if (child < s_level_arity(ilevel) - 1) then
          ycomm_sendbuf(offset + 7) = recovered(row0 + 5)
          ycomm_sendbuf(offset + 8) = recovered(row0 + 6)
        else
          ycomm_sendbuf(offset + 7) = next1
          ycomm_sendbuf(offset + 8) = next2
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_schur_pack_recovered_values")
  end subroutine ys_schur_pack_recovered_values

  subroutine ys_schur_unpack_parent_values(ilevel)
    implicit none
    integer, intent(in) :: ilevel
    integer(C_INT) :: src, local_line, k, line_first, global_line, dst_line, offset

    call roctxPush("ys_schur_unpack_parent_values")
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(ilevel, ycomm_recvbuf, s_values, s_level_arity, s_level_line_first, s_level_line_count, &
    !$omp& s_value_recv_displs, s_level_prev_first) &
    !$omp private(src, local_line, k, line_first, global_line, dst_line, offset)
    do src = 0, s_level_arity(ilevel) - 1
      do local_line = 1, s_level_line_count(src + 1, ilevel)
        do k = 1, YS_SCHUR_VALUE_WIDTH
          line_first = s_level_line_first(src + 1, ilevel)
          global_line = line_first + local_line - 1
          dst_line = global_line - s_level_prev_first(ilevel) + 1
          offset = s_value_recv_displs(src + 1, ilevel) + (local_line - 1)*YS_SCHUR_VALUE_WIDTH
          s_values(k, dst_line, ilevel - 1) = ycomm_recvbuf(offset + k)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_schur_unpack_parent_values")
  end subroutine ys_schur_unpack_parent_values

  subroutine ys_schur_split_range(rank, nitems, nranks, first_item, item_count)
    implicit none
    integer(C_INT), intent(in) :: rank, nitems, nranks
    integer(C_INT), intent(out) :: first_item, item_count
    integer(C_INT) :: base_count, remainder

    base_count = nitems/nranks
    remainder = mod(nitems, nranks)
    item_count = base_count
    if (rank < remainder) item_count = item_count + 1_C_INT
    first_item = rank*base_count + min(rank, remainder) + 1_C_INT
  end subroutine ys_schur_split_range

#ifdef HAVE_MPI
  subroutine ys_schur_report_comm_stats(label, comm, send_elems, recv_elems, elapsed)
    implicit none
    character(*), intent(in) :: label
    type(MPI_Comm), intent(in) :: comm
    integer(C_INT), intent(in) :: send_elems, recv_elems
    real(C_DOUBLE), intent(in) :: elapsed
    integer :: comm_rank, comm_size, ierr_local
    real(C_DOUBLE) :: local_send_bytes, local_recv_bytes

    local_send_bytes = 16.0_C_DOUBLE*real(send_elems, C_DOUBLE)
    local_recv_bytes = 16.0_C_DOUBLE*real(recv_elems, C_DOUBLE)
    call MPI_Comm_rank(comm, comm_rank, ierr_local)
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_rank y-Schur comm stats failed"
    call MPI_Comm_size(comm, comm_size, ierr_local)
    if (ierr_local /= MPI_SUCCESS) error stop "MPI_Comm_size y-Schur comm stats failed"
    write (*, '("Y_SCHUR_COMM_STATS label=",a," ipy=",i0," comm_rank=",i0," comm_size=",i0,'// &
           '" send_MB=",f12.6," recv_MB=",f12.6," time_ms=",f12.6)') &
      trim(label), ipy, comm_rank, comm_size, local_send_bytes/1.0e6_C_DOUBLE, local_recv_bytes/1.0e6_C_DOUBLE, &
      elapsed*1.0e3_C_DOUBLE
  end subroutine ys_schur_report_comm_stats
#endif

  subroutine ys_schur_factor_banded_complex_fixed(a, n)
    implicit none
    integer(C_INT), intent(in) :: n
    complex(C_DOUBLE_COMPLEX), intent(inout) :: a(YS_SCHUR_MAX_ROWS, 2*YS_SCHUR_BW + 1)
    integer(C_INT) :: i, j, t, last
    complex(C_DOUBLE_COMPLEX) :: piv, factor

    do i = 1, n
      piv = a(i, YS_SCHUR_BW + 1)
      last = min(YS_SCHUR_BW, n - i)
      do j = 1, last
        factor = a(i + j, YS_SCHUR_BW + 1 - j)/piv
        a(i + j, YS_SCHUR_BW + 1 - j) = factor
        do t = 1, last
          a(i + j, YS_SCHUR_BW + 1 + t - j) = &
            a(i + j, YS_SCHUR_BW + 1 + t - j) - factor*a(i, YS_SCHUR_BW + 1 + t)
        end do
      end do
    end do
  end subroutine ys_schur_factor_banded_complex_fixed

  subroutine ys_schur_solve_factored_banded_complex_multi_fixed(rhs, a, n, nrhs)
    implicit none
    integer(C_INT), intent(in) :: n, nrhs
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(YS_SCHUR_MAX_ROWS, 5)
    complex(C_DOUBLE_COMPLEX), intent(in) :: a(YS_SCHUR_MAX_ROWS, 2*YS_SCHUR_BW + 1)
    integer(C_INT) :: i, j, irhs
    complex(C_DOUBLE_COMPLEX) :: factor, piv

    do i = 1, n
      do j = max(1_C_INT, i - YS_SCHUR_BW), i - 1
        factor = a(i, YS_SCHUR_BW + 1 + j - i)
        do irhs = 1, nrhs
          rhs(i, irhs) = rhs(i, irhs) - factor*rhs(j, irhs)
        end do
      end do
    end do

    do i = n, 1, -1
      do j = i + 1, min(n, i + YS_SCHUR_BW)
        factor = a(i, YS_SCHUR_BW + 1 + j - i)
        do irhs = 1, nrhs
          rhs(i, irhs) = rhs(i, irhs) - factor*rhs(j, irhs)
        end do
      end do
      piv = a(i, YS_SCHUR_BW + 1)
      do irhs = 1, nrhs
        rhs(i, irhs) = rhs(i, irhs)/piv
      end do
    end do
  end subroutine ys_schur_solve_factored_banded_complex_multi_fixed

end module y_schur_solver
