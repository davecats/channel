#!/usr/bin/env bash
set -euo pipefail

NP="${NP:-8}"
ITERS="${ITERS:-20}"
WARMUP="${WARMUP:-3}"
ACTIVE_N="${ACTIVE_N:-0}"
NLINES="${NLINES:-0}"
BATCHES="${BATCHES:-4}"
STORAGE_MODE="${STORAGE_MODE:-workspace}"
KERNEL_MODE="${KERNEL_MODE:-real}"
DETAIL="${DETAIL:-0}"
SWEEP_BATCHES="${SWEEP_BATCHES:-0}"
NX="${NX:-255}"
NY="${NY:-250}"
NZ="${NZ:-257}"
NPXZ="${NPXZ:-1}"
NPY="${NPY:-${NP}}"
OUT_DIR="${OUT_DIR:-omp_y_lu_pipeline_repro_$(date +%Y%m%d_%H%M%S)}"
PROFILE="${PROFILE:-0}"
TRACE="${TRACE:-nvtx,cuda}"
NVHPC_MODULE="${NVHPC_MODULE:-toolkits/nvhpc/25.5}"
MPIF90="${MPIF90:-mpifort}"
MPIRUN="${MPIRUN:-mpirun}"
MPIRUN_ARGS="${MPIRUN_ARGS:-}"
GPU_FLAGS="${GPU_FLAGS:--gpu=cc80}"
read -r -a MPIRUN_EXTRA_ARGS <<< "${MPIRUN_ARGS}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

if [[ -f /etc/profile.d/lmod.sh ]]; then
  # shellcheck disable=SC1091
  . /etc/profile.d/lmod.sh
  module load "${NVHPC_MODULE}" || true
fi

cat > bench_y_lu_pipeline_autotune.f90 <<'F90'
program bench_y_lu_pipeline_autotune
  use, intrinsic :: iso_c_binding
  use mpi_f08
  implicit none

  integer(C_SIZE_T), parameter :: WORKSPACE_ALIGNMENT = 256_C_SIZE_T
  integer(C_INT), parameter :: TAG_FACTOR = 8410_C_INT
  integer(C_INT), parameter :: TAG_FORWARD = 8510_C_INT
  integer(C_INT), parameter :: TAG_BACKWARD = 8610_C_INT
  integer(C_INT), parameter :: TAG_HALO_LOW = 8710_C_INT
  integer(C_INT), parameter :: TAG_HALO_HIGH = 8711_C_INT

  integer :: ierr, rank, nranks, arg_count
  integer(C_INT) :: active_n, nlines, nbatches, batch_max_lines
  integer(C_INT) :: nx_arg, ny_arg, nz_arg, npxz_arg, npy_arg, ipy
  integer :: iters, warmup, detail_arg
  character(64) :: arg
  character(64) :: kernel_mode
  logical :: detail_enabled, detail_active
  integer :: detail_count, detail_capacity
  character(len=32), allocatable :: detail_phase(:)
  integer(C_INT), allocatable :: detail_batch(:), detail_peer(:), detail_elems(:)
  real(C_DOUBLE), allocatable :: detail_t0(:), detail_t1(:)
  real(C_DOUBLE) :: time_offset

  integer(C_INT8_T), allocatable, target :: workspace_raw(:)
  integer(C_SIZE_T) :: workspace_base_index
  logical :: use_workspace_storage
  character(64) :: storage_mode
  complex(C_DOUBLE_COMPLEX), pointer, contiguous :: matrix_store(:), rhs_store(:)
  complex(C_DOUBLE_COMPLEX), pointer, contiguous :: ds(:), dl(:), d(:), du(:), dw(:), x(:)
  complex(C_DOUBLE_COMPLEX), allocatable, target :: ycomm_sendbuf(:), ycomm_recvbuf(:)
  complex(C_DOUBLE_COMPLEX), allocatable, target :: left_halo(:), right_halo(:)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

  iters = 20
  warmup = 3
  active_n = 64_C_INT
  nlines = 263680_C_INT
  nbatches = 4_C_INT
  nx_arg = 255_C_INT
  ny_arg = 250_C_INT
  nz_arg = 257_C_INT
  npxz_arg = 1_C_INT
  npy_arg = -1_C_INT
  storage_mode = "workspace"
  kernel_mode = "real"
  detail_arg = 0
  detail_enabled = .false.
  detail_active = .false.

  arg_count = command_argument_count()
  if (arg_count >= 1) then
    call get_command_argument(1, arg)
    read(arg, *) iters
  end if
  if (arg_count >= 2) then
    call get_command_argument(2, arg)
    read(arg, *) warmup
  end if
  if (arg_count >= 3) then
    call get_command_argument(3, arg)
    read(arg, *) active_n
  end if
  if (arg_count >= 4) then
    call get_command_argument(4, arg)
    read(arg, *) nlines
  end if
  if (arg_count >= 5) then
    call get_command_argument(5, arg)
    read(arg, *) nbatches
  end if
  if (arg_count >= 6) then
    call get_command_argument(6, storage_mode)
    storage_mode = adjustl(storage_mode)
  end if
  if (arg_count >= 7) then
    call get_command_argument(7, arg)
    read(arg, *) nx_arg
  end if
  if (arg_count >= 8) then
    call get_command_argument(8, arg)
    read(arg, *) ny_arg
  end if
  if (arg_count >= 9) then
    call get_command_argument(9, arg)
    read(arg, *) nz_arg
  end if
  if (arg_count >= 10) then
    call get_command_argument(10, arg)
    read(arg, *) npxz_arg
  end if
  if (arg_count >= 11) then
    call get_command_argument(11, arg)
    read(arg, *) npy_arg
  end if
  if (arg_count >= 12) then
    call get_command_argument(12, kernel_mode)
    kernel_mode = adjustl(kernel_mode)
  end if
  if (arg_count >= 13) then
    call get_command_argument(13, arg)
    read(arg, *) detail_arg
    detail_enabled = detail_arg /= 0
  end if
  if (npy_arg < 1_C_INT) npy_arg = nranks/npxz_arg

  if (npxz_arg < 1_C_INT) error stop "NPXZ must be >= 1"
  if (npy_arg < 1_C_INT) error stop "NPY must be >= 1"
  if (npxz_arg*npy_arg /= nranks) error stop "NP must equal NPXZ*NPY"
  if (mod(nx_arg + 1_C_INT, npxz_arg) /= 0_C_INT) error stop "NPXZ must divide NX+1"

  ipy = rank/npxz_arg
  if (active_n <= 0_C_INT) then
    active_n = ((ipy + 1_C_INT)*(ny_arg - 1_C_INT))/npy_arg - &
               (1_C_INT + ipy*(ny_arg - 1_C_INT)/npy_arg) + 1_C_INT
  end if
  if (nlines <= 0_C_INT) nlines = ((nx_arg + 1_C_INT)/npxz_arg)*(2_C_INT*nz_arg + 1_C_INT)

  if (active_n < 4_C_INT) error stop "ACTIVE_N must be >= 4"
  if (nlines < 1_C_INT) error stop "NLINES must be >= 1"
  if (nbatches < 1_C_INT) error stop "BATCHES must be >= 1"
  if (trim(kernel_mode) /= "real" .and. trim(kernel_mode) /= "dummy" .and. trim(kernel_mode) /= "none") &
    error stop "KERNEL_MODE must be real, dummy, or none"
  nbatches = min(nbatches, nlines)
  batch_max_lines = (nlines + nbatches - 1_C_INT)/nbatches
  use_workspace_storage = trim(storage_mode) /= "direct"
  detail_capacity = max(256, 80*nbatches + 128)
  allocate(detail_phase(detail_capacity), detail_batch(detail_capacity), detail_peer(detail_capacity), &
           detail_elems(detail_capacity), detail_t0(detail_capacity), detail_t1(detail_capacity))
  call calibrate_time_offset()

  call allocate_solver_like_storage()
  call initialize_system()

  if (rank == 0) then
    print '(a)', "OpenMP target y-line pipelined LU autotune probe"
    print '(a,i0,a,i0,a,i0)', "active_n=", active_n, " nlines=", nlines, " batches=", nbatches
    print '(a,i0,a,i0,a,i0,a,i0,a,i0)', "nx=", nx_arg, " ny=", ny_arg, " nz=", nz_arg, " npxz=", npxz_arg, " npy=", npy_arg
    print '(a,a)', "storage_mode=", trim(storage_mode)
    print '(a,i0,a,i0)', "iters=", iters, " warmup=", warmup
    print '(a,a,a,l1)', "kernel_mode=", trim(kernel_mode), " detail=", detail_enabled
    print '(a,i0,a,i0)', "factor_batch_bytes=", 6_C_INT*batch_max_lines*16_C_INT, &
      " solve_batch_bytes=", 2_C_INT*batch_max_lines*16_C_INT
    print '(a)', "rank, total_us, factor_wait_us, forward_wait_us, backward_wait_us, factor_kernel_us, forward_kernel_us, backward_kernel_us, pack_us, send_wait_us, halo_us"
  end if

  call run_loops(warmup, .false.)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call run_loops(iters, .true.)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (use_workspace_storage) then
    !$omp target exit data map(release: matrix_store, rhs_store)
    !$omp target exit data map(delete: workspace_raw)
    deallocate(workspace_raw)
  else
    !$omp target exit data map(delete: matrix_store, rhs_store)
    deallocate(matrix_store, rhs_store)
  end if
  !$omp target exit data map(delete: ycomm_sendbuf, ycomm_recvbuf, left_halo, right_halo)
  deallocate(ycomm_sendbuf, ycomm_recvbuf, left_halo, right_halo)
  deallocate(detail_phase, detail_batch, detail_peer, detail_elems, detail_t0, detail_t1)
  call MPI_Finalize(ierr)

contains
  subroutine allocate_solver_like_storage()
    integer(C_SIZE_T) :: nall, workspace_bytes, alloc_bytes, offset, count
    integer(C_INT) :: factor_stride, solve_stride, factor_total, solve_total

    nall = int(active_n, C_SIZE_T)*int(nlines, C_SIZE_T)
    if (use_workspace_storage) then
      call solver_workspace_bytes(active_n, nlines, nranks, workspace_bytes)
      alloc_bytes = workspace_bytes + WORKSPACE_ALIGNMENT
      allocate(workspace_raw(alloc_bytes))
      call choose_base_index(workspace_raw, workspace_base_index)
      !$omp target enter data map(alloc: workspace_raw)

      offset = 0_C_SIZE_T
      call bind_complex(offset, 5_C_SIZE_T*nall, matrix_store)
      call bind_complex(offset, nall, rhs_store)
      count = int(nlines, C_SIZE_T)
      call skip_complex(offset, count)
      call skip_complex(offset, count)
      call skip_complex(offset, count)
      call skip_complex(offset, count)
      call skip_real(offset, 5_C_SIZE_T*count)
      call skip_real(offset, 5_C_SIZE_T*count)
      call skip_real(offset, 5_C_SIZE_T*count)
      call skip_real(offset, 5_C_SIZE_T*count)
      call skip_complex(offset, count)
      call skip_complex(offset, count)
      call skip_real(offset, 4_C_SIZE_T*count)
      call skip_real(offset, 4_C_SIZE_T*count)
      call skip_complex(offset, nall)
      call skip_complex(offset, nall)
      call skip_complex(offset, nall)
      call skip_complex(offset, nall)
      call skip_complex(offset, nall)
      call skip_complex(offset, nall)
      if (nranks > 1) then
        call skip_complex(offset, 20_C_SIZE_T*count)
        call skip_complex(offset, 2_C_SIZE_T*count)
        call skip_complex(offset, 2_C_SIZE_T*count)
        call skip_complex(offset, 4_C_SIZE_T*int(nranks, C_SIZE_T)*count)
      end if
    else
      allocate(matrix_store(int(5_C_SIZE_T*nall)))
      allocate(rhs_store(int(nall)))
    end if

    ds(1:int(nall)) => matrix_store(1:int(nall))
    dl(1:int(nall)) => matrix_store(int(nall) + 1:2*int(nall))
    d(1:int(nall)) => matrix_store(2*int(nall) + 1:3*int(nall))
    du(1:int(nall)) => matrix_store(3*int(nall) + 1:4*int(nall))
    dw(1:int(nall)) => matrix_store(4*int(nall) + 1:5*int(nall))
    x(1:int(nall)) => rhs_store(1:int(nall))

    if (use_workspace_storage) then
      !$omp target enter data map(to: matrix_store, rhs_store)
    else
      !$omp target enter data map(alloc: matrix_store, rhs_store)
    end if

    factor_stride = 6_C_INT*batch_max_lines
    solve_stride = 2_C_INT*batch_max_lines
    factor_total = factor_stride*nbatches
    solve_total = solve_stride*nbatches
    allocate(ycomm_sendbuf(max(1_C_INT, factor_total + solve_total)))
    allocate(ycomm_recvbuf(max(1_C_INT, factor_total + 2_C_INT*solve_total)))
    allocate(left_halo(max(1_C_INT, 2_C_INT*nlines)))
    allocate(right_halo(max(1_C_INT, 2_C_INT*nlines)))
    !$omp target enter data map(alloc: ycomm_sendbuf, ycomm_recvbuf, left_halo, right_halo)
  end subroutine allocate_solver_like_storage

  subroutine bind_complex(offset_bytes, count_complex, target)
    integer(C_SIZE_T), intent(inout) :: offset_bytes
    integer(C_SIZE_T), intent(in) :: count_complex
    complex(C_DOUBLE_COMPLEX), pointer, contiguous, intent(out) :: target(:)
    complex(C_DOUBLE_COMPLEX), pointer :: tmp(:)
    type(C_PTR) :: ptr

    ptr = c_loc(workspace_raw(workspace_base_index + offset_bytes))
    call c_f_pointer(ptr, tmp, [int(count_complex)])
    target(1:int(count_complex)) => tmp
    offset_bytes = align_offset(offset_bytes + count_complex*int(C_SIZEOF((0.0_C_DOUBLE, 0.0_C_DOUBLE)), C_SIZE_T))
  end subroutine bind_complex

  subroutine initialize_system()
    integer(C_INT64_T) :: i
    integer(C_INT64_T) :: nall
    nall = int(active_n, C_INT64_T)*int(nlines, C_INT64_T)

    !$omp target teams distribute parallel do default(none) shared(ds, dl, d, du, dw, x, nall) private(i)
    do i = 1_C_INT64_T, nall
      ds(i) = (-0.010d0, 0.001d0)
      dl(i) = (-0.050d0, 0.002d0)
      d(i) = (1.250d0, 0.010d0)
      du(i) = (-0.040d0, -0.001d0)
      dw(i) = (-0.008d0, 0.001d0)
      x(i) = cmplx(0.125d0 + 1.0d-9*real(i, C_DOUBLE), -0.25d0 + 2.0d-9*real(i, C_DOUBLE), C_DOUBLE)
    end do
    !$omp end target teams distribute parallel do
  end subroutine initialize_system

  subroutine calibrate_time_offset()
    integer :: r
    integer(C_INT) :: token
    real(C_DOUBLE) :: t0, t1, remote_t, offset_r
    type(MPI_Status) :: status

    time_offset = 0.0_C_DOUBLE
    token = 0_C_INT
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    do r = 1, nranks - 1
      if (rank == 0) then
        t0 = MPI_Wtime()
        call MPI_Send(token, 1, MPI_INTEGER, r, 9101, MPI_COMM_WORLD, ierr)
        call MPI_Recv(remote_t, 1, MPI_DOUBLE_PRECISION, r, 9102, MPI_COMM_WORLD, status, ierr)
        t1 = MPI_Wtime()
        offset_r = 0.5d0*(t0 + t1) - remote_t
        call MPI_Send(offset_r, 1, MPI_DOUBLE_PRECISION, r, 9103, MPI_COMM_WORLD, ierr)
      else if (rank == r) then
        call MPI_Recv(token, 1, MPI_INTEGER, 0, 9101, MPI_COMM_WORLD, status, ierr)
        remote_t = MPI_Wtime()
        call MPI_Send(remote_t, 1, MPI_DOUBLE_PRECISION, 0, 9102, MPI_COMM_WORLD, ierr)
        call MPI_Recv(time_offset, 1, MPI_DOUBLE_PRECISION, 0, 9103, MPI_COMM_WORLD, status, ierr)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end do
  end subroutine calibrate_time_offset

  subroutine run_loops(nloop, report)
    integer, intent(in) :: nloop
    logical, intent(in) :: report
    integer :: loop, print_rank
    real(C_DOUBLE) :: total_us, factor_wait_us, forward_wait_us, backward_wait_us
    real(C_DOUBLE) :: factor_kernel_us, forward_kernel_us, backward_kernel_us
    real(C_DOUBLE) :: pack_us, send_wait_us, halo_us
    real(C_DOUBLE) :: accum(10), local(10)

    accum = 0.0_C_DOUBLE
    do loop = 1, nloop
      call initialize_system()
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      detail_active = detail_enabled .and. report .and. loop == 1
      if (detail_active) detail_count = 0
      call solve_once(local)
      accum = accum + local
      if (detail_active) call print_detail_events()
      detail_active = .false.
    end do

    if (.not. report) return
    if (nloop > 0) accum = accum/real(nloop, C_DOUBLE)

    total_us = accum(1)*1.0d6
    factor_wait_us = accum(2)*1.0d6
    forward_wait_us = accum(3)*1.0d6
    backward_wait_us = accum(4)*1.0d6
    factor_kernel_us = accum(5)*1.0d6
    forward_kernel_us = accum(6)*1.0d6
    backward_kernel_us = accum(7)*1.0d6
    pack_us = accum(8)*1.0d6
    send_wait_us = accum(9)*1.0d6
    halo_us = accum(10)*1.0d6

    do print_rank = 0, nranks - 1
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (rank == print_rank) then
        write(*,'(i0,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3)') &
          rank, total_us, factor_wait_us, forward_wait_us, backward_wait_us, factor_kernel_us, &
          forward_kernel_us, backward_kernel_us, pack_us, send_wait_us, halo_us
      end if
    end do
  end subroutine run_loops

  subroutine record_event(phase, batch, peer, elem_count, t_start, t_end)
    character(*), intent(in) :: phase
    integer(C_INT), intent(in) :: batch, peer, elem_count
    real(C_DOUBLE), intent(in) :: t_start, t_end
    if (.not. detail_active) return
    if (detail_count >= detail_capacity) return
    detail_count = detail_count + 1
    detail_phase(detail_count) = phase
    detail_batch(detail_count) = batch
    detail_peer(detail_count) = peer
    detail_elems(detail_count) = elem_count
    detail_t0(detail_count) = t_start + time_offset
    detail_t1(detail_count) = t_end + time_offset
  end subroutine record_event

  subroutine print_detail_events()
    integer :: print_rank, i
    do print_rank = 0, nranks - 1
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (rank == print_rank) then
        if (rank == 0) then
          write(*,'(a)') "detail, rank, phase, batch, peer, elems, start_s, end_s, dt_us"
        end if
        do i = 1, detail_count
          write(*,'(a,i0,", ",a,", ",i0,", ",i0,", ",i0,", ",f18.9,", ",f18.9,", ",f12.3)') &
            "detail, ", rank, trim(detail_phase(i)), detail_batch(i), detail_peer(i), detail_elems(i), &
            detail_t0(i), detail_t1(i), (detail_t1(i) - detail_t0(i))*1.0d6
        end do
      end if
    end do
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end subroutine print_detail_events

  subroutine solve_once(timing)
    real(C_DOUBLE), intent(out) :: timing(10)
    integer(C_INT) :: factor_stride, solve_stride, factor_total, solve_total
    integer(C_INT) :: batch, first_line, line_count, factor_count, solve_count
    integer(C_INT) :: factor_offset, forward_offset, backward_offset
    integer(C_INT) :: next_backward
    type(MPI_Request), allocatable :: factor_send_req(:), forward_send_req(:), backward_recv_req(:), backward_send_req(:)
    type(MPI_Request) :: factor_recv_req, forward_recv_req
    type(MPI_Status) :: status
    real(C_DOUBLE) :: t0, t1, t2

    timing = 0.0_C_DOUBLE
    t0 = MPI_Wtime()
    factor_stride = 6_C_INT*batch_max_lines
    solve_stride = 2_C_INT*batch_max_lines
    factor_total = factor_stride*nbatches
    solve_total = solve_stride*nbatches
    allocate(factor_send_req(nbatches), forward_send_req(nbatches), backward_recv_req(nbatches), backward_send_req(nbatches))
    factor_send_req = MPI_REQUEST_NULL
    forward_send_req = MPI_REQUEST_NULL
    backward_recv_req = MPI_REQUEST_NULL
    backward_send_req = MPI_REQUEST_NULL

    if (rank < nranks - 1) then
      do batch = 1_C_INT, nbatches
        call batch_range(batch, first_line, line_count)
        solve_count = 2_C_INT*line_count
        backward_offset = factor_total + solve_total + (batch - 1_C_INT)*solve_stride
        !$omp target data use_device_addr(ycomm_recvbuf)
        t1 = MPI_Wtime()
        call MPI_Irecv(ycomm_recvbuf(backward_offset + 1), solve_count, MPI_DOUBLE_COMPLEX, &
          rank + 1, TAG_BACKWARD + batch, MPI_COMM_WORLD, backward_recv_req(batch), ierr)
        t2 = MPI_Wtime()
        !$omp end target data
        call record_event("backward_irecv_post", batch, rank + 1, solve_count, t1, t2)
      end do
    end if

    next_backward = 1_C_INT
    do batch = 1_C_INT, nbatches
      call batch_range(batch, first_line, line_count)
      factor_count = 6_C_INT*line_count
      solve_count = 2_C_INT*line_count
      factor_offset = (batch - 1_C_INT)*factor_stride
      forward_offset = factor_total + (batch - 1_C_INT)*solve_stride
      if (rank > 0) then
        !$omp target data use_device_addr(ycomm_recvbuf)
        t1 = MPI_Wtime()
        call MPI_Irecv(ycomm_recvbuf(factor_offset + 1), factor_count, MPI_DOUBLE_COMPLEX, &
          rank - 1, TAG_FACTOR + batch, MPI_COMM_WORLD, factor_recv_req, ierr)
        t2 = MPI_Wtime()
        call record_event("factor_irecv_post", batch, rank - 1, factor_count, t1, t2)
        t1 = MPI_Wtime()
        call MPI_Irecv(ycomm_recvbuf(forward_offset + 1), solve_count, MPI_DOUBLE_COMPLEX, &
          rank - 1, TAG_FORWARD + batch, MPI_COMM_WORLD, forward_recv_req, ierr)
        t2 = MPI_Wtime()
        call record_event("forward_irecv_post", batch, rank - 1, solve_count, t1, t2)
        !$omp end target data
        t1 = MPI_Wtime()
        call MPI_Wait(factor_recv_req, status, ierr)
        t2 = MPI_Wtime()
        timing(2) = timing(2) + (t2 - t1)
        call record_event("factor_recv_wait", batch, rank - 1, factor_count, t1, t2)
        t1 = MPI_Wtime()
        call apply_factor_continuation(first_line, line_count, factor_offset)
        t2 = MPI_Wtime()
        timing(5) = timing(5) + (t2 - t1)
        call record_event("factor_continue_kernel", batch, rank, line_count, t1, t2)
      end if

      t1 = MPI_Wtime()
      call factor_batch(first_line, line_count, rank < nranks - 1)
      t2 = MPI_Wtime()
      timing(5) = timing(5) + (t2 - t1)
      call record_event("factor_kernel", batch, rank, line_count, t1, t2)

      if (rank < nranks - 1) then
        t1 = MPI_Wtime()
        call pack_factor_state(first_line, line_count, factor_offset)
        t2 = MPI_Wtime()
        timing(8) = timing(8) + (t2 - t1)
        call record_event("pack_factor", batch, rank + 1, factor_count, t1, t2)
        !$omp target data use_device_addr(ycomm_sendbuf)
        t1 = MPI_Wtime()
        call MPI_Isend(ycomm_sendbuf(factor_offset + 1), factor_count, MPI_DOUBLE_COMPLEX, &
          rank + 1, TAG_FACTOR + batch, MPI_COMM_WORLD, factor_send_req(batch), ierr)
        t2 = MPI_Wtime()
        !$omp end target data
        call record_event("factor_isend_post", batch, rank + 1, factor_count, t1, t2)
      end if

      if (rank > 0) then
        t1 = MPI_Wtime()
        call MPI_Wait(forward_recv_req, status, ierr)
        t2 = MPI_Wtime()
        timing(3) = timing(3) + (t2 - t1)
        call record_event("forward_recv_wait", batch, rank - 1, solve_count, t1, t2)
        t1 = MPI_Wtime()
        call forward_batch_continue(first_line, line_count, forward_offset)
      else
        t1 = MPI_Wtime()
        call forward_batch(first_line, line_count)
      end if
      t2 = MPI_Wtime()
      timing(6) = timing(6) + (t2 - t1)
      call record_event("forward_kernel", batch, rank, line_count, t1, t2)

      if (rank < nranks - 1) then
        t1 = MPI_Wtime()
        call pack_forward_state(first_line, line_count, forward_offset)
        t2 = MPI_Wtime()
        timing(8) = timing(8) + (t2 - t1)
        call record_event("pack_forward", batch, rank + 1, solve_count, t1, t2)
        !$omp target data use_device_addr(ycomm_sendbuf)
        t1 = MPI_Wtime()
        call MPI_Isend(ycomm_sendbuf(forward_offset + 1), solve_count, MPI_DOUBLE_COMPLEX, &
          rank + 1, TAG_FORWARD + batch, MPI_COMM_WORLD, forward_send_req(batch), ierr)
        t2 = MPI_Wtime()
        !$omp end target data
        call record_event("forward_isend_post", batch, rank + 1, solve_count, t1, t2)
      end if

      call drain_backward_batches(batch, .false., next_backward, factor_total, solve_stride, solve_total, &
        forward_send_req, backward_recv_req, backward_send_req, timing)
    end do
    call drain_backward_batches(nbatches, .true., next_backward, factor_total, solve_stride, solve_total, &
      forward_send_req, backward_recv_req, backward_send_req, timing)

    if (rank < nranks - 1) then
      t1 = MPI_Wtime()
      call MPI_Waitall(nbatches, factor_send_req, MPI_STATUSES_IGNORE, ierr)
      call MPI_Waitall(nbatches, forward_send_req, MPI_STATUSES_IGNORE, ierr)
      t2 = MPI_Wtime()
      timing(9) = timing(9) + (t2 - t1)
      call record_event("forward_waitall", 0_C_INT, rank + 1, factor_total + solve_total, t1, t2)
    end if
    if (rank > 0) then
      t1 = MPI_Wtime()
      call MPI_Waitall(nbatches, backward_send_req, MPI_STATUSES_IGNORE, ierr)
      t2 = MPI_Wtime()
      timing(9) = timing(9) + (t2 - t1)
      call record_event("backward_waitall", 0_C_INT, rank - 1, solve_total, t1, t2)
    end if

    t1 = MPI_Wtime()
    call exchange_halos()
    t2 = MPI_Wtime()
    timing(10) = timing(10) + (t2 - t1)
    call record_event("halo_exchange", 0_C_INT, -1_C_INT, 4_C_INT*nlines, t1, t2)
    timing(1) = MPI_Wtime() - t0
    deallocate(factor_send_req, forward_send_req, backward_recv_req, backward_send_req)
  end subroutine solve_once

  subroutine drain_backward_batches(completed_batch, blocking, next_backward, factor_total, solve_stride, solve_total, &
                                    forward_send_req, backward_recv_req, backward_send_req, timing)
    integer(C_INT), intent(in) :: completed_batch, factor_total, solve_stride, solve_total
    logical, intent(in) :: blocking
    integer(C_INT), intent(inout) :: next_backward
    type(MPI_Request), intent(inout) :: forward_send_req(:), backward_recv_req(:), backward_send_req(:)
    real(C_DOUBLE), intent(inout) :: timing(10)
    integer(C_INT) :: first_line, line_count, solve_count, forward_offset, backward_offset
    real(C_DOUBLE) :: t1, t2
    type(MPI_Status) :: status
    logical :: ready

    do while (next_backward <= completed_batch)
      call batch_range(next_backward, first_line, line_count)
      solve_count = 2_C_INT*line_count
      forward_offset = factor_total + (next_backward - 1_C_INT)*solve_stride
      backward_offset = factor_total + solve_total + (next_backward - 1_C_INT)*solve_stride

      if (rank < nranks - 1) then
        if (blocking) then
          t1 = MPI_Wtime()
          call MPI_Wait(backward_recv_req(next_backward), status, ierr)
          t2 = MPI_Wtime()
          timing(4) = timing(4) + (t2 - t1)
          call record_event("backward_recv_wait", next_backward, rank + 1, solve_count, t1, t2)
          ready = .true.
        else
          call MPI_Test(backward_recv_req(next_backward), ready, status, ierr)
        end if
        if (.not. ready) exit
        t1 = MPI_Wtime()
        call backward_batch_continue(first_line, line_count, backward_offset)
        t2 = MPI_Wtime()
        timing(7) = timing(7) + (t2 - t1)
        call record_event("backward_kernel", next_backward, rank, line_count, t1, t2)
      else
        t1 = MPI_Wtime()
        call backward_batch(first_line, line_count)
        t2 = MPI_Wtime()
        timing(7) = timing(7) + (t2 - t1)
        call record_event("backward_kernel", next_backward, rank, line_count, t1, t2)
      end if

      if (rank > 0) then
        if (rank < nranks - 1) then
          t1 = MPI_Wtime()
          call MPI_Wait(forward_send_req(next_backward), status, ierr)
          t2 = MPI_Wtime()
          timing(9) = timing(9) + (t2 - t1)
          call record_event("forward_send_wait", next_backward, rank + 1, solve_count, t1, t2)
        end if
        t1 = MPI_Wtime()
        call pack_backward_state(first_line, line_count, forward_offset)
        t2 = MPI_Wtime()
        timing(8) = timing(8) + (t2 - t1)
        call record_event("pack_backward", next_backward, rank - 1, solve_count, t1, t2)
        !$omp target data use_device_addr(ycomm_sendbuf)
        t1 = MPI_Wtime()
        call MPI_Isend(ycomm_sendbuf(forward_offset + 1), solve_count, MPI_DOUBLE_COMPLEX, &
          rank - 1, TAG_BACKWARD + next_backward, MPI_COMM_WORLD, backward_send_req(next_backward), ierr)
        t2 = MPI_Wtime()
        !$omp end target data
        call record_event("backward_isend_post", next_backward, rank - 1, solve_count, t1, t2)
      end if
      next_backward = next_backward + 1_C_INT
    end do
  end subroutine drain_backward_batches

  subroutine batch_range(batch, first_line, line_count)
    integer(C_INT), intent(in) :: batch
    integer(C_INT), intent(out) :: first_line, line_count
    integer(C_INT) :: base_count, remainder
    base_count = nlines/nbatches
    remainder = mod(nlines, nbatches)
    line_count = base_count
    if (batch <= remainder) line_count = line_count + 1_C_INT
    first_line = (batch - 1_C_INT)*base_count + min(batch - 1_C_INT, remainder) + 1_C_INT
  end subroutine batch_range

  subroutine apply_factor_continuation(first_line, line_count, recv_offset)
    integer(C_INT), intent(in) :: first_line, line_count, recv_offset
    integer(C_INT) :: local_line, iline
    integer(C_INT64_T) :: q
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_line_kernel(first_line, line_count)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(ds, dl, d, du, dw, ycomm_recvbuf, first_line, line_count, active_n, nlines, recv_offset) private(local_line, iline, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      q = int(recv_offset + 6_C_INT*(local_line - 1_C_INT), C_INT64_T)
      call factor_continue(ds, dl, d, du, nlines, iline, active_n, ycomm_recvbuf(q + 1_C_INT64_T), &
        ycomm_recvbuf(q + 2_C_INT64_T), ycomm_recvbuf(q + 3_C_INT64_T), ycomm_recvbuf(q + 4_C_INT64_T), &
        ycomm_recvbuf(q + 5_C_INT64_T), ycomm_recvbuf(q + 6_C_INT64_T))
    end do
    !$omp end target teams distribute parallel do
  end subroutine apply_factor_continuation

  subroutine factor_batch(first_line, line_count, open_right)
    integer(C_INT), intent(in) :: first_line, line_count
    logical, intent(in) :: open_right
    integer(C_INT) :: local_line, iline
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_line_kernel(first_line, line_count)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(ds, dl, d, du, dw, first_line, line_count, active_n, nlines, open_right) private(local_line, iline)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      if (open_right) then
        call factor_open_right(ds, dl, d, du, dw, nlines, iline, active_n)
      else
        call factor_full(ds, dl, d, du, dw, nlines, iline, active_n)
      end if
    end do
    !$omp end target teams distribute parallel do
  end subroutine factor_batch

  subroutine forward_batch(first_line, line_count)
    integer(C_INT), intent(in) :: first_line, line_count
    integer(C_INT) :: local_line, iline
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_line_kernel(first_line, line_count)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(x, ds, dl, first_line, line_count, active_n, nlines) private(local_line, iline)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      call forward_full(x, ds, dl, nlines, iline, active_n)
    end do
    !$omp end target teams distribute parallel do
  end subroutine forward_batch

  subroutine forward_batch_continue(first_line, line_count, recv_offset)
    integer(C_INT), intent(in) :: first_line, line_count, recv_offset
    integer(C_INT) :: local_line, iline
    integer(C_INT64_T) :: q
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_line_kernel(first_line, line_count)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(x, ds, dl, ycomm_recvbuf, first_line, line_count, active_n, nlines, recv_offset) private(local_line, iline, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      q = int(recv_offset + 2_C_INT*(local_line - 1_C_INT), C_INT64_T)
      call forward_continue(x, ds, dl, nlines, iline, active_n, ycomm_recvbuf(q + 1_C_INT64_T), ycomm_recvbuf(q + 2_C_INT64_T))
    end do
    !$omp end target teams distribute parallel do
  end subroutine forward_batch_continue

  subroutine backward_batch(first_line, line_count)
    integer(C_INT), intent(in) :: first_line, line_count
    integer(C_INT) :: local_line, iline
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_line_kernel(first_line, line_count)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(x, d, du, dw, first_line, line_count, active_n, nlines) private(local_line, iline)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      call backward_full(x, d, du, dw, nlines, iline, active_n)
    end do
    !$omp end target teams distribute parallel do
  end subroutine backward_batch

  subroutine backward_batch_continue(first_line, line_count, recv_offset)
    integer(C_INT), intent(in) :: first_line, line_count, recv_offset
    integer(C_INT) :: local_line, iline
    integer(C_INT64_T) :: q
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_line_kernel(first_line, line_count)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(x, d, du, dw, ycomm_recvbuf, first_line, line_count, active_n, nlines, recv_offset) private(local_line, iline, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      q = int(recv_offset + 2_C_INT*(local_line - 1_C_INT), C_INT64_T)
      call backward_continue(x, d, du, dw, nlines, iline, active_n, ycomm_recvbuf(q + 1_C_INT64_T), ycomm_recvbuf(q + 2_C_INT64_T))
    end do
    !$omp end target teams distribute parallel do
  end subroutine backward_batch_continue

  subroutine dummy_line_kernel(first_line, line_count)
    integer(C_INT), intent(in) :: first_line, line_count
    integer(C_INT) :: local_line, iline
    integer(C_INT64_T) :: p
    !$omp target teams distribute parallel do default(none) shared(x, first_line, line_count, nlines) private(local_line, iline, p)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      p = int(iline, C_INT64_T)
      x(p) = x(p) + (0.0d0, 0.0d0)
    end do
    !$omp end target teams distribute parallel do
  end subroutine dummy_line_kernel

  subroutine pack_factor_state(first_line, line_count, send_offset)
    integer(C_INT), intent(in) :: first_line, line_count, send_offset
    integer(C_INT) :: local_line, iline
    integer(C_INT64_T) :: p0, p1, q
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      !$omp target teams distribute parallel do default(none) shared(ycomm_sendbuf, line_count, send_offset) private(local_line, q)
      do local_line = 1_C_INT, line_count
        q = int(send_offset + 6_C_INT*(local_line - 1_C_INT), C_INT64_T)
        ycomm_sendbuf(q + 1_C_INT64_T) = (1.0d0, 0.0d0)
        ycomm_sendbuf(q + 2_C_INT64_T) = (0.0d0, 0.0d0)
        ycomm_sendbuf(q + 3_C_INT64_T) = (0.0d0, 0.0d0)
        ycomm_sendbuf(q + 4_C_INT64_T) = (1.0d0, 0.0d0)
        ycomm_sendbuf(q + 5_C_INT64_T) = (0.0d0, 0.0d0)
        ycomm_sendbuf(q + 6_C_INT64_T) = (0.0d0, 0.0d0)
      end do
      !$omp end target teams distribute parallel do
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(d, du, dw, ycomm_sendbuf, first_line, line_count, active_n, nlines, send_offset) private(local_line, iline, p0, p1, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      p0 = int(iline, C_INT64_T) + int(active_n - 2_C_INT, C_INT64_T)*int(nlines, C_INT64_T)
      p1 = p0 + int(nlines, C_INT64_T)
      q = int(send_offset + 6_C_INT*(local_line - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = d(p0)
      ycomm_sendbuf(q + 2_C_INT64_T) = du(p0)
      ycomm_sendbuf(q + 3_C_INT64_T) = dw(p0)
      ycomm_sendbuf(q + 4_C_INT64_T) = d(p1)
      ycomm_sendbuf(q + 5_C_INT64_T) = du(p1)
      ycomm_sendbuf(q + 6_C_INT64_T) = dw(p1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine pack_factor_state

  subroutine pack_forward_state(first_line, line_count, send_offset)
    integer(C_INT), intent(in) :: first_line, line_count, send_offset
    integer(C_INT) :: local_line, iline
    integer(C_INT64_T) :: p0, p1, q
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_pack_solve(line_count, send_offset)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(x, ycomm_sendbuf, first_line, line_count, active_n, nlines, send_offset) private(local_line, iline, p0, p1, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      p0 = int(iline, C_INT64_T) + int(active_n - 2_C_INT, C_INT64_T)*int(nlines, C_INT64_T)
      p1 = p0 + int(nlines, C_INT64_T)
      q = int(send_offset + 2_C_INT*(local_line - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = x(p0)
      ycomm_sendbuf(q + 2_C_INT64_T) = x(p1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine pack_forward_state

  subroutine pack_backward_state(first_line, line_count, send_offset)
    integer(C_INT), intent(in) :: first_line, line_count, send_offset
    integer(C_INT) :: local_line, iline
    integer(C_INT64_T) :: p0, p1, q
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_pack_solve(line_count, send_offset)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(x, ycomm_sendbuf, first_line, line_count, nlines, send_offset) private(local_line, iline, p0, p1, q)
    do local_line = 1_C_INT, line_count
      iline = first_line + local_line - 1_C_INT
      p0 = int(iline, C_INT64_T)
      p1 = p0 + int(nlines, C_INT64_T)
      q = int(send_offset + 2_C_INT*(local_line - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = x(p0)
      ycomm_sendbuf(q + 2_C_INT64_T) = x(p1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine pack_backward_state

  subroutine dummy_pack_solve(line_count, send_offset)
    integer(C_INT), intent(in) :: line_count, send_offset
    integer(C_INT) :: local_line
    integer(C_INT64_T) :: q
    !$omp target teams distribute parallel do default(none) shared(ycomm_sendbuf, line_count, send_offset) private(local_line, q)
    do local_line = 1_C_INT, line_count
      q = int(send_offset + 2_C_INT*(local_line - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = (0.0d0, 0.0d0)
      ycomm_sendbuf(q + 2_C_INT64_T) = (0.0d0, 0.0d0)
    end do
    !$omp end target teams distribute parallel do
  end subroutine dummy_pack_solve

  subroutine exchange_halos()
    type(MPI_Request) :: req(4)
    integer :: nreq
    nreq = 0
    if (rank > 0) then
      nreq = nreq + 1
      !$omp target data use_device_addr(left_halo)
      call MPI_Irecv(left_halo(1), 2_C_INT*nlines, MPI_DOUBLE_COMPLEX, rank - 1, TAG_HALO_LOW, MPI_COMM_WORLD, req(nreq), ierr)
      !$omp end target data
    end if
    if (rank < nranks - 1) then
      nreq = nreq + 1
      !$omp target data use_device_addr(right_halo)
      call MPI_Irecv(right_halo(1), 2_C_INT*nlines, MPI_DOUBLE_COMPLEX, rank + 1, TAG_HALO_HIGH, MPI_COMM_WORLD, req(nreq), ierr)
      !$omp end target data
    end if
    if (rank > 0) then
      call pack_lower_halo(0_C_INT)
      nreq = nreq + 1
      !$omp target data use_device_addr(ycomm_sendbuf)
      call MPI_Isend(ycomm_sendbuf(1), 2_C_INT*nlines, MPI_DOUBLE_COMPLEX, rank - 1, TAG_HALO_HIGH, MPI_COMM_WORLD, req(nreq), ierr)
      !$omp end target data
    end if
    if (rank < nranks - 1) then
      call pack_upper_halo(2_C_INT*nlines)
      nreq = nreq + 1
      !$omp target data use_device_addr(ycomm_sendbuf)
      call MPI_Isend(ycomm_sendbuf(2_C_INT*nlines + 1), 2_C_INT*nlines, MPI_DOUBLE_COMPLEX, rank + 1, TAG_HALO_LOW, MPI_COMM_WORLD, req(nreq), ierr)
      !$omp end target data
    end if
    if (nreq > 0) call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
  end subroutine exchange_halos

  subroutine pack_lower_halo(send_offset)
    integer(C_INT), intent(in) :: send_offset
    integer(C_INT) :: iline
    integer(C_INT64_T) :: q
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_pack_solve(nlines, send_offset)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(x, ycomm_sendbuf, nlines, send_offset) private(iline, q)
    do iline = 1_C_INT, nlines
      q = int(send_offset + 2_C_INT*(iline - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = x(iline)
      ycomm_sendbuf(q + 2_C_INT64_T) = x(iline + nlines)
    end do
    !$omp end target teams distribute parallel do
  end subroutine pack_lower_halo

  subroutine pack_upper_halo(send_offset)
    integer(C_INT), intent(in) :: send_offset
    integer(C_INT) :: iline
    integer(C_INT64_T) :: p0, p1, q
    if (trim(kernel_mode) == "none") return
    if (trim(kernel_mode) == "dummy") then
      call dummy_pack_solve(nlines, send_offset)
      return
    end if
    !$omp target teams distribute parallel do default(none) shared(x, ycomm_sendbuf, active_n, nlines, send_offset) private(iline, p0, p1, q)
    do iline = 1_C_INT, nlines
      p0 = int(iline, C_INT64_T) + int(active_n - 2_C_INT, C_INT64_T)*int(nlines, C_INT64_T)
      p1 = p0 + int(nlines, C_INT64_T)
      q = int(send_offset + 2_C_INT*(iline - 1_C_INT), C_INT64_T)
      ycomm_sendbuf(q + 1_C_INT64_T) = x(p0)
      ycomm_sendbuf(q + 2_C_INT64_T) = x(p1)
    end do
    !$omp end target teams distribute parallel do
  end subroutine pack_upper_halo

  subroutine solver_workspace_bytes(active_n_in, nlines_in, nranks_in, nbytes)
    integer(C_INT), intent(in) :: active_n_in, nlines_in, nranks_in
    integer(C_SIZE_T), intent(out) :: nbytes
    integer(C_SIZE_T) :: nall, count, offset
    nall = int(active_n_in, C_SIZE_T)*int(nlines_in, C_SIZE_T)
    count = int(nlines_in, C_SIZE_T)
    offset = 0_C_SIZE_T
    call skip_complex(offset, 5_C_SIZE_T*nall)
    call skip_complex(offset, nall)
    call skip_complex(offset, count)
    call skip_complex(offset, count)
    call skip_complex(offset, count)
    call skip_complex(offset, count)
    call skip_real(offset, 5_C_SIZE_T*count)
    call skip_real(offset, 5_C_SIZE_T*count)
    call skip_real(offset, 5_C_SIZE_T*count)
    call skip_real(offset, 5_C_SIZE_T*count)
    call skip_complex(offset, count)
    call skip_complex(offset, count)
    call skip_real(offset, 4_C_SIZE_T*count)
    call skip_real(offset, 4_C_SIZE_T*count)
    call skip_complex(offset, nall)
    call skip_complex(offset, nall)
    call skip_complex(offset, nall)
    call skip_complex(offset, nall)
    call skip_complex(offset, nall)
    call skip_complex(offset, nall)
    if (nranks_in > 1) then
      call skip_complex(offset, 20_C_SIZE_T*count)
      call skip_complex(offset, 2_C_SIZE_T*count)
      call skip_complex(offset, 2_C_SIZE_T*count)
      call skip_complex(offset, 4_C_SIZE_T*int(nranks_in, C_SIZE_T)*count)
    end if
    nbytes = offset
  end subroutine solver_workspace_bytes

  subroutine skip_complex(offset, count)
    integer(C_SIZE_T), intent(inout) :: offset
    integer(C_SIZE_T), intent(in) :: count
    offset = align_offset(offset + count*int(C_SIZEOF((0.0_C_DOUBLE, 0.0_C_DOUBLE)), C_SIZE_T))
  end subroutine skip_complex

  subroutine skip_real(offset, count)
    integer(C_SIZE_T), intent(inout) :: offset
    integer(C_SIZE_T), intent(in) :: count
    offset = align_offset(offset + count*int(C_SIZEOF(0.0_C_DOUBLE), C_SIZE_T))
  end subroutine skip_real

  integer(C_SIZE_T) function align_offset(offset)
    integer(C_SIZE_T), intent(in) :: offset
    integer(C_SIZE_T) :: rem
    rem = modulo(offset, WORKSPACE_ALIGNMENT)
    if (rem == 0_C_SIZE_T) then
      align_offset = offset
    else
      align_offset = offset + WORKSPACE_ALIGNMENT - rem
    end if
  end function align_offset

  subroutine choose_base_index(raw, base_index)
    integer(C_INT8_T), target, intent(in) :: raw(:)
    integer(C_SIZE_T), intent(out) :: base_index
    integer(C_INTPTR_T) :: addr
    integer(C_SIZE_T) :: rem
    addr = transfer(c_loc(raw(1)), addr)
    rem = modulo(int(addr, C_SIZE_T), WORKSPACE_ALIGNMENT)
    if (rem == 0_C_SIZE_T) then
      base_index = 1_C_SIZE_T
    else
      base_index = 1_C_SIZE_T + WORKSPACE_ALIGNMENT - rem
    end if
  end subroutine choose_base_index

  subroutine factor_full(ds, dl, d, du, dw, stride, first, n)
    !$omp declare target
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i
    integer(C_INT64_T) :: p, p1, p2
    complex(C_DOUBLE_COMPLEX) :: factor
    if (n <= 0) return
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
    d(p1) = 1.0d0/d(p1)
  end subroutine factor_full

  subroutine factor_open_right(ds, dl, d, du, dw, stride, first, n)
    !$omp declare target
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i
    integer(C_INT64_T) :: p, p1, p2
    complex(C_DOUBLE_COMPLEX) :: factor
    if (n <= 0) return
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
    du(p1) = du(p1) - factor*dw(p)
    d(p1) = 1.0d0/d(p1)
  end subroutine factor_open_right

  subroutine factor_continue(ds, dl, d, du, stride, first, n, prev0_d, prev0_du, prev0_dw, prev1_d, prev1_du, prev1_dw)
    !$omp declare target
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
  end subroutine factor_continue

  subroutine forward_full(rhs, ds, dl, stride, first, n)
    !$omp declare target
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: ds(:), dl(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i
    integer(C_INT64_T) :: p
    if (n >= 2) then
      p = int(first, C_INT64_T) + int(stride, C_INT64_T)
      rhs(p) = rhs(p) - dl(p)*rhs(first)
    end if
    do i = 2, n - 1
      p = int(first, C_INT64_T) + int(i, C_INT64_T)*int(stride, C_INT64_T)
      rhs(p) = rhs(p) - ds(p)*rhs(int(first, C_INT64_T) + int(i - 2, C_INT64_T)*int(stride, C_INT64_T)) - &
        dl(p)*rhs(int(first, C_INT64_T) + int(i - 1, C_INT64_T)*int(stride, C_INT64_T))
    end do
  end subroutine forward_full

  subroutine forward_continue(rhs, ds, dl, stride, first, n, prev0_rhs, prev1_rhs)
    !$omp declare target
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: ds(:), dl(:)
    integer(C_INT), intent(in) :: stride, first, n
    complex(C_DOUBLE_COMPLEX), intent(in) :: prev0_rhs, prev1_rhs
    integer(C_INT) :: i
    integer(C_INT64_T) :: p0, p1, p
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
  end subroutine forward_continue

  subroutine backward_full(rhs, d, du, dw, stride, first, n)
    !$omp declare target
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i
    integer(C_INT64_T) :: p
    p = int(first, C_INT64_T) + int(n - 1, C_INT64_T)*int(stride, C_INT64_T)
    rhs(p) = rhs(p)*d(p)
    p = int(first, C_INT64_T) + int(n - 2, C_INT64_T)*int(stride, C_INT64_T)
    rhs(p) = (rhs(p) - du(p)*rhs(p + int(stride, C_INT64_T)))*d(p)
    do i = n - 3, 0, -1
      p = int(first, C_INT64_T) + int(i, C_INT64_T)*int(stride, C_INT64_T)
      rhs(p) = (rhs(p) - du(p)*rhs(p + int(stride, C_INT64_T)) - dw(p)*rhs(p + 2_C_INT64_T*int(stride, C_INT64_T)))*d(p)
    end do
  end subroutine backward_full

  subroutine backward_continue(rhs, d, du, dw, stride, first, n, next0_rhs, next1_rhs)
    !$omp declare target
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    complex(C_DOUBLE_COMPLEX), intent(in) :: next0_rhs, next1_rhs
    integer(C_INT) :: i
    integer(C_INT64_T) :: p, p_last, p_prev
    p_last = int(first, C_INT64_T) + int(n - 1, C_INT64_T)*int(stride, C_INT64_T)
    rhs(p_last) = (rhs(p_last) - du(p_last)*next0_rhs - dw(p_last)*next1_rhs)*d(p_last)
    p_prev = p_last - int(stride, C_INT64_T)
    rhs(p_prev) = (rhs(p_prev) - du(p_prev)*rhs(p_last) - dw(p_prev)*next0_rhs)*d(p_prev)
    do i = n - 3, 0, -1
      p = int(first, C_INT64_T) + int(i, C_INT64_T)*int(stride, C_INT64_T)
      rhs(p) = (rhs(p) - du(p)*rhs(p + int(stride, C_INT64_T)) - dw(p)*rhs(p + 2_C_INT64_T*int(stride, C_INT64_T)))*d(p)
    end do
  end subroutine backward_continue
end program bench_y_lu_pipeline_autotune
F90

"${MPIF90}" -O2 -mp=gpu ${GPU_FLAGS} -Minfo=mp -o bench_y_lu_pipeline_autotune bench_y_lu_pipeline_autotune.f90

{
  echo "host=$(hostname -f 2>/dev/null || hostname)"
  echo "date=$(date -Is)"
  echo "np=${NP}"
  echo "iters=${ITERS}"
  echo "warmup=${WARMUP}"
  echo "active_n=${ACTIVE_N}"
  echo "nlines=${NLINES}"
  echo "nx=${NX}"
  echo "ny=${NY}"
  echo "nz=${NZ}"
  echo "npxz=${NPXZ}"
  echo "npy=${NPY}"
  echo "batches=${BATCHES}"
  echo "storage_mode=${STORAGE_MODE}"
  echo "kernel_mode=${KERNEL_MODE}"
  echo "detail=${DETAIL}"
  echo "sweep_batches=${SWEEP_BATCHES}"
  echo "mpif90=$("${MPIF90}" --version 2>&1 | head -1)"
  echo "mpirun=$("${MPIRUN}" --version 2>&1 | head -1)"
  echo "mpirun_args=${MPIRUN_ARGS}"
  echo "gpu_flags=${GPU_FLAGS}"
  echo "profile=${PROFILE}"
  nvidia-smi -L 2>/dev/null || true
} > run_info.txt

if [[ "${PROFILE}" == "1" ]]; then
  "${MPIRUN}" "${MPIRUN_EXTRA_ARGS[@]}" -np "${NP}" bash -lc '
    if [[ -f /etc/profile.d/lmod.sh ]]; then
      . /etc/profile.d/lmod.sh
      module load "'"${NVHPC_MODULE}"'"
    fi
    host=$(hostname -s)
    rank=${OMPI_COMM_WORLD_RANK:-${PMIX_RANK:-${PMI_RANK:-${SLURM_PROCID:-0}}}}
    out="nsys_lu_${host}_rank_${rank}"
    nsys profile --trace="'"${TRACE}"'" --sample=none --cpuctxsw=none --stats=true \
      --force-overwrite=true --export=sqlite -o "${out}" \
      ./bench_y_lu_pipeline_autotune "'"${ITERS}"'" "'"${WARMUP}"'" "'"${ACTIVE_N}"'" "'"${NLINES}"'" "'"${BATCHES}"'" "'"${STORAGE_MODE}"'" "'"${NX}"'" "'"${NY}"'" "'"${NZ}"'" "'"${NPXZ}"'" "'"${NPY}"'" "'"${KERNEL_MODE}"'" "'"${DETAIL}"'"
  ' > stdout.log 2> nsys.log
else
  if [[ "${SWEEP_BATCHES}" == "1" ]]; then
    : > stdout.log
    : > stderr.log
    for b in 1 2 4 8 16; do
      {
        echo "=== batches=${b} kernel_mode=${KERNEL_MODE} storage_mode=${STORAGE_MODE} ==="
        "${MPIRUN}" "${MPIRUN_EXTRA_ARGS[@]}" -np "${NP}" ./bench_y_lu_pipeline_autotune "${ITERS}" "${WARMUP}" "${ACTIVE_N}" "${NLINES}" "${b}" "${STORAGE_MODE}" "${NX}" "${NY}" "${NZ}" "${NPXZ}" "${NPY}" "${KERNEL_MODE}" "${DETAIL}"
      } >> stdout.log 2>> stderr.log
    done
  else
    "${MPIRUN}" "${MPIRUN_EXTRA_ARGS[@]}" -np "${NP}" ./bench_y_lu_pipeline_autotune "${ITERS}" "${WARMUP}" "${ACTIVE_N}" "${NLINES}" "${BATCHES}" "${STORAGE_MODE}" "${NX}" "${NY}" "${NZ}" "${NPXZ}" "${NPY}" "${KERNEL_MODE}" "${DETAIL}" \
      > stdout.log 2> stderr.log
  fi
fi

echo "Wrote ${OUT_DIR}"
