#!/usr/bin/env bash
set -euo pipefail

NP="${NP:-8}"
ITERS="${ITERS:-50}"
WARMUP="${WARMUP:-5}"
MIN_BYTES="${MIN_BYTES:-512}"
MAX_BYTES="${MAX_BYTES:-67108864}"
COMM_MODE="${COMM_MODE:-nonblocking}"
PIPELINE_BATCHES="${PIPELINE_BATCHES:-4}"
OUT_DIR="${OUT_DIR:-omp_device_mpi_pingpong_$(date +%Y%m%d_%H%M%S)}"
PROFILE="${PROFILE:-0}"
TRACE="${TRACE:-nvtx,cuda}"
NVHPC_MODULE="${NVHPC_MODULE:-toolkit/nvidia-hpc-sdk/25.3}"
MPIF90="${MPIF90:-mpifort}"
MPIRUN="${MPIRUN:-mpirun}"
GPU_FLAGS="${GPU_FLAGS:--gpu=cc80}"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

if [[ -f /etc/profile.d/lmod.sh ]]; then
  # shellcheck disable=SC1091
  . /etc/profile.d/lmod.sh
  module load "${NVHPC_MODULE}" || true
fi

cat > omp_device_mpi_pingpong.f90 <<'F90'
program omp_device_mpi_pingpong
  use, intrinsic :: iso_c_binding
  use mpi_f08
  implicit none

  integer :: ierr, rank, nranks, peer, src_rank, dst_rank, warmup, iters, arg_count
  integer :: pipeline_batches
  integer(C_SIZE_T) :: min_bytes, max_bytes, nbytes, alloc_bytes
  integer(C_INT64_T) :: nelems, i
  integer(C_INT64_T), allocatable, target :: sendbuf(:), recvbuf(:)
  real(C_DOUBLE) :: t0, t1, elapsed, bw_gbs
  logical :: even_rank
  character(64) :: arg, comm_mode
  type(MPI_Status) :: status
  type(MPI_Request) :: req(2)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

  iters = 50
  warmup = 5
  min_bytes = 512_C_SIZE_T
  max_bytes = 67108864_C_SIZE_T
  comm_mode = "nonblocking"
  pipeline_batches = 4

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
    read(arg, *) min_bytes
  end if
  if (arg_count >= 4) then
    call get_command_argument(4, arg)
    read(arg, *) max_bytes
  end if
  if (arg_count >= 5) then
    call get_command_argument(5, comm_mode)
    comm_mode = adjustl(comm_mode)
  end if
  if (arg_count >= 6) then
    call get_command_argument(6, arg)
    read(arg, *) pipeline_batches
  end if

  alloc_bytes = max_bytes
  if (trim(comm_mode) == "pipeline") then
    alloc_bytes = max_bytes*int(max(1, pipeline_batches), C_SIZE_T)
  end if
  nelems = max(1_C_INT64_T, int((alloc_bytes + 7_C_SIZE_T)/8_C_SIZE_T, C_INT64_T))
  allocate(sendbuf(nelems), recvbuf(nelems))

  !$omp target enter data map(alloc: sendbuf, recvbuf)
  !$omp target teams distribute parallel do default(none) shared(sendbuf, recvbuf, nelems) private(i)
  do i = 1_C_INT64_T, nelems
    sendbuf(i) = i
    recvbuf(i) = 0_C_INT64_T
  end do
  !$omp end target teams distribute parallel do

  if (rank == 0) then
    print '(a)', "OpenMP target MPI device-buffer ping-pong"
    print '(a,i0,a,i0)', "iters=", iters, " warmup=", warmup
    print '(a,a)', "comm_mode=", trim(comm_mode)
    if (trim(comm_mode) == "pipeline") then
      print '(a,i0)', "pipeline_batches=", pipeline_batches
      print '(a)', "bytes, batches, rank, prev_rank, next_rank, total_us, recv_wait_us, send_wait_us, pack_us, chain_GBps"
    else
      print '(a)', "bytes, src_rank, dst_rank, active_rank, one_way_us, bidirectional_GBps"
    end if
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (trim(comm_mode) == "pipeline") then
    call run_pipeline_sweep()
  else
    call run_pair_sweep()
  end if

  !$omp target exit data map(delete: sendbuf, recvbuf)
  deallocate(sendbuf, recvbuf)
  call MPI_Finalize(ierr)

contains
  subroutine run_pair_sweep()
    nbytes = min_bytes
    do while (nbytes <= max_bytes)
      nelems = max(1_C_INT64_T, int((nbytes + 7_C_SIZE_T)/8_C_SIZE_T, C_INT64_T))

      do src_rank = 0, nranks - 2
        do dst_rank = src_rank + 1, nranks - 1
          peer = -1
          even_rank = .false.
          if (rank == src_rank) then
            peer = dst_rank
            even_rank = .true.
          else if (rank == dst_rank) then
            peer = src_rank
          end if

          call exchange_loop(warmup, nbytes, nelems)
          call MPI_Barrier(MPI_COMM_WORLD, ierr)
          t0 = MPI_Wtime()
          call exchange_loop(iters, nbytes, nelems)
          call MPI_Barrier(MPI_COMM_WORLD, ierr)
          t1 = MPI_Wtime()

          if (peer >= 0) then
            elapsed = (t1 - t0)/real(max(1, iters), C_DOUBLE)
            bw_gbs = 2.0_C_DOUBLE*real(nbytes, C_DOUBLE)/(elapsed*1.0e9_C_DOUBLE)
            write (*, '(i0,", ",i0,", ",i0,", ",i0,", ",f12.3,", ",f12.3)') &
              int(nbytes), src_rank, dst_rank, rank, 0.5_C_DOUBLE*elapsed*1.0e6_C_DOUBLE, bw_gbs
          end if
        end do
      end do

      nbytes = nbytes*2_C_SIZE_T
    end do
  end subroutine run_pair_sweep

  subroutine run_pipeline_sweep()
    integer :: print_rank
    real(C_DOUBLE) :: avg_total, avg_recv, avg_send, avg_pack, chain_gbs

    nbytes = min_bytes
    do while (nbytes <= max_bytes)
      nelems = max(1_C_INT64_T, int((nbytes + 7_C_SIZE_T)/8_C_SIZE_T, C_INT64_T))

      call pipeline_loop(warmup, nbytes, nelems, avg_total, avg_recv, avg_send, avg_pack)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call pipeline_loop(iters, nbytes, nelems, avg_total, avg_recv, avg_send, avg_pack)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if (avg_total > 0.0_C_DOUBLE) then
        chain_gbs = real(nbytes, C_DOUBLE)*real(max(1, pipeline_batches), C_DOUBLE)/(avg_total*1.0e9_C_DOUBLE)
      else
        chain_gbs = 0.0_C_DOUBLE
      end if

      do print_rank = 0, nranks - 1
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        if (rank == print_rank) then
          write (*, '(i0,", ",i0,", ",i0,", ",i0,", ",i0,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3,", ",f12.3)') &
            int(nbytes), pipeline_batches, rank, merge(rank - 1, -1, rank > 0), merge(rank + 1, -1, rank < nranks - 1), &
            avg_total*1.0e6_C_DOUBLE, avg_recv*1.0e6_C_DOUBLE, avg_send*1.0e6_C_DOUBLE, &
            avg_pack*1.0e6_C_DOUBLE, chain_gbs
        end if
      end do
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      nbytes = nbytes*2_C_SIZE_T
    end do
  end subroutine run_pipeline_sweep

  subroutine exchange_loop(nloop, nbytes, nelems)
    integer, intent(in) :: nloop
    integer(C_SIZE_T), intent(in) :: nbytes
    integer(C_INT64_T), intent(in) :: nelems
    integer :: k, count

    if (peer < 0) return
    count = int(nelems)
    do k = 1, nloop
      !$omp target data use_device_addr(sendbuf, recvbuf)
      select case (trim(comm_mode))
      case ("blocking")
        if (even_rank) then
          call MPI_Send(sendbuf, count, MPI_INTEGER8, peer, 1001, MPI_COMM_WORLD, ierr)
          call MPI_Recv(recvbuf, count, MPI_INTEGER8, peer, 1002, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_Recv(recvbuf, count, MPI_INTEGER8, peer, 1001, MPI_COMM_WORLD, status, ierr)
          call MPI_Send(sendbuf, count, MPI_INTEGER8, peer, 1002, MPI_COMM_WORLD, ierr)
        end if
      case ("nonblocking")
        call MPI_Irecv(recvbuf, count, MPI_INTEGER8, peer, 2001, MPI_COMM_WORLD, req(1), ierr)
        call MPI_Isend(sendbuf, count, MPI_INTEGER8, peer, 2001, MPI_COMM_WORLD, req(2), ierr)
        call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
      case default
        if (rank == 0) print *, "unknown comm mode: ", trim(comm_mode)
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end select
      !$omp end target data
    end do
  end subroutine exchange_loop

  subroutine pipeline_loop(nloop, nbytes, nelems, avg_total, avg_recv, avg_send, avg_pack)
    integer, intent(in) :: nloop
    integer(C_SIZE_T), intent(in) :: nbytes
    integer(C_INT64_T), intent(in) :: nelems
    real(C_DOUBLE), intent(out) :: avg_total, avg_recv, avg_send, avg_pack

    integer :: b, k, count, prev_rank, next_rank
    integer(C_INT64_T) :: offset, i
    real(C_DOUBLE) :: t_start, t_end, t_wait0, t_wait1, t_pack0, t_pack1
    real(C_DOUBLE) :: sum_total, sum_recv, sum_send, sum_pack
    type(MPI_Request), allocatable :: recv_req(:), send_req(:)

    count = int(nelems)
    prev_rank = rank - 1
    next_rank = rank + 1
    sum_total = 0.0_C_DOUBLE
    sum_recv = 0.0_C_DOUBLE
    sum_send = 0.0_C_DOUBLE
    sum_pack = 0.0_C_DOUBLE
    allocate(recv_req(max(1, pipeline_batches)), send_req(max(1, pipeline_batches)))

    do k = 1, nloop
      recv_req = MPI_REQUEST_NULL
      send_req = MPI_REQUEST_NULL
      t_start = MPI_Wtime()

      !$omp target data use_device_addr(sendbuf, recvbuf)
      if (rank > 0) then
        do b = 1, pipeline_batches
          offset = int(b - 1, C_INT64_T)*nelems
          call MPI_Irecv(recvbuf(offset + 1_C_INT64_T), count, MPI_INTEGER8, prev_rank, &
            3000 + b, MPI_COMM_WORLD, recv_req(b), ierr)
        end do
      end if

      do b = 1, pipeline_batches
        offset = int(b - 1, C_INT64_T)*nelems
        if (rank > 0) then
          t_wait0 = MPI_Wtime()
          call MPI_Wait(recv_req(b), MPI_STATUS_IGNORE, ierr)
          t_wait1 = MPI_Wtime()
          sum_recv = sum_recv + (t_wait1 - t_wait0)
        end if

        t_pack0 = MPI_Wtime()
        !$omp target teams distribute parallel do default(none) firstprivate(offset, rank, b) &
        !$omp& shared(sendbuf, recvbuf, nelems) private(i)
        do i = 1_C_INT64_T, nelems
          sendbuf(offset + i) = recvbuf(offset + i) + int(rank + b, C_INT64_T) + i
        end do
        !$omp end target teams distribute parallel do
        t_pack1 = MPI_Wtime()
        sum_pack = sum_pack + (t_pack1 - t_pack0)

        if (rank < nranks - 1) then
          call MPI_Isend(sendbuf(offset + 1_C_INT64_T), count, MPI_INTEGER8, next_rank, &
            3000 + b, MPI_COMM_WORLD, send_req(b), ierr)
        end if
      end do

      if (rank < nranks - 1) then
        t_wait0 = MPI_Wtime()
        call MPI_Waitall(pipeline_batches, send_req, MPI_STATUSES_IGNORE, ierr)
        t_wait1 = MPI_Wtime()
        sum_send = sum_send + (t_wait1 - t_wait0)
      end if
      !$omp end target data

      t_end = MPI_Wtime()
      sum_total = sum_total + (t_end - t_start)
    end do

    avg_total = sum_total/real(max(1, nloop), C_DOUBLE)
    avg_recv = sum_recv/real(max(1, nloop), C_DOUBLE)
    avg_send = sum_send/real(max(1, nloop), C_DOUBLE)
    avg_pack = sum_pack/real(max(1, nloop), C_DOUBLE)
    deallocate(recv_req, send_req)
  end subroutine pipeline_loop
end program omp_device_mpi_pingpong
F90

"${MPIF90}" -O2 -mp=gpu ${GPU_FLAGS} -Minfo=mp -o omp_device_mpi_pingpong omp_device_mpi_pingpong.f90

{
  echo "host=$(hostname -f 2>/dev/null || hostname)"
  echo "date=$(date -Is)"
  echo "np=${NP}"
  echo "mpif90=$("${MPIF90}" --version 2>&1 | head -1)"
  echo "mpirun=$("${MPIRUN}" --version 2>&1 | head -1)"
  echo "gpu_flags=${GPU_FLAGS}"
  echo "profile=${PROFILE}"
  echo "comm_mode=${COMM_MODE}"
  echo "pipeline_batches=${PIPELINE_BATCHES}"
  nvidia-smi -L 2>/dev/null || true
} > run_info.txt

if [[ "${PROFILE}" == "1" ]]; then
  "${MPIRUN}" -np "${NP}" bash -lc '
    if [[ -f /etc/profile.d/lmod.sh ]]; then
      . /etc/profile.d/lmod.sh
      module load "'"${NVHPC_MODULE}"'"
    fi
    host=$(hostname -s)
    rank=${OMPI_COMM_WORLD_RANK:-${PMIX_RANK:-${PMI_RANK:-${SLURM_PROCID:-0}}}}
    out="nsys_${host}_rank_${rank}"
    nsys profile --trace="'"${TRACE}"'" --sample=none --cpuctxsw=none --stats=true \
      --force-overwrite=true --export=sqlite -o "${out}" \
      ./omp_device_mpi_pingpong "'"${ITERS}"'" "'"${WARMUP}"'" "'"${MIN_BYTES}"'" "'"${MAX_BYTES}"'" "'"${COMM_MODE}"'" "'"${PIPELINE_BATCHES}"'"
  ' > stdout.log 2> nsys.log
else
  "${MPIRUN}" -np "${NP}" ./omp_device_mpi_pingpong "${ITERS}" "${WARMUP}" "${MIN_BYTES}" "${MAX_BYTES}" "${COMM_MODE}" "${PIPELINE_BATCHES}" \
    > stdout.log 2> stderr.log
fi

echo "Wrote ${OUT_DIR}"
