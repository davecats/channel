#include "header.h"

program schur_alltoallv_repro
  use, intrinsic :: iso_c_binding
  use mpi_f08
  use omp_lib
  implicit none

  integer(C_INT), parameter :: row_width = 20_C_INT
  integer(C_INT), parameter :: value_width = 8_C_INT

  integer :: ierr, rank, nprocs, comm_size
  integer :: num_dev, dev
  integer(C_INT) :: nx, ny, nz, npy, npxz, ipy, ipxz
  integer(C_INT) :: nxpp, nxB, nlines_z, nlines
  integer(C_INT) :: ny0, nyN, active_n
  integer(C_INT) :: arity, owned_first, owned_count, max_line_count
  integer(C_INT) :: row_send_elems, row_recv_elems, value_send_elems, value_recv_elems
  integer(C_INT) :: send_capacity, recv_capacity
  integer(C_INT) :: iters, warmup, wrap_use_device_addr
  integer(C_INT) :: use_alltoall
  integer(C_INT) :: iter, dest, src, local_line, irow, k, child, row0
  integer(C_INT) :: line_count, first_line, global_line, offset
  integer(C_INT) :: row_send_offset, row_recv_offset, value_recv_offset
  integer :: irank
  integer, allocatable, target :: row_send_counts(:), row_send_displs(:), row_recv_counts(:), row_recv_displs(:)
  integer, allocatable, target :: value_send_counts(:), value_send_displs(:), value_recv_counts(:), value_recv_displs(:)
  integer(C_INT), allocatable, target :: line_first(:), line_counts(:)
  complex(C_DOUBLE_COMPLEX), pointer :: ycomm_sendbuf(:), ycomm_recvbuf(:)
  complex(C_DOUBLE_COMPLEX), allocatable, target :: reduced_rows_send(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, target :: reduced_rhs(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, target :: left_interface_values(:, :), right_interface_values(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, target :: s_rows(:, :, :), s_values(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, target :: s_recover_basis(:, :, :, :)
  type(C_PTR) :: ycomm_sendptr, ycomm_recvptr
  type(MPI_Comm) :: comm_y
  real(C_DOUBLE) :: t0, t1, total_rows, total_values
  complex(C_DOUBLE_COMPLEX) :: prev1, prev2, next1, next2
  complex(C_DOUBLE_COMPLEX) :: recovered(4*16)
  character(len=32) :: exchange_mode, alloc_mode

#if defined(HAVE_HIP)
  interface
    subroutine roctxRangePushA(message) bind(C, name="roctxRangePushA")
      import :: C_CHAR
      character(C_CHAR) :: message(*)
    end subroutine roctxRangePushA

    subroutine roctxRangePop() bind(C, name="roctxRangePop")
    end subroutine roctxRangePop

    integer(C_INT) function hipMalloc(ptr, bytes) bind(C, name="hipMalloc")
      import :: C_INT, C_PTR, C_SIZE_T
      type(C_PTR), intent(out) :: ptr
      integer(C_SIZE_T), value :: bytes
    end function hipMalloc

    integer(C_INT) function hipFree(ptr) bind(C, name="hipFree")
      import :: C_INT, C_PTR
      type(C_PTR), value :: ptr
    end function hipFree
  end interface
#endif

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  nx = 511_C_INT
  ny = 300_C_INT
  nz = 513_C_INT
  npy = 4_C_INT
  iters = 20_C_INT
  warmup = 2_C_INT
  wrap_use_device_addr = 0_C_INT
  exchange_mode = "alltoallv"
  alloc_mode = "omp"
  call getenv_int("NX", nx)
  call getenv_int("NY", ny)
  call getenv_int("NZ", nz)
  call getenv_int("NPY", npy)
  call getenv_int("ITERS", iters)
  call getenv_int("WARMUP", warmup)
  call getenv_int("WRAP_USE_DEVICE_ADDR", wrap_use_device_addr)
  call getenv_string("EXCHANGE", exchange_mode)
  call getenv_string("ALLOC", alloc_mode)
  call lowercase(exchange_mode)
  call lowercase(alloc_mode)
  select case (trim(exchange_mode))
  case ("alltoallv")
    use_alltoall = 0_C_INT
  case ("alltoall")
    use_alltoall = 1_C_INT
  case default
    if (rank == 0) print *, "EXCHANGE must be alltoallv or alltoall"
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end select
  select case (trim(alloc_mode))
  case ("omp", "hip")
  case default
    if (rank == 0) print *, "ALLOC must be omp or hip"
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end select

  if (npy < 1_C_INT .or. mod(int(nprocs, C_INT), npy) /= 0_C_INT) then
    if (rank == 0) print *, "nprocs must be a multiple of NPY"
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
  npxz = int(nprocs, C_INT)/npy
  ipy = int(rank, C_INT)/npxz
  ipxz = mod(int(rank, C_INT), npxz)
  if (npxz /= 1_C_INT) then
    if (rank == 0) print *, "This reproducer is intended for np=npy=4; got npxz=", npxz
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  num_dev = omp_get_num_devices()
  if (num_dev > 0) then
    dev = mod(rank, num_dev)
    call omp_set_default_device(dev)
  else
    dev = omp_get_initial_device()
  end if
#else
  num_dev = 0
  dev = omp_get_initial_device()
#endif

  call MPI_Comm_split(MPI_COMM_WORLD, int(ipxz), int(ipy), comm_y, ierr)
  call MPI_Comm_size(comm_y, comm_size, ierr)
  if (comm_size /= int(npy)) then
    if (rank == 0) print *, "unexpected y communicator size", comm_size
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  nxpp = nx + 1_C_INT
  nxB = nxpp
  nlines_z = 2_C_INT*nz + 1_C_INT
  nlines = nxB*nlines_z
  ny0 = 1_C_INT + ipy*(ny - 1_C_INT)/npy
  nyN = (ipy + 1_C_INT)*(ny - 1_C_INT)/npy
  active_n = nyN - ny0 + 1_C_INT
  arity = npy

  allocate(row_send_counts(arity), row_send_displs(arity), row_recv_counts(arity), row_recv_displs(arity))
  allocate(value_send_counts(arity), value_send_displs(arity), value_recv_counts(arity), value_recv_displs(arity))
  allocate(line_first(arity), line_counts(arity))

  call split_range(ipy, nlines, arity, owned_first, owned_count)
  row_send_offset = 0_C_INT
  row_recv_offset = 0_C_INT
  value_recv_offset = 0_C_INT
  max_line_count = 0_C_INT
  do irank = 0, int(arity) - 1
    call split_range(int(irank, C_INT), nlines, arity, first_line, line_count)
    line_first(irank + 1) = first_line
    line_counts(irank + 1) = line_count
    max_line_count = max(max_line_count, line_count)
    row_send_counts(irank + 1) = int(row_width*line_count)
    row_recv_counts(irank + 1) = int(row_width*owned_count)
    row_send_displs(irank + 1) = int(row_send_offset)
    row_recv_displs(irank + 1) = int(row_recv_offset)
    value_send_counts(irank + 1) = int(value_width*owned_count)
    value_send_displs(irank + 1) = int(value_width*owned_count*irank)
    value_recv_counts(irank + 1) = int(value_width*line_count)
    value_recv_displs(irank + 1) = int(value_recv_offset)
    row_send_offset = row_send_offset + row_width*line_count
    row_recv_offset = row_recv_offset + row_width*owned_count
    value_recv_offset = value_recv_offset + value_width*line_count
  end do
  row_send_elems = row_send_offset
  row_recv_elems = row_recv_offset
  value_send_elems = value_width*arity*owned_count
  value_recv_elems = value_recv_offset
  if (use_alltoall /= 0_C_INT .and. mod(nlines, arity) /= 0_C_INT) then
    if (rank == 0) print *, "EXCHANGE=alltoall requires nlines divisible by NPY"
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
  send_capacity = max(row_send_elems, value_send_elems)
  recv_capacity = max(row_recv_elems, value_recv_elems)

  call allocate_ycomm_buffers(send_capacity, recv_capacity, trim(alloc_mode), ycomm_sendptr, ycomm_recvptr, &
                              ycomm_sendbuf, ycomm_recvbuf)

  allocate(reduced_rows_send(row_width, nlines))
  allocate(reduced_rhs(4_C_INT*npy, nlines))
  allocate(left_interface_values(2, nlines), right_interface_values(2, nlines))
  allocate(s_rows(row_width, nlines, 1))
  allocate(s_values(value_width, nlines, 1))
  allocate(s_recover_basis(4_C_INT*arity, 5, nlines, 1))
  !$omp target enter data map(to: line_first, line_counts, row_recv_displs, value_send_displs, value_recv_displs)
  !$omp target enter data map(alloc: reduced_rows_send, reduced_rhs, left_interface_values, right_interface_values)
  !$omp target enter data map(alloc: s_rows, s_values, s_recover_basis)

  if (rank == 0) then
    print *, "schur_alltoallv_repro"
    print *, "nx=", nx, " ny=", ny, " nz=", nz, " np=", nprocs, " npy=", npy
    print *, "nlines_z=", nlines_z, " nxB=", nxB, " nlines=", nlines
    print *, "row send elems=", row_send_elems, " recv elems=", row_recv_elems
    print *, "value send elems=", value_send_elems, " recv elems=", value_recv_elems
    print *, "row bytes per rank=", int(row_send_elems, C_INT64_T)*16_C_INT64_T
    print *, "value bytes per rank=", int(value_send_elems, C_INT64_T)*16_C_INT64_T
    print *, "iters=", iters, " warmup=", warmup, " wrap_use_device_addr=", wrap_use_device_addr
    print *, "exchange=", trim(exchange_mode), " alloc=", trim(alloc_mode)
  end if
  print *, "Rank", rank, "ipy=", ipy, "active_n=", active_n, "device=", dev, "num_dev=", num_dev

  total_rows = 0.0_C_DOUBLE
  total_values = 0.0_C_DOUBLE

  do iter = 1, warmup + iters
    call roctx_push("ys_schur_pack_leaf_rows")
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(arity, max_line_count, line_counts, line_first, reduced_rows_send, ycomm_sendbuf) &
    !$omp private(dest, local_line, irow, first_line, line_count, global_line, offset)
    do dest = 0, arity - 1
      do local_line = 1, max_line_count
        do irow = 1, row_width
          line_count = line_counts(dest + 1)
          if (local_line <= line_count) then
            first_line = line_first(dest + 1)
            global_line = first_line + local_line - 1_C_INT
            offset = row_width*(first_line - 1_C_INT) + (local_line - 1_C_INT)*row_width + irow
            reduced_rows_send(irow, global_line) = cmplx(dble(irow), dble(global_line), kind=C_DOUBLE)
            ycomm_sendbuf(offset) = reduced_rows_send(irow, global_line)
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctx_pop("ys_schur_pack_leaf_rows")

    call MPI_Barrier(comm_y, ierr)
    if (use_alltoall /= 0_C_INT) then
      call roctx_push("MPI_Alltoall ys_schur_rows")
    else
      call roctx_push("MPI_Alltoallv ys_schur_rows")
    end if
    t0 = MPI_Wtime()
    if (use_alltoall /= 0_C_INT) then
      call alltoall_complex(ycomm_sendbuf, row_send_elems, row_width*owned_count, &
                            ycomm_recvbuf, row_recv_elems, row_width*owned_count, &
                            comm_y, wrap_use_device_addr, ierr)
    else
      call alltoallv_complex(ycomm_sendbuf, row_send_elems, row_send_counts, row_send_displs, &
                             ycomm_recvbuf, row_recv_elems, row_recv_counts, row_recv_displs, &
                             comm_y, wrap_use_device_addr, ierr)
    end if
    t1 = MPI_Wtime()
    if (use_alltoall /= 0_C_INT) then
      call roctx_pop("MPI_Alltoall ys_schur_rows")
    else
      call roctx_pop("MPI_Alltoallv ys_schur_rows")
    end if
    if (iter > warmup) total_rows = total_rows + (t1 - t0)

    call roctx_push("ys_schur_compose_level")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(owned_count, arity, row_recv_displs, ycomm_recvbuf, s_rows, s_recover_basis) &
    !$omp private(local_line, child, k, row0, offset)
    do local_line = 1, owned_count
      do child = 0, arity - 1
        offset = int(row_recv_displs(child + 1), C_INT) + (local_line - 1_C_INT)*row_width
        row0 = 4_C_INT*child
        do k = 1, 4_C_INT
          s_rows(k, local_line, 1) = ycomm_recvbuf(offset + k)
          s_recover_basis(row0 + k, 1, local_line, 1) = ycomm_recvbuf(offset + k)
          s_recover_basis(row0 + k, 2, local_line, 1) = ycomm_recvbuf(offset + 4_C_INT + k)
          s_recover_basis(row0 + k, 3, local_line, 1) = ycomm_recvbuf(offset + 8_C_INT + k)
          s_recover_basis(row0 + k, 4, local_line, 1) = ycomm_recvbuf(offset + 12_C_INT + k)
          s_recover_basis(row0 + k, 5, local_line, 1) = ycomm_recvbuf(offset + 16_C_INT + k)
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctx_pop("ys_schur_compose_level")

    call roctx_push("ys_schur_seed_root_values")
    !$omp target teams distribute parallel do collapse(2) default(none) shared(owned_count, s_values, s_rows) private(local_line, k)
    do local_line = 1, owned_count
      do k = 1, value_width
        if (k <= 4_C_INT) then
          s_values(k, local_line, 1) = s_rows(k, local_line, 1)
        else
          s_values(k, local_line, 1) = (0.0d0, 0.0d0)
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctx_pop("ys_schur_seed_root_values")

    call roctx_push("ys_schur_pack_recovered_values")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(owned_count, arity, value_send_displs, s_recover_basis, s_values, ycomm_sendbuf) &
    !$omp private(local_line, child, row0, k, offset, prev1, prev2, next1, next2, recovered)
    do local_line = 1, owned_count
      prev1 = s_values(5, local_line, 1)
      prev2 = s_values(6, local_line, 1)
      next1 = s_values(7, local_line, 1)
      next2 = s_values(8, local_line, 1)
      do k = 1, 4_C_INT*arity
        recovered(k) = s_recover_basis(k, 1, local_line, 1) + &
                       s_recover_basis(k, 2, local_line, 1)*prev1 + &
                       s_recover_basis(k, 3, local_line, 1)*prev2 + &
                       s_recover_basis(k, 4, local_line, 1)*next1 + &
                       s_recover_basis(k, 5, local_line, 1)*next2
      end do
      do child = 0, arity - 1
        row0 = 4_C_INT*child
        offset = int(value_send_displs(child + 1), C_INT) + (local_line - 1_C_INT)*value_width
        ycomm_sendbuf(offset + 1) = recovered(row0 + 1)
        ycomm_sendbuf(offset + 2) = recovered(row0 + 2)
        ycomm_sendbuf(offset + 3) = recovered(row0 + 3)
        ycomm_sendbuf(offset + 4) = recovered(row0 + 4)
        if (child > 0_C_INT) then
          ycomm_sendbuf(offset + 5) = recovered(row0 - 1)
          ycomm_sendbuf(offset + 6) = recovered(row0)
        else
          ycomm_sendbuf(offset + 5) = prev1
          ycomm_sendbuf(offset + 6) = prev2
        end if
        if (child < arity - 1_C_INT) then
          ycomm_sendbuf(offset + 7) = recovered(row0 + 5)
          ycomm_sendbuf(offset + 8) = recovered(row0 + 6)
        else
          ycomm_sendbuf(offset + 7) = next1
          ycomm_sendbuf(offset + 8) = next2
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctx_pop("ys_schur_pack_recovered_values")

    call MPI_Barrier(comm_y, ierr)
    if (use_alltoall /= 0_C_INT) then
      call roctx_push("MPI_Alltoall ys_schur_values")
    else
      call roctx_push("MPI_Alltoallv ys_schur_values")
    end if
    t0 = MPI_Wtime()
    if (use_alltoall /= 0_C_INT) then
      call alltoall_complex(ycomm_sendbuf, value_send_elems, value_width*owned_count, &
                            ycomm_recvbuf, value_recv_elems, value_width*owned_count, &
                            comm_y, wrap_use_device_addr, ierr)
    else
      call alltoallv_complex(ycomm_sendbuf, value_send_elems, value_send_counts, value_send_displs, &
                             ycomm_recvbuf, value_recv_elems, value_recv_counts, value_recv_displs, &
                             comm_y, wrap_use_device_addr, ierr)
    end if
    t1 = MPI_Wtime()
    if (use_alltoall /= 0_C_INT) then
      call roctx_pop("MPI_Alltoall ys_schur_values")
    else
      call roctx_pop("MPI_Alltoallv ys_schur_values")
    end if
    if (iter > warmup) total_values = total_values + (t1 - t0)

    call roctx_push("ys_schur_unpack_leaf_values")
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(arity, max_line_count, line_counts, line_first, value_recv_displs, ycomm_recvbuf, &
    !$omp& reduced_rhs, left_interface_values, right_interface_values, ipy) &
    !$omp private(src, local_line, k, first_line, line_count, global_line, offset)
    do src = 0, arity - 1
      do local_line = 1, max_line_count
        do k = 1, value_width
          line_count = line_counts(src + 1)
          if (local_line <= line_count) then
            first_line = line_first(src + 1)
            global_line = first_line + local_line - 1_C_INT
            offset = int(value_recv_displs(src + 1), C_INT) + (local_line - 1_C_INT)*value_width
            if (k <= 4_C_INT) reduced_rhs(4_C_INT*ipy + k, global_line) = ycomm_recvbuf(offset + k)
            if (k == 5_C_INT) left_interface_values(1, global_line) = ycomm_recvbuf(offset + k)
            if (k == 6_C_INT) left_interface_values(2, global_line) = ycomm_recvbuf(offset + k)
            if (k == 7_C_INT) right_interface_values(1, global_line) = ycomm_recvbuf(offset + k)
            if (k == 8_C_INT) right_interface_values(2, global_line) = ycomm_recvbuf(offset + k)
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctx_pop("ys_schur_unpack_leaf_values")
  end do

  if (iters > 0_C_INT) then
    print *, "Rank", rank, "avg rows "//trim(exchange_mode)//" s=", total_rows/dble(iters), &
             "avg values "//trim(exchange_mode)//" s=", total_values/dble(iters)
  end if

  !$omp target exit data map(delete: s_rows, s_values, s_recover_basis)
  !$omp target exit data map(delete: reduced_rows_send, reduced_rhs, left_interface_values, right_interface_values)
  !$omp target exit data map(delete: line_first, line_counts, row_recv_displs, value_send_displs, value_recv_displs)
  call free_ycomm_buffers(trim(alloc_mode), ycomm_sendptr, ycomm_recvptr, ycomm_sendbuf, ycomm_recvbuf)
  call MPI_Comm_free(comm_y, ierr)
  call MPI_Finalize(ierr)

contains

  subroutine getenv_int(name, value)
    character(len=*), intent(in) :: name
    integer(C_INT), intent(inout) :: value
    character(len=64) :: raw
    integer :: length, status

    call get_environment_variable(name, raw, length=length, status=status)
    if (status == 0 .and. length > 0) read(raw(1:length), *) value
  end subroutine getenv_int

  subroutine getenv_string(name, value)
    character(len=*), intent(in) :: name
    character(len=*), intent(inout) :: value
    character(len=64) :: raw
    integer :: length, status

    call get_environment_variable(name, raw, length=length, status=status)
    if (status == 0 .and. length > 0) value = raw(1:min(length, len(value)))
  end subroutine getenv_string

  subroutine lowercase(value)
    character(len=*), intent(inout) :: value
    integer :: i, c

    do i = 1, len_trim(value)
      c = iachar(value(i:i))
      if (c >= iachar("A") .and. c <= iachar("Z")) value(i:i) = achar(c + iachar("a") - iachar("A"))
    end do
  end subroutine lowercase

  subroutine split_range(range_rank, nitems, nranks, first_item, item_count)
    integer(C_INT), intent(in) :: range_rank, nitems, nranks
    integer(C_INT), intent(out) :: first_item, item_count
    integer(C_INT) :: base_count, remainder

    base_count = nitems/nranks
    remainder = mod(nitems, nranks)
    item_count = base_count
    if (range_rank < remainder) item_count = item_count + 1_C_INT
    first_item = range_rank*base_count + min(range_rank, remainder) + 1_C_INT
  end subroutine split_range

  subroutine allocate_ycomm_buffers(send_elems, recv_elems, alloc_kind, sendptr, recvptr, sendbuf, recvbuf)
    integer(C_INT), intent(in) :: send_elems, recv_elems
    character(len=*), intent(in) :: alloc_kind
    type(C_PTR), intent(out) :: sendptr, recvptr
    complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: sendbuf(:), recvbuf(:)
#if defined(HAVE_HIP)
    integer(C_INT) :: hip_status
#endif

#if defined(HAVE_HIP)
    if (trim(alloc_kind) == "hip") then
      hip_status = hipMalloc(sendptr, int(max(1_C_INT, send_elems), C_SIZE_T)*int(16, C_SIZE_T))
      if (hip_status /= 0_C_INT) error stop "hipMalloc send failed"
      hip_status = hipMalloc(recvptr, int(max(1_C_INT, recv_elems), C_SIZE_T)*int(16, C_SIZE_T))
      if (hip_status /= 0_C_INT) error stop "hipMalloc recv failed"
    else
      sendptr = omp_target_alloc(int(max(1_C_INT, send_elems), C_SIZE_T)*int(16, C_SIZE_T), omp_get_default_device())
      recvptr = omp_target_alloc(int(max(1_C_INT, recv_elems), C_SIZE_T)*int(16, C_SIZE_T), omp_get_default_device())
    end if
    if (.not. c_associated(sendptr) .or. .not. c_associated(recvptr)) error stop "omp_target_alloc failed"
    call c_f_pointer(sendptr, sendbuf, [max(1_C_INT, send_elems)])
    call c_f_pointer(recvptr, recvbuf, [max(1_C_INT, recv_elems)])
#else
    sendptr = C_NULL_PTR
    recvptr = C_NULL_PTR
    allocate(sendbuf(max(1_C_INT, send_elems)))
    allocate(recvbuf(max(1_C_INT, recv_elems)))
    !$omp target enter data map(alloc: sendbuf, recvbuf)
#endif
  end subroutine allocate_ycomm_buffers

  subroutine free_ycomm_buffers(alloc_kind, sendptr, recvptr, sendbuf, recvbuf)
    character(len=*), intent(in) :: alloc_kind
    type(C_PTR), intent(inout) :: sendptr, recvptr
    complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: sendbuf(:), recvbuf(:)
#if defined(HAVE_HIP)
    integer(C_INT) :: hip_status
#endif

#if defined(HAVE_HIP)
    if (trim(alloc_kind) == "hip") then
      if (c_associated(sendptr)) then
        hip_status = hipFree(sendptr)
        if (hip_status /= 0_C_INT) error stop "hipFree send failed"
      end if
      if (c_associated(recvptr)) then
        hip_status = hipFree(recvptr)
        if (hip_status /= 0_C_INT) error stop "hipFree recv failed"
      end if
    else
      if (c_associated(sendptr)) call omp_target_free(sendptr, omp_get_default_device())
      if (c_associated(recvptr)) call omp_target_free(recvptr, omp_get_default_device())
    end if
    sendptr = C_NULL_PTR
    recvptr = C_NULL_PTR
#else
    !$omp target exit data map(delete: sendbuf, recvbuf)
    deallocate(sendbuf, recvbuf)
#endif
  end subroutine free_ycomm_buffers

  subroutine alltoall_complex(sendbuf, send_elems, send_count, recvbuf, recv_elems, recv_count, &
                              comm, wrap_use_device_addr, ierr_out)
    complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: sendbuf(:), recvbuf(:)
    integer(C_INT), intent(in) :: send_elems, recv_elems, send_count, recv_count, wrap_use_device_addr
    type(MPI_Comm), intent(in) :: comm
    integer, intent(out) :: ierr_out

    if (wrap_use_device_addr /= 0_C_INT) then
      !$omp target data use_device_addr(sendbuf, recvbuf)
      call MPI_Alltoall(sendbuf(1:send_elems), send_count, MPI_DOUBLE_COMPLEX, &
                        recvbuf(1:recv_elems), recv_count, MPI_DOUBLE_COMPLEX, comm, ierr_out)
      !$omp end target data
    else
      call MPI_Alltoall(sendbuf(1:send_elems), send_count, MPI_DOUBLE_COMPLEX, &
                        recvbuf(1:recv_elems), recv_count, MPI_DOUBLE_COMPLEX, comm, ierr_out)
    end if
    if (ierr_out /= MPI_SUCCESS) error stop "MPI_Alltoall failed"
  end subroutine alltoall_complex

  subroutine alltoallv_complex(sendbuf, send_elems, send_counts, send_displs, recvbuf, recv_elems, recv_counts, recv_displs, &
                               comm, wrap_use_device_addr, ierr_out)
    complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: sendbuf(:), recvbuf(:)
    integer(C_INT), intent(in) :: send_elems, recv_elems, wrap_use_device_addr
    integer, intent(in) :: send_counts(:), send_displs(:), recv_counts(:), recv_displs(:)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(out) :: ierr_out

    if (wrap_use_device_addr /= 0_C_INT) then
      !$omp target data use_device_addr(sendbuf, recvbuf)
      call MPI_Alltoallv(sendbuf(1:send_elems), send_counts, send_displs, MPI_DOUBLE_COMPLEX, &
                         recvbuf(1:recv_elems), recv_counts, recv_displs, MPI_DOUBLE_COMPLEX, comm, ierr_out)
      !$omp end target data
    else
      call MPI_Alltoallv(sendbuf(1:send_elems), send_counts, send_displs, MPI_DOUBLE_COMPLEX, &
                         recvbuf(1:recv_elems), recv_counts, recv_displs, MPI_DOUBLE_COMPLEX, comm, ierr_out)
    end if
    if (ierr_out /= MPI_SUCCESS) error stop "MPI_Alltoallv failed"
  end subroutine alltoallv_complex

  subroutine roctx_push(name)
    character(len=*), intent(in) :: name
#if defined(HAVE_HIP)
    character(kind=C_CHAR, len=len_trim(name) + 1) :: cname
    cname = trim(name)//C_NULL_CHAR
    call roctxRangePushA(cname)
#endif
  end subroutine roctx_push

  subroutine roctx_pop(name)
    character(len=*), intent(in) :: name
#if defined(HAVE_HIP)
    call roctxRangePop()
#endif
  end subroutine roctx_pop

end program schur_alltoallv_repro
