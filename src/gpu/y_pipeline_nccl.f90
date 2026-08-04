#include "build_options.h"

module y_pipeline_nccl
  use, intrinsic :: iso_c_binding
  use mpi_f08
  use env_options, only: env_text, lowercase
  implicit none
  private

  integer(c_int), parameter, public :: CHANNEL_COMM_BACKEND_AUTO = 0_c_int
  integer(c_int), parameter, public :: CHANNEL_COMM_BACKEND_MPI = 1_c_int
  integer(c_int), parameter, public :: CHANNEL_COMM_BACKEND_NCCL = 2_c_int

  integer(c_int), save :: channel_comm_override = CHANNEL_COMM_BACKEND_AUTO

  ! Cache of communication contexts keyed by communicator size.  Both the x
  ! transpose and the y pipeline need one, and the autotuner walks several
  ! decompositions within a single run, so contexts are kept per size rather
  ! than rebuilt each time a size comes back into use.
  type, public :: channel_comm_cache
    type(c_ptr), allocatable :: by_size(:)
  end type channel_comm_cache

  public :: channel_comm_cache_reserve, channel_comm_cache_finalize
  public :: channel_comm_use_nccl, channel_comm_available
  public :: channel_comm_backend_from_env, channel_comm_backend_name
  public :: channel_comm_set_backend_override, channel_comm_clear_backend_override
  public :: channel_comm_alltoall_complex
  public :: channel_comm_p2p_ensure, channel_comm_send, channel_comm_recv, channel_comm_sendrecv
  public :: channel_comm_context_reset

#ifdef HAVE_NCCL
  interface
    function channel_nccl_get_unique_id(id_bytes) bind(c, name="channel_nccl_get_unique_id")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: channel_nccl_get_unique_id
      type(c_ptr), value :: id_bytes
    end function channel_nccl_get_unique_id

    function channel_nccl_init(nranks, rank, id_bytes) bind(c, name="channel_nccl_init")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: channel_nccl_init
      integer(c_int), value :: nranks, rank
      type(c_ptr), value :: id_bytes
    end function channel_nccl_init

    function channel_nccl_finalize() bind(c, name="channel_nccl_finalize")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: channel_nccl_finalize
    end function channel_nccl_finalize

    function channel_nccl_sendrecv(sendbuf, send_elems, send_peer, recvbuf, recv_elems, recv_peer) &
      bind(c, name="channel_nccl_sendrecv")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: channel_nccl_sendrecv
      type(c_ptr), value :: sendbuf, recvbuf
      integer(c_size_t), value :: send_elems, recv_elems
      integer(c_int), value :: send_peer, recv_peer
    end function channel_nccl_sendrecv

    function channel_nccl_context_create(nranks, rank, id_bytes, ctx) bind(c, name="channel_nccl_context_create")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: channel_nccl_context_create
      integer(c_int), value :: nranks, rank
      type(c_ptr), value :: id_bytes
      type(c_ptr) :: ctx
    end function channel_nccl_context_create

    function channel_nccl_context_alltoall(ctx, sendbuf, recvbuf, count_elems) bind(c, name="channel_nccl_context_alltoall")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: channel_nccl_context_alltoall
      type(c_ptr), value :: ctx, sendbuf, recvbuf
      integer(c_size_t), value :: count_elems
    end function channel_nccl_context_alltoall

    function channel_nccl_context_sendrecv(ctx, sendbuf, send_elems, send_peer, recvbuf, recv_elems, recv_peer) &
      bind(c, name="channel_nccl_context_sendrecv")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: channel_nccl_context_sendrecv
      type(c_ptr), value :: ctx, sendbuf, recvbuf
      integer(c_size_t), value :: send_elems, recv_elems
      integer(c_int), value :: send_peer, recv_peer
    end function channel_nccl_context_sendrecv

    function channel_nccl_context_destroy(ctx) bind(c, name="channel_nccl_context_destroy")
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: channel_nccl_context_destroy
      type(c_ptr), value :: ctx
    end function channel_nccl_context_destroy
  end interface
#endif

contains
  logical function channel_comm_available()
#ifdef HAVE_NCCL
    channel_comm_available = .true.
#else
    channel_comm_available = .false.
#endif
  end function channel_comm_available

  integer(c_int) function channel_comm_backend_from_env()
    character(len=32) :: value

    channel_comm_backend_from_env = CHANNEL_COMM_BACKEND_AUTO
    if (.not. env_text("CHANNEL_COMM", value)) return
    select case (lowercase(trim(value)))
    case ("auto")
      channel_comm_backend_from_env = CHANNEL_COMM_BACKEND_AUTO
    case ("mpi")
      channel_comm_backend_from_env = CHANNEL_COMM_BACKEND_MPI
    case ("nccl", "rccl")
      channel_comm_backend_from_env = CHANNEL_COMM_BACKEND_NCCL
    case default
      error stop "CHANNEL_COMM must be mpi, auto, or nccl"
    end select
  end function channel_comm_backend_from_env

  character(8) function channel_comm_backend_name(backend)
    integer(c_int), intent(in) :: backend

    select case (backend)
    case (CHANNEL_COMM_BACKEND_MPI)
      channel_comm_backend_name = "mpi"
    case (CHANNEL_COMM_BACKEND_NCCL)
      channel_comm_backend_name = "nccl"
    case (CHANNEL_COMM_BACKEND_AUTO)
      channel_comm_backend_name = "auto"
    case default
      channel_comm_backend_name = "invalid"
    end select
  end function channel_comm_backend_name

  subroutine channel_comm_set_backend_override(backend)
    integer(c_int), intent(in) :: backend

    select case (backend)
    case (CHANNEL_COMM_BACKEND_AUTO, CHANNEL_COMM_BACKEND_MPI)
      channel_comm_override = backend
    case (CHANNEL_COMM_BACKEND_NCCL)
      if (.not. channel_comm_available()) &
        error stop "CHANNEL_COMM=nccl requested, but this build has no NCCL/RCCL support"
      channel_comm_override = backend
    case default
      error stop "invalid CHANNEL_COMM backend override"
    end select
  end subroutine channel_comm_set_backend_override

  subroutine channel_comm_clear_backend_override()
    channel_comm_override = CHANNEL_COMM_BACKEND_AUTO
  end subroutine channel_comm_clear_backend_override

  logical function channel_comm_use_nccl()
    integer(c_int) :: backend

    backend = channel_comm_override
    if (backend == CHANNEL_COMM_BACKEND_AUTO) backend = channel_comm_backend_from_env()
    select case (backend)
    case (CHANNEL_COMM_BACKEND_NCCL)
      if (.not. channel_comm_available()) &
        error stop "CHANNEL_COMM=nccl requested, but this build has no NCCL/RCCL support"
      channel_comm_use_nccl = .true.
    case (CHANNEL_COMM_BACKEND_AUTO, CHANNEL_COMM_BACKEND_MPI)
      channel_comm_use_nccl = .false.
    case default
      error stop "invalid CHANNEL_COMM backend"
    end select
  end function channel_comm_use_nccl

  ! Makes sure cache%by_size(comm_size) exists; the slot itself starts null and
  ! is filled in lazily by whoever first needs a live context for it.
  subroutine channel_comm_cache_reserve(cache, comm_size)
    type(channel_comm_cache), intent(inout) :: cache
    integer(c_int), intent(in) :: comm_size
    type(c_ptr), allocatable :: grown(:)
    integer :: old_upper, new_upper

    if (.not. allocated(cache%by_size)) then
      allocate (cache%by_size(0:max(1_c_int, comm_size)))
      cache%by_size = c_null_ptr
      return
    end if
    if (ubound(cache%by_size, 1) >= comm_size) return

    old_upper = ubound(cache%by_size, 1)
    new_upper = max(int(comm_size), 2*old_upper)
    allocate (grown(0:new_upper))
    grown = c_null_ptr
    grown(0:old_upper) = cache%by_size(0:old_upper)
    call move_alloc(grown, cache%by_size)
  end subroutine channel_comm_cache_reserve

  subroutine channel_comm_cache_finalize(cache)
    type(channel_comm_cache), intent(inout) :: cache
    integer :: i

    if (.not. allocated(cache%by_size)) return
    do i = lbound(cache%by_size, 1), ubound(cache%by_size, 1)
      call channel_comm_context_reset(cache%by_size(i))
    end do
    deallocate (cache%by_size)
  end subroutine channel_comm_cache_finalize

  subroutine channel_comm_p2p_ensure(comm, comm_ctx)
    type(MPI_Comm), intent(in) :: comm
    type(c_ptr), intent(inout) :: comm_ctx
#ifdef HAVE_NCCL
    call channel_comm_ensure_context(comm, comm_ctx)
#else
    error stop "NCCL/RCCL point-to-point backend requested, but this build has no support"
#endif
  end subroutine channel_comm_p2p_ensure

  subroutine channel_comm_context_reset(ctx)
    type(c_ptr), intent(inout) :: ctx
#ifdef HAVE_NCCL
    integer(c_int) :: status

    if (.not. c_associated(ctx)) return
    status = channel_nccl_context_destroy(ctx)
    if (status /= 0_c_int) error stop "channel_nccl_context_destroy failed"
    ctx = c_null_ptr
#else
    ctx = c_null_ptr
#endif
  end subroutine channel_comm_context_reset

  subroutine channel_comm_ensure_context(comm, ctx)
    type(MPI_Comm), intent(in) :: comm
    type(c_ptr), intent(inout) :: ctx
#ifdef HAVE_NCCL
    integer(c_int8_t), target :: id_bytes(128)
    integer(c_int) :: nranks, rank, status
    integer :: ierr

    if (c_associated(ctx)) return

    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, nranks, ierr)
    if (rank == 0_c_int) then
      status = channel_nccl_get_unique_id(c_loc(id_bytes))
      if (status /= 0_c_int) error stop "channel_nccl_get_unique_id failed"
    end if
    call MPI_Bcast(id_bytes, 128, MPI_BYTE, 0, comm, ierr)
    status = channel_nccl_context_create(nranks, rank, c_loc(id_bytes), ctx)
    if (status /= 0_c_int) error stop "channel_nccl_context_create failed"
#else
    error stop "NCCL/RCCL collective backend requested, but this build has no support"
#endif
  end subroutine channel_comm_ensure_context

  subroutine channel_comm_alltoall_complex(sendbuf, recvbuf, count, comm, comm_ctx, request)
    complex(c_double_complex), intent(in), target, contiguous :: sendbuf(:)
    complex(c_double_complex), intent(out), target, contiguous :: recvbuf(:)
    integer(c_int), intent(in), value :: count
    type(MPI_Comm), intent(in) :: comm
    type(c_ptr), intent(inout) :: comm_ctx
    type(MPI_Request), intent(inout), optional :: request
    integer :: ierr

    if (channel_comm_use_nccl()) then
#ifdef HAVE_NCCL
      call channel_comm_ensure_context(comm, comm_ctx)
#ifndef HAVE_HIP
      !$omp target data use_device_addr(sendbuf, recvbuf)
#endif
      call channel_nccl_alltoall(comm_ctx, c_loc(sendbuf(1)), c_loc(recvbuf(1)), count)
#ifndef HAVE_HIP
      !$omp end target data
#endif
      if (present(request)) request = MPI_REQUEST_NULL
#else
      error stop "CHANNEL_COMM=nccl requested, but this build has no NCCL/RCCL support"
#endif
    else
#ifndef HAVE_HIP
      !$omp target data use_device_addr(sendbuf, recvbuf)
#endif
      if (present(request)) then
        call MPI_Ialltoall(sendbuf, int(count), MPI_DOUBLE_COMPLEX, &
                           recvbuf, int(count), MPI_DOUBLE_COMPLEX, comm, request, ierr)
      else
        call MPI_Alltoall(sendbuf, int(count), MPI_DOUBLE_COMPLEX, &
                          recvbuf, int(count), MPI_DOUBLE_COMPLEX, comm, ierr)
      end if
#ifndef HAVE_HIP
      !$omp end target data
#endif
      if (ierr /= MPI_SUCCESS) error stop "channel_comm_alltoall_complex MPI failed"
    end if
  end subroutine channel_comm_alltoall_complex

  subroutine channel_nccl_alltoall(ctx, sendptr, recvptr, count)
    type(c_ptr), intent(in), value :: ctx, sendptr, recvptr
    integer(c_int), intent(in), value :: count
#ifdef HAVE_NCCL
    integer(c_int) :: status

    status = channel_nccl_context_alltoall(ctx, sendptr, recvptr, int(count, c_size_t))
    if (status /= 0_c_int) error stop "channel_nccl_alltoall failed"
#else
    error stop "NCCL/RCCL collective backend requested, but this build has no support"
#endif
  end subroutine channel_nccl_alltoall

  subroutine channel_comm_send(comm_ctx, sendptr, count, peer)
    type(c_ptr), intent(in), value :: comm_ctx, sendptr
    integer(c_int), intent(in), value :: count, peer
#ifdef HAVE_NCCL
    integer(c_int) :: status

    status = channel_nccl_context_sendrecv(comm_ctx, sendptr, int(count, c_size_t), peer, &
                                           c_null_ptr, 0_c_size_t, -1_c_int)
    if (status /= 0_c_int) error stop "channel_comm_send failed"
#else
    error stop "NCCL/RCCL point-to-point backend requested, but this build has no support"
#endif
  end subroutine channel_comm_send

  subroutine channel_comm_recv(comm_ctx, recvptr, count, peer)
    type(c_ptr), intent(in), value :: comm_ctx, recvptr
    integer(c_int), intent(in), value :: count, peer
#ifdef HAVE_NCCL
    integer(c_int) :: status

    status = channel_nccl_context_sendrecv(comm_ctx, c_null_ptr, 0_c_size_t, -1_c_int, &
                                           recvptr, int(count, c_size_t), peer)
    if (status /= 0_c_int) error stop "channel_comm_recv failed"
#else
    error stop "NCCL/RCCL point-to-point backend requested, but this build has no support"
#endif
  end subroutine channel_comm_recv

  subroutine channel_comm_sendrecv(comm_ctx, sendptr, send_count, send_peer, recvptr, recv_count, recv_peer)
    type(c_ptr), intent(in), value :: comm_ctx, sendptr, recvptr
    integer(c_int), intent(in), value :: send_count, send_peer, recv_count, recv_peer
#ifdef HAVE_NCCL
    integer(c_int) :: status

    status = channel_nccl_context_sendrecv(comm_ctx, sendptr, int(send_count, c_size_t), send_peer, &
                                           recvptr, int(recv_count, c_size_t), recv_peer)
    if (status /= 0_c_int) error stop "channel_comm_sendrecv failed"
#else
    error stop "NCCL/RCCL point-to-point backend requested, but this build has no support"
#endif
  end subroutine channel_comm_sendrecv
end module y_pipeline_nccl
