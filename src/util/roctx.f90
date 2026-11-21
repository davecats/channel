#include "build_options.h"

module roctx

  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use mpi_f08, only: MPI_Wtime, MPI_Comm, MPI_COMM_WORLD, MPI_Comm_rank, MPI_Comm_size, MPI_Send, MPI_Recv, &
                     MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE          ! MPI_Wtime: the one clock of all traces
  implicit none

  private
  integer :: n = 0

  public :: roctxpush, roctxpop

  ! ---------------------------------------------------------------------------------------------------
  ! Step trace (CHANNEL_TRACE_STEP=n, CHANNEL_TRACE_DEPTH=d; roctx_step_begin / roctx_step_end in the time-step
  ! loop): while rt_on, every push / pop pair becomes an event (start, end, label id, depth) in fixed host buffers
  ! -- no allocation, and nothing but one logical test per call when off.  Ranges up to depth rt_sync_depth (1 = the
  ! outermost open range) synchronise the device at push and pop, so that the outer phases are exact although GPU
  ! work (hipFFT, nowait regions) may be asynchronous; deeper ranges are host brackets only.  Labels are interned
  ! in rt_labels.
  ! ---------------------------------------------------------------------------------------------------
  integer, parameter :: RT_MAX = 200000, RT_MAXLAB = 512, RT_LABLEN = 64, RT_STACK = 64
  logical, save :: rt_on = .false.
  integer, save :: rt_sync_depth = 2, rt_n = 0, rt_lost = 0, rt_nlab = 0
  double precision, save :: rt_t0(RT_MAX), rt_t1(RT_MAX)
  integer, save :: rt_lab(RT_MAX), rt_depth(RT_MAX)
  character(len=RT_LABLEN), save :: rt_labels(RT_MAXLAB)
  integer, save :: rt_llen(RT_MAXLAB) = 0
  integer, save :: rt_stk_lab(RT_STACK), rt_d = 0
  double precision, save :: rt_stk_t(RT_STACK)
  public :: roctx_step_begin, roctx_step_end
  ! The traced step (0: none; -1: the environment not read yet), its start, and whether it is being recorded.
  integer, save :: rt_step = -1, rt_step_depth = 2
  logical, save :: rt_step_active = .false.
  double precision, save :: rt_step_t0 = 0.0d0
  ! 1 when this rank's last pipelined y-solve ran the mirrored sweep of the twisted pipeline (set by the y solve;
  ! written to the trace header).
  integer, save, public :: rt_mirror = 0

#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
  interface
    integer(c_int) function nvtxrangepush(message) bind(c, name="nvtxRangePushA")
      use iso_c_binding, only: c_char, c_int
      implicit none
      character(c_char) :: message(*)
    end function nvtxrangepush

    integer(c_int) function nvtxrangepop() bind(c, name="nvtxRangePop")
      use iso_c_binding, only: c_int
      implicit none
    end function nvtxrangepop
  end interface
#elif defined(HAVE_HIP)
  interface
    subroutine roctxrangepush(message) bind(c, name="roctxRangePushA")
      use iso_c_binding, only: c_char
      implicit none
      character(c_char) :: message(*)
    end subroutine roctxrangepush

    subroutine roctxrangepop() bind(c, name="roctxRangePop")
      implicit none
    end subroutine roctxrangepop
  end interface
#endif
#if defined(HAVE_HIP)
  interface
    ! REVIEW: linked from the HIP runtime the build already uses (hipFFT, omp offload)
    integer(c_int) function hipdevicesynchronize() bind(c, name="hipDeviceSynchronize")
      use iso_c_binding, only: c_int
      implicit none
    end function hipdevicesynchronize
  end interface
#endif

contains

  double precision function rt_now()
    rt_now = MPI_Wtime()
  end function rt_now

  subroutine rt_device_sync()
#if defined(HAVE_HIP)
    integer(c_int) :: ignored

    ignored = hipDeviceSynchronize()
#else
    !$omp taskwait
#endif
  end subroutine rt_device_sync

  subroutine roctx_trace_start(sync_depth)
    integer, intent(in) :: sync_depth

    rt_sync_depth = sync_depth
    rt_n = 0; rt_lost = 0; rt_d = 0
    call rt_device_sync()
    rt_on = .true.
  end subroutine roctx_trace_start

  subroutine roctx_trace_stop()
    call rt_device_sync()
    rt_on = .false.
  end subroutine roctx_trace_stop

  ! ---- whole time step (CHANNEL_TRACE_STEP=n): called by the time-step loop of the driver around step n.
  ! begin: synchronises the device and starts the recorder.
  ! end:   synchronises, clock offsets to world rank 0, writes ttrace_<world rank>.txt with
  !   # header (rank, ipy, ny0, nyN, offset, rtt, lost events, mirror, the label table `# label <id> <text>`, the note
  !   that the trace is synchronising at roctx depth <= CHANNEL_TRACE_DEPTH), then  t0_us t1_us code a b  lines:
  !   code 300 = a roctx range (a = label id, b = depth, 1 = the outermost open range).  Times on world rank 0's
  !   clock, relative to its start of the step.  Collective over MPI_COMM_WORLD.
  subroutine roctx_step_begin(step)
    use env_options, only: env_int
    integer, intent(in) :: step
    integer(c_int) :: v

    if (rt_step < 0) then
      rt_step = 0
      if (env_int("CHANNEL_TRACE_STEP", v)) rt_step = max(0, int(v))
      if (env_int("CHANNEL_TRACE_DEPTH", v)) rt_step_depth = max(0, int(v))
    end if
    if (rt_step == 0 .or. step /= rt_step) return
    call roctx_trace_start(rt_step_depth)
    rt_step_active = .true.
    rt_step_t0 = rt_now()
  end subroutine roctx_step_begin

  ! ipy, npy, ny0, nyN: this rank's place in the y decomposition, for the header.
  subroutine roctx_step_end(ipy, npy, ny0, nyN)
    integer, intent(in) :: ipy, npy, ny0, nyN
    integer :: wr, wn, ierr, k, u, ios
    double precision :: buf(3)
    character(len=64) :: fn

    if (.not. rt_step_active) return
    call roctx_trace_stop()
    rt_step_active = .false.
    call MPI_Comm_rank(MPI_COMM_WORLD, wr, ierr); call MPI_Comm_size(MPI_COMM_WORLD, wn, ierr)
    call rt_clock_sync(MPI_COMM_WORLD, wr, wn, rt_step_t0, buf)
    write (fn, '(a,i0,a)') "ttrace_", wr, ".txt"
    open (newunit=u, file=trim(fn), status='replace', action='write', iostat=ios)
    if (ios == 0) then
      write (u, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,f12.3,a,f12.3,a,i0,a,i0)') "# rank ", wr, " nranks ", wn, " ipy ", ipy, &
        " npy ", npy, " ny0 ", ny0, " nyN ", nyN, " offset_us ", buf(1)*1.0d6, " rtt_us ", buf(2)*1.0d6, " lost ", rt_lost, &
        " mirror ", rt_mirror
      write (u, '(a,i0)') "# synchronizing trace: the device is synchronised at every roctx push/pop of depth <= ", rt_step_depth
      write (u, '(a)') "# t0_us t1_us code a b   (300: roctx range, a = label id, b = depth)"
      do k = 1, rt_nlab
        write (u, '(a,i0,1x,a)') "# label ", k, trim(rt_labels(k))
      end do
      do k = 1, rt_n
        write (u, '(f13.3,1x,f13.3,1x,i0,1x,i0,1x,i0)') (rt_t0(k) - buf(1) - buf(3))*1.0d6, (rt_t1(k) - buf(1) - buf(3))*1.0d6, &
          300, rt_lab(k), rt_depth(k)
      end do
      close (u)
    end if
  end subroutine roctx_step_end

  ! Clock offset of this rank to rank 0 of `comm` (ping-pong, ROUNDS rounds one rank after the other, the round of
  ! minimal round-trip time; offset = t_remote - (t_send + t_recv)/2).  buf = (offset, rtt, ref) where ref is rank 0's
  ! `ref` as sent to all.  Collective over comm.
  subroutine rt_clock_sync(comm, myrank, nrank, ref, buf)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: myrank, nrank
    double precision, intent(in) :: ref
    double precision, intent(out) :: buf(3)
    integer, parameter :: TAG = 9300, ROUNDS = 30
    integer :: r, k, ierr
    double precision :: ta, tb, tr, best, off

    off = 0.0d0; best = 0.0d0; buf = [0.0d0, 0.0d0, ref]
    if (nrank <= 1) return
    if (myrank == 0) then
      do r = 1, nrank - 1
        best = huge(1.0d0)
        do k = 1, ROUNDS
          ta = MPI_Wtime()
          call MPI_Send(ta, 1, MPI_DOUBLE_PRECISION, r, TAG, comm, ierr)
          call MPI_Recv(tr, 1, MPI_DOUBLE_PRECISION, r, TAG, comm, MPI_STATUS_IGNORE, ierr)
          tb = MPI_Wtime()
          if (tb - ta < best) then
            best = tb - ta; off = tr - 0.5d0*(ta + tb)
          end if
        end do
        buf = [off, best, ref]
        call MPI_Send(buf, 3, MPI_DOUBLE_PRECISION, r, TAG, comm, ierr)
      end do
      buf = [0.0d0, 0.0d0, ref]
    else
      do k = 1, ROUNDS
        call MPI_Recv(ta, 1, MPI_DOUBLE_PRECISION, 0, TAG, comm, MPI_STATUS_IGNORE, ierr)
        tr = MPI_Wtime()
        call MPI_Send(tr, 1, MPI_DOUBLE_PRECISION, 0, TAG, comm, ierr)
      end do
      call MPI_Recv(buf, 3, MPI_DOUBLE_PRECISION, 0, TAG, comm, MPI_STATUS_IGNORE, ierr)
    end if
  end subroutine rt_clock_sync

  integer function rt_intern(name)
    character(len=*), intent(in) :: name
    integer :: i, l

    l = len_trim(name)
    do i = 1, rt_nlab
      if (rt_llen(i) == l) then
        if (rt_labels(i)(1:min(l, RT_LABLEN)) == name(1:min(l, RT_LABLEN))) then
          rt_intern = i
          return
        end if
      end if
    end do
    if (rt_nlab >= RT_MAXLAB) then
      rt_intern = RT_MAXLAB                  ! table full: the last label stands for all further ones
      return
    end if
    rt_nlab = rt_nlab + 1
    rt_labels(rt_nlab) = name(1:min(l, RT_LABLEN)); rt_llen(rt_nlab) = l
    rt_intern = rt_nlab
  end function rt_intern

  subroutine roctxPush(name)
    character(len=*), intent(in) :: name
#if (defined(HAVE_CUDA) && defined(HAVE_NVTX)) || defined(HAVE_HIP)
    ! Only built on profiled builds: the null-terminated copy is an automatic
    ! character temporary, and these calls sit in the y-solve inner loop.
    character(kind=c_char, len=len_trim(name) + 1) :: cname
#endif
#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
    integer(c_int) :: ignored
#endif

    n = n + 1
    if (rt_on) then
      if (rt_d < RT_STACK) then
        rt_d = rt_d + 1
        if (rt_d <= rt_sync_depth) call rt_device_sync()
        rt_stk_lab(rt_d) = rt_intern(name); rt_stk_t(rt_d) = rt_now()
      end if
    end if
#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
    cname = trim(name)//c_null_char
    ignored = nvtxRangePush(cname)
#elif defined(HAVE_HIP)
    cname = trim(name)//c_null_char
    call roctxRangePush(cname)
#endif
  end subroutine roctxPush

  subroutine roctxPop(name)
    character(len=*), intent(in) :: name
#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
    integer(c_int) :: ignored
#endif
    double precision :: t1

    n = n - 1
    if (n < 0) then
      print *, "invalid pop for: ", name
      return
    end if
    if (rt_on .and. rt_d > 0) then
      if (rt_d <= rt_sync_depth) call rt_device_sync()
      t1 = rt_now()
      if (rt_n < RT_MAX) then
        rt_n = rt_n + 1
        rt_t0(rt_n) = rt_stk_t(rt_d); rt_t1(rt_n) = t1; rt_lab(rt_n) = rt_stk_lab(rt_d); rt_depth(rt_n) = rt_d
      else
        rt_lost = rt_lost + 1
      end if
      rt_d = rt_d - 1
    end if
#if defined(HAVE_CUDA) && defined(HAVE_NVTX)
    ignored = nvtxRangePop()
#elif defined(HAVE_HIP)
    call roctxRangePop()
#endif
  end subroutine roctxPop

end module roctx
