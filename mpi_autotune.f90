#include "header.h"

module mpi_autotune
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use mpi_transpose, only: init_MPI, free_MPI, nxB, nzB, ny0, nyN, fft_transpose_is_local, &
                           repack_zTOx_local, pack_zTOx, alltoall, unpack_zTOx, &
                           repack_xTOz_local, pack_xTOz, unpack_xTOz, sendbuf, recvbuf
#ifdef HAVE_CUDA
  use ffts, only: init_cufft, free_fft, VVdz, VVdx, rVVdx
#elif defined(HAVE_HIP)
  use ffts, only: init_hipfft, free_fft, VVdz, VVdx, rVVdx
#elif defined(HAVE_FFTW)
  use ffts, only: init_fft, free_fft
#endif
  use y_line_solvers, only: ys_prepare_assembled_workspace, ys_release_workspace, &
                            ys_solve_endpoint_schur, ys_solve_pipelined_lu, ys_gpsv_ds, ys_gpsv_dl, &
                            ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x
  use y_schur_solver, only: ys_schur_default_pass_counts, YS_SCHUR_EXCHANGE_AUTO, &
                            YS_SCHUR_EXCHANGE_ALLTOALL, YS_SCHUR_EXCHANGE_ALLGATHER
  use byte_workspace, only: workspace_finalize
  use mpi_f08
  implicit none
  private
  integer(C_INT), parameter :: MAXP = 16_C_INT
  integer(C_INT), parameter :: ARITY(5) = [2_C_INT, 3_C_INT, 4_C_INT, 6_C_INT, 8_C_INT]
  integer(C_INT), parameter, public :: Y_SOLVER_SCHUR = 1_C_INT
  integer(C_INT), parameter, public :: Y_SOLVER_PIPELINED_LU = 2_C_INT
  real(C_DOUBLE), parameter :: RK_SUBSTEPS = 3.0_C_DOUBLE
  real(C_DOUBLE), parameter :: Y_SOLVE_CHECK_TOL = 1.0e-6_C_DOUBLE
  integer(C_INT), save, public :: mpi_autotune_selected_y_solver = Y_SOLVER_SCHUR
  integer(C_INT), save, public :: mpi_autotune_selected_y_batches = 0_C_INT
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(autotune_exact_value)
#endif
  public :: configure_mpi_decomposition, mpi_autotune_has_pass_sequence

contains

  subroutine configure_mpi_decomposition(nxpp, nxd, nzd, nz, ny, nphi, overlapping, requested_npy, &
                                         npy_out, npxz_out, passes, exchange)
    integer(C_INT), intent(in) :: nxpp, nxd, nzd, nz, ny, nphi, requested_npy
    logical, intent(in) :: overlapping
    integer(C_INT), intent(out) :: npy_out, npxz_out, exchange
    integer(C_INT), allocatable, intent(out) :: passes(:)
    character(256) :: text, tmp
    integer :: ierr, rank, nranks, status, length, io, i, mode, nforced
    integer(C_INT) :: env_npy, env_npxz, forced(MAXP), best_pass(MAXP), best_npass, best_npxz, best_npy, best_exchange
    integer(C_INT) :: best_y_solver, best_y_batches
    integer(C_INT), allocatable :: node(:)
    logical :: has_npy, has_npxz, has_passes, has_exchange, manual, applied, found, ran_scan
    logical :: has_y_solver
    real(C_DOUBLE) :: best_score

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

    npy_out = max(1_C_INT, requested_npy); npxz_out = 1_C_INT
    call get_environment_variable("CHANNEL_NPY", text, length, status); has_npy = (status == 0); env_npy = 0_C_INT
    if (has_npy) read (text(:length), *, iostat=io) env_npy
    call get_environment_variable("CHANNEL_NPXZ", text, length, status); has_npxz = (status == 0); env_npxz = 0_C_INT
    if (has_npxz) read (text(:length), *, iostat=io) env_npxz
    call get_environment_variable("CHANNEL_Y_SOLVER", text, length, status); has_y_solver = (status == 0)
    mpi_autotune_selected_y_solver = Y_SOLVER_SCHUR
    mpi_autotune_selected_y_batches = 0_C_INT
    if (has_y_solver) then
      select case (adjustl(trim(text(:length))))
      case ("auto", "AUTO", "default", "DEFAULT")
        has_y_solver = .false.
      case ("schur", "SCHUR")
        mpi_autotune_selected_y_solver = Y_SOLVER_SCHUR
      case ("pipelined_lu", "PIPELINED_LU", "pipelined-lu", "PIPELINED-LU")
        mpi_autotune_selected_y_solver = Y_SOLVER_PIPELINED_LU
      case default
        call abort_msg(rank, "invalid CHANNEL_Y_SOLVER")
      end select
    end if
    if (has_npy .and. env_npy < 1_C_INT) call abort_msg(rank, "CHANNEL_NPY must be >= 1")
    if (has_npxz .and. env_npxz < 1_C_INT) call abort_msg(rank, "CHANNEL_NPXZ must be >= 1")
    if (has_npy .and. has_npxz) then
      if (env_npy*env_npxz /= int(nranks, C_INT)) call abort_msg(rank, "CHANNEL_NPY * CHANNEL_NPXZ must equal rank count")
      npy_out = env_npy; npxz_out = env_npxz
    else if (has_npy) then
      if (mod(int(nranks, C_INT), env_npy) /= 0_C_INT) call abort_msg(rank, "rank count must be divisible by CHANNEL_NPY")
      npy_out = env_npy; npxz_out = int(nranks, C_INT)/env_npy
    else if (has_npxz) then
      if (mod(int(nranks, C_INT), env_npxz) /= 0_C_INT) call abort_msg(rank, "rank count must be divisible by CHANNEL_NPXZ")
      npxz_out = env_npxz; npy_out = int(nranks, C_INT)/env_npxz
    else
      if (mod(int(nranks, C_INT), npy_out) /= 0_C_INT) call abort_msg(rank, "requested npy must divide rank count")
      npxz_out = int(nranks, C_INT)/npy_out
    end if

    forced = 0_C_INT; nforced = 0
    call get_environment_variable("CHANNEL_Y_SCHUR_PASSES", text, length, status); has_passes = (status == 0)
    if (has_passes) then
      tmp = text
      do i = 1, length
        if (index(",;:xX", tmp(i:i)) > 0) tmp(i:i) = " "
      end do
      read (tmp, *, iostat=io) forced
      do i = 1, MAXP
        if (forced(i) == 0_C_INT) exit
        nforced = nforced + 1
      end do
    end if
    call get_environment_variable("CHANNEL_Y_SCHUR_GLOBAL_EXCHANGE", text, length, status)
    if (status /= 0) call get_environment_variable("CHANNEL_Y_SCHUR_EXCHANGE", text, length, status)
    has_exchange = (status == 0); exchange = YS_SCHUR_EXCHANGE_AUTO
    if (has_exchange) then
      select case (adjustl(trim(text(:length))))
      case ("auto", "AUTO", "default", "DEFAULT"); exchange = YS_SCHUR_EXCHANGE_AUTO
      case ("alltoall", "ALLTOALL", "alltoallv", "ALLTOALLV"); exchange = YS_SCHUR_EXCHANGE_ALLTOALL
      case ("allgather", "ALLGATHER", "allgatherv", "ALLGATHERV"); exchange = YS_SCHUR_EXCHANGE_ALLGATHER
      case default; call abort_msg(rank, "invalid CHANNEL_Y_SCHUR_EXCHANGE")
      end select
    end if

    call get_environment_variable("CHANNEL_MPI_AUTOTUNE", text, length, status); mode = 1
    if (status == 0) then
      select case (adjustl(trim(text(:length))))
      case ("0", "false", "FALSE", "off", "OFF", "no", "NO"); mode = 0
      case ("report", "REPORT"); mode = 2
      end select
    end if
    manual = has_npy .or. has_npxz .or. has_passes .or. has_exchange .or. has_y_solver
    found = .false.; applied = .false.; ran_scan = .false.; best_score = huge(0.0_C_DOUBLE)
    if ((mode == 1 .and. .not. manual) .or. mode == 2) then
      ran_scan = .true.
      call node_ids(node)
      call scan(int(nranks, C_INT), nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, best_score, &
                best_npxz, best_npy, best_pass, best_npass, best_exchange, best_y_solver, best_y_batches, found)
      if (rank == 0 .and. found) &
        call print_config("MPI autotune recommendation", best_npxz, best_npy, best_pass, best_npass, best_exchange)
      if (found .and. mode == 1 .and. .not. manual) then
        npxz_out = best_npxz; npy_out = best_npy; exchange = best_exchange
        mpi_autotune_selected_y_solver = best_y_solver
        mpi_autotune_selected_y_batches = best_y_batches
        applied = .true.
      else if (.not. found .and. rank == 0) then
        print *, "Warning: MPI autotune found no valid candidates; keeping configured decomposition."
      end if
      deallocate (node)
    end if
    if (ran_scan) then
      call workspace_finalize()
    end if

    if (applied) then
      allocate (passes(best_npass)); if (best_npass > 0) passes = best_pass(1:best_npass)
    else if (has_passes) then
      allocate (passes(nforced)); if (nforced > 0) passes = forced(1:nforced)
      if (.not. mpi_autotune_has_pass_sequence(npy_out, passes)) &
        call abort_msg(rank, "CHANNEL_Y_SCHUR_PASSES product/arity is invalid")
    else
      call ys_schur_default_pass_counts(npy_out, passes)
    end if
    if (rank == 0 .and. applied) &
      call print_config("MPI autotune selected", npxz_out, npy_out, passes, int(size(passes), C_INT), exchange)
    if (rank == 0 .and. .not. applied) &
      call print_config("MPI decomposition configured", npxz_out, npy_out, passes, int(size(passes), C_INT), exchange)
  end subroutine configure_mpi_decomposition

  subroutine scan(nranks, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, best_score, &
                  best_npxz, best_npy, best_pass, best_npass, best_exchange, best_y_solver, best_y_batches, found)
    integer(C_INT), intent(in) :: nranks, nxpp, nxd, nzd, nz, ny, nphi, node(0:)
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(out) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange, best_y_solver, best_y_batches
    logical, intent(inout) :: found
    integer(C_INT) :: npxz, path(MAXP)
    path = 1_C_INT
    do npxz = 1_C_INT, nranks
      if (mod(nranks, npxz) == 0_C_INT) &
        call gen(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
                 nranks/npxz, path, 0_C_INT, best_score, best_npxz, best_npy, best_pass, &
                 best_npass, best_exchange, best_y_solver, best_y_batches, found)
    end do
  end subroutine scan

  recursive subroutine gen(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
                           remaining, path, npass, best_score, best_npxz, best_npy, best_pass, &
                           best_npass, best_exchange, best_y_solver, best_y_batches, found)
    integer(C_INT), intent(in) :: nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, node(0:), remaining, path(MAXP), npass
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(inout) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange
    integer(C_INT), intent(inout) :: best_y_solver, best_y_batches
    logical, intent(inout) :: found
    integer(C_INT) :: i, next_path(MAXP)
    if (remaining == 1_C_INT) then
      call try_candidate(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, path, &
                         npass, YS_SCHUR_EXCHANGE_ALLTOALL, best_score, best_npxz, best_npy, &
                         best_pass, best_npass, best_exchange, best_y_solver, best_y_batches, found)
      if (npass > 0_C_INT) then
        if (path(npass) < 4_C_INT) &
          call try_candidate(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, path, &
                             npass, YS_SCHUR_EXCHANGE_ALLGATHER, best_score, best_npxz, best_npy, &
                             best_pass, best_npass, best_exchange, best_y_solver, best_y_batches, found)
      end if
      return
    end if
    if (npass >= MAXP) return
    do i = 1_C_INT, int(size(ARITY), C_INT)
      if (mod(remaining, ARITY(i)) /= 0_C_INT) cycle
      next_path = path; next_path(npass + 1_C_INT) = ARITY(i)
      call gen(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
               remaining/ARITY(i), next_path, npass + 1_C_INT, best_score, best_npxz, &
               best_npy, best_pass, best_npass, best_exchange, best_y_solver, best_y_batches, found)
    end do
  end subroutine gen

  subroutine try_candidate(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, path, &
                           npass, exchange, best_score, best_npxz, best_npy, best_pass, &
                           best_npass, best_exchange, best_y_solver, best_y_batches, found)
    integer(C_INT), intent(in) :: nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, node(0:), path(MAXP), npass, exchange
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(inout) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange
    integer(C_INT), intent(inout) :: best_y_solver, best_y_batches
    logical, intent(inout) :: found
    integer :: ierr, rank
    real(C_DOUBLE) :: xz_forward_cost, xz_back_cost, y_cost, score
    real(C_DOUBLE) :: y_error
    integer(C_INT) :: ib, batches, nlines, batch_candidates(5)
    logical :: y_ok
    if (.not. valid(nranks, npxz, nxpp, nzd, nz, ny, node, path, npass)) return
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call time_xz_sweep(nxpp, nxd, nzd, nz, ny, nphi, overlapping, nranks/npxz, xz_forward_cost, xz_back_cost)

    call score_y_backend(Y_SOLVER_SCHUR, 0_C_INT)
    if (nranks/npxz > 1_C_INT) then
      nlines = (nxpp/npxz)*(2_C_INT*nz + 1_C_INT)
      batch_candidates = [1_C_INT, 2_C_INT, 4_C_INT, 8_C_INT, 16_C_INT]
      do ib = 1_C_INT, int(size(batch_candidates), C_INT)
        batches = batch_candidates(ib)
        if (batches <= nlines) call score_y_backend(Y_SOLVER_PIPELINED_LU, batches)
      end do
    end if

  contains
    subroutine score_y_backend(y_solver, y_batches)
      integer(C_INT), intent(in) :: y_solver, y_batches

      y_cost = time_y(nxpp, nzd, nz, ny, nphi, overlapping, nranks/npxz, path, npass, exchange, &
                      y_solver, y_batches, y_ok, y_error)
      if (.not. y_ok) then
        if (rank == 0) write (*, '(*(g0,1x))') "MPI autotune rejected: npxz=", npxz, &
          "npy=", nranks/npxz, "passes=", trim(pass_string(path, npass)), &
          "exchange=", trim(exchange_string(exchange)), "y_solver=", trim(y_solver_string(y_solver)), &
          "y_batches=", y_batches, "y_correctness_error=", y_error
        return
      end if
      score = RK_SUBSTEPS*(xz_forward_cost*real(3_C_INT + nphi, C_DOUBLE) + &
                           xz_back_cost*real(6_C_INT + 3_C_INT*nphi, C_DOUBLE) + &
                           y_cost*real(3_C_INT + nphi, C_DOUBLE))
      if (rank == 0) write (*, '(*(g0,1x))') "MPI autotune tested: npxz=", npxz, &
        "npy=", nranks/npxz, "passes=", trim(pass_string(path, npass)), &
        "exchange=", trim(exchange_string(exchange)), "xz_forward_ms=", 1d3*xz_forward_cost, &
        "xz_back_ms=", 1d3*xz_back_cost, "y_solver=", trim(y_solver_string(y_solver)), &
        "y_batches=", y_batches, "y_ms=", 1d3*y_cost, "y_correctness_error=", y_error, &
        "score_timestep_ms=", 1d3*score
      if (.not. found .or. score < best_score) then
        found = .true.; best_score = score; best_npxz = npxz; best_npy = nranks/npxz
        best_pass = path; best_npass = npass; best_exchange = exchange
        best_y_solver = y_solver; best_y_batches = y_batches
      end if
    end subroutine score_y_backend
  end subroutine try_candidate

  logical function valid(nranks, npxz, nxpp, nzd, nz, ny, node, path, npass)
    integer(C_INT), intent(in) :: nranks, npxz, nxpp, nzd, nz, ny, node(0:), path(MAXP), npass
    integer(C_INT) :: npy, ipy, ipxz, prev, level, y0, yN
    valid = .false.; npy = nranks/npxz
    if (mod(nxpp, npxz) /= 0_C_INT .or. mod(nzd, npxz) /= 0_C_INT) return
    if (npy > 1_C_INT) then
      do ipy = 0_C_INT, npy - 1_C_INT
        y0 = 1_C_INT + ipy*(ny - 1_C_INT)/npy
        yN = (ipy + 1_C_INT)*(ny - 1_C_INT)/npy
        if (yN - y0 + 1_C_INT < 4_C_INT) return
      end do
    end if
    if (npxz > 1_C_INT) then
      do ipy = 0_C_INT, npy - 1_C_INT
        if (.not. local(node, ipy*npxz, npxz, 1_C_INT)) return
      end do
    end if
    do ipxz = 0_C_INT, npxz - 1_C_INT; if (.not. clean(node, ipxz, npy, npxz)) return; end do
    prev = (nxpp/npxz)*(2_C_INT*nz + 1_C_INT)
    do level = 1_C_INT, npass
      if (mod(prev, path(level)) /= 0_C_INT .or. .not. clean_schur_level(npxz, npy, node, path, level)) return
      prev = prev/path(level)
    end do
    valid = .true.
  end function valid

  logical function clean_schur_level(npxz, npy, node, path, level)
    integer(C_INT), intent(in) :: npxz, npy, node(0:), path(MAXP), level
    integer(C_INT) :: span, child_span, parent, pos, ipxz, i
    span = 1_C_INT
    do i = 1_C_INT, level - 1_C_INT; span = span*path(i); end do
    child_span = span; span = span*path(level); clean_schur_level = .true.
    do ipxz = 0_C_INT, npxz - 1_C_INT; do parent = 0_C_INT, npy/span - 1_C_INT; do pos = 0_C_INT, child_span - 1_C_INT
        if (.not. clean_schur_group(npxz, node, ipxz, parent*span + pos, child_span, path(level))) then
          clean_schur_level = .false.; return
        end if
      end do; end do; end do
  end function clean_schur_level

  logical function clean_schur_group(npxz, node, ipxz, first_ipy, child_span, arity)
    integer(C_INT), intent(in) :: npxz, node(0:), ipxz, first_ipy, child_span, arity
    integer(C_INT) :: i, j, group_node(8)
    do i = 1_C_INT, arity; group_node(i) = node((first_ipy + (i - 1_C_INT)*child_span)*npxz + ipxz); end do
    clean_schur_group = .true.
    do i = 2_C_INT, arity; if (group_node(i) /= group_node(1)) clean_schur_group = .false.; end do
    if (clean_schur_group) return
    clean_schur_group = .true.
    do i = 1_C_INT, arity
      do j = i + 1_C_INT, arity
        if (group_node(i) == group_node(j)) clean_schur_group = .false.
      end do
    end do
  end function clean_schur_group

  logical function local(node, first, count, stride)
    integer(C_INT), intent(in) :: node(0:), first, count, stride
    integer(C_INT) :: i
    local = .true.
    do i = 1_C_INT, count - 1_C_INT; if (node(first + i*stride) /= node(first)) local = .false.; end do
  end function local

  logical function clean(node, first, count, stride)
    integer(C_INT), intent(in) :: node(0:), first, count, stride
    integer(C_INT) :: i, j
    clean = .true.; if (local(node, first, count, stride)) return
    do i = 0_C_INT, count - 1_C_INT; do j = i + 1_C_INT, count - 1_C_INT
        if (node(first + i*stride) == node(first + j*stride)) clean = .false.
      end do; end do
  end function clean

  subroutine time_xz_sweep(nxpp, nxd, nzd, nz, ny, nphi, overlapping, npy, forward_cost, back_cost)
    integer(C_INT), intent(in) :: nxpp, nxd, nzd, nz, ny, nphi, npy
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(out) :: forward_cost, back_cost
#ifdef HAVE_FFTW
    complex(C_DOUBLE_COMPLEX), pointer :: VVdz(:, :, :, :), VVdx(:, :, :, :)
    real(C_DOUBLE), pointer :: rVVdx(:, :, :, :)
#endif
    type(MPI_Request) :: request
    type(MPI_Status) :: status
    integer :: ierr, iter, nrepeat
    real(C_DOUBLE) :: t0, forward_elapsed, back_elapsed, global_elapsed

    call init_MPI(nxpp, nz, ny, nzd, nphi, overlapping, npy, .true.)
#ifdef HAVE_CUDA
    call init_cufft(nxd, nxB, nzd, nzB, nphi, overlapping)
#elif defined(HAVE_HIP)
    call init_hipfft(nxd, nxB, nzd, nzB, nphi, overlapping)
#elif defined(HAVE_FFTW)
    call init_fft(VVdz, VVdx, rVVdx, nxd, nxB, nzd, nzB, nphi, overlapping)
#endif

    nrepeat = tune_repeats()
    forward_elapsed = 0.0_C_DOUBLE
    back_elapsed = 0.0_C_DOUBLE
    do iter = 0, nrepeat
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      t0 = MPI_Wtime()

      if (fft_transpose_is_local) then
        call repack_zTOx_local(VVdz(:, :, :, 1), VVdx(:, :, :, 1), ny)
      else
        call pack_zTOx(VVdz(:, :, :, 1), sendbuf(:, 1), ny)
        call alltoall(sendbuf(:, 1), recvbuf(:, 1), request, "zTOx autotune_xz_sweep")
        call MPI_Wait(request, status, ierr)
        call unpack_zTOx(recvbuf(:, 1), VVdx(:, :, :, 1), ny)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (iter > 0) forward_elapsed = forward_elapsed + MPI_Wtime() - t0

      t0 = MPI_Wtime()
      if (fft_transpose_is_local) then
        call repack_xTOz_local(VVdx(:, :, :, 1), VVdz(:, :, :, 1), ny)
      else
        call pack_xTOz(VVdx(:, :, :, 1), sendbuf(:, 1), ny)
        call alltoall(sendbuf(:, 1), recvbuf(:, 1), request, "xTOz autotune_xz_sweep")
        call MPI_Wait(request, status, ierr)
        call unpack_xTOz(recvbuf(:, 1), VVdz(:, :, :, 1), ny)
      end if

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (iter > 0) back_elapsed = back_elapsed + MPI_Wtime() - t0
    end do

    forward_elapsed = forward_elapsed/real(max(1, nrepeat), C_DOUBLE)
    back_elapsed = back_elapsed/real(max(1, nrepeat), C_DOUBLE)
    call MPI_Allreduce(forward_elapsed, global_elapsed, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    forward_cost = global_elapsed
    call MPI_Allreduce(back_elapsed, global_elapsed, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    back_cost = global_elapsed

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    call free_fft()
#elif defined(HAVE_FFTW)
    call free_fft(VVdz, VVdx, rVVdx)
#endif
    call free_MPI()
  end subroutine time_xz_sweep

  real(C_DOUBLE) function time_y(nxpp, nzd, nz, ny, nphi, overlapping, npy, path, npass, exchange, &
                                 y_solver, y_batches, ok, max_error)
    integer(C_INT), intent(in) :: nxpp, nzd, nz, ny, nphi, npy, path(MAXP), npass, exchange
    integer(C_INT), intent(in) :: y_solver, y_batches
    logical, intent(in) :: overlapping
    logical, intent(out) :: ok
    real(C_DOUBLE), intent(out) :: max_error
    integer(C_INT), allocatable :: passes(:)
    integer(C_INT) :: nlines

    call init_MPI(nxpp, nz, ny, nzd, nphi, overlapping, npy, .true.)
    allocate (passes(npass))
    if (npass > 0_C_INT) passes = path(1:npass)
    nlines = nxB*(2_C_INT*nz + 1_C_INT)
    call time_y_endpoint_solve(ny, nz, ny0, nyN, 1_C_INT, nlines, passes, exchange, &
                               tune_repeats(), y_solver, y_batches, time_y, ok, max_error)
    deallocate (passes)
    call free_MPI()
  end function time_y

  subroutine time_y_endpoint_solve(ny, nz, row_start, row_end, line_start, nlines, &
                                   passes, exchange, repeats, y_solver, y_batches, elapsed, ok, max_error)
    integer(C_INT), intent(in) :: ny, nz, row_start, row_end, line_start, nlines
    integer(C_INT), intent(in) :: passes(:), exchange
    integer, intent(in) :: repeats
    integer(C_INT), intent(in) :: y_solver, y_batches
    real(C_DOUBLE), intent(out) :: elapsed
    logical, intent(out) :: ok
    real(C_DOUBLE), intent(out) :: max_error
    complex(C_DOUBLE_COMPLEX), allocatable :: dst(:, :, :)
    integer(C_INT) :: active_n, nx_count
    integer(C_INT64_T) :: t0, t1, rate
    integer :: ierr, iter, bad_local, bad_global
    real(C_DOUBLE) :: local_elapsed, local_error, global_error

    call ys_prepare_assembled_workspace(ny, nz, row_start, row_end, line_start, nlines, passes, exchange)
    active_n = row_end - row_start + 1_C_INT
    nx_count = nlines/(2_C_INT*nz + 1_C_INT)
    allocate (dst(active_n + 4_C_INT, 2_C_INT*nz + 1_C_INT, nx_count))
    !$omp target enter data map(alloc: dst)

    call system_clock(count_rate=rate)
    local_elapsed = 0.0_C_DOUBLE
    max_error = 0.0_C_DOUBLE
    ok = .true.
    do iter = 0, repeats
      call seed_y_endpoint_system(ny, row_start, active_n, nlines)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call system_clock(t0)
      select case (y_solver)
      case (Y_SOLVER_PIPELINED_LU)
        call ys_solve_pipelined_lu(dst, y_batches)
      case (Y_SOLVER_SCHUR)
        call ys_solve_endpoint_schur(dst, .true.)
      case default
        error stop "unknown autotune y solver"
      end select
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call system_clock(t1)
      !$omp target update from(dst)
      call check_y_endpoint_solution(dst, row_start, row_end, nlines, local_error, bad_local)
      call MPI_Allreduce(local_error, global_error, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(bad_local, bad_global, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      max_error = max(max_error, global_error)
      if (bad_global /= 0 .or. global_error > Y_SOLVE_CHECK_TOL) then
        ok = .false.
        exit
      end if
      if (iter > 0) local_elapsed = local_elapsed + real(t1 - t0, C_DOUBLE)/real(rate, C_DOUBLE)
    end do

    if (ok) then
      local_elapsed = local_elapsed/real(max(1, repeats), C_DOUBLE)
      call MPI_Allreduce(local_elapsed, elapsed, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    else
      elapsed = huge(0.0_C_DOUBLE)
    end if

    !$omp target exit data map(delete: dst)
    deallocate (dst)
    call ys_release_workspace(.true.)
  end subroutine time_y_endpoint_solve

  subroutine seed_y_endpoint_system(ny, row_start, active_n, nlines)
    integer(C_INT), intent(in) :: ny, row_start, active_n, nlines
    complex(C_DOUBLE_COMPLEX) :: coeff(-2:2), rhs
    integer(C_INT) :: irow, iline, p, iy, col

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, ny, row_start, active_n, nlines) &
    !$omp private(irow, iline, p, iy, col, coeff, rhs)
    do irow = 1_C_INT, active_n
      do iline = 1_C_INT, nlines
        iy = row_start + irow - 1_C_INT
        p = (irow - 1_C_INT)*nlines + iline
        coeff(:) = (0.0d0, 0.0d0)
        coeff(0) = (1.25d0, 0.0d0)
        if (iy > 2_C_INT) coeff(-2) = (-0.015d0, 0.0d0)
        if (iy > 1_C_INT) coeff(-1) = (-0.08d0, 0.0d0)
        if (iy < ny - 1_C_INT) coeff(1) = (-0.08d0, 0.0d0)
        if (iy < ny - 2_C_INT) coeff(2) = (-0.015d0, 0.0d0)

        rhs = (0.0d0, 0.0d0)
        do col = -2_C_INT, 2_C_INT
          rhs = rhs + coeff(col)*autotune_exact_value(iy + col, iline)
        end do

        ys_gpsv_ds(p) = coeff(-2)
        ys_gpsv_dl(p) = coeff(-1)
        ys_gpsv_d(p) = coeff(0)
        ys_gpsv_du(p) = coeff(1)
        ys_gpsv_dw(p) = coeff(2)
        ys_gpsv_x(p) = rhs
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine seed_y_endpoint_system

  subroutine check_y_endpoint_solution(dst, row_start, row_end, nlines, local_error, bad)
    complex(C_DOUBLE_COMPLEX), intent(in) :: dst(:, :, :)
    integer(C_INT), intent(in) :: row_start, row_end, nlines
    real(C_DOUBLE), intent(out) :: local_error
    integer, intent(out) :: bad
    complex(C_DOUBLE_COMPLEX) :: expected, got
    integer(C_INT) :: active_n, nlines_z, nx_count, ix, iz_index, irow, iline, iy

    active_n = row_end - row_start + 1_C_INT
    nlines_z = size(dst, 2, kind=C_INT)
    nx_count = nlines/nlines_z
    local_error = 0.0_C_DOUBLE
    bad = 0
    do ix = 1_C_INT, nx_count
      do iz_index = 1_C_INT, nlines_z
        iline = (ix - 1_C_INT)*nlines_z + iz_index
        do irow = 1_C_INT, active_n
          iy = row_start + irow - 1_C_INT
          got = dst(irow + 2_C_INT, iz_index, ix)
          if (.not. finite_complex(got)) then
            bad = 1
            local_error = huge(0.0_C_DOUBLE)
            return
          end if
          expected = autotune_exact_value(iy, iline)
          local_error = max(local_error, abs(got - expected))
        end do
      end do
    end do
  end subroutine check_y_endpoint_solution

  pure logical function finite_complex(value)
    complex(C_DOUBLE_COMPLEX), intent(in) :: value
    finite_complex = ieee_is_finite(real(value, C_DOUBLE)) .and. ieee_is_finite(aimag(value))
  end function finite_complex

  pure complex(C_DOUBLE_COMPLEX) function autotune_exact_value(iy, iline)
    integer(C_INT), intent(in) :: iy, iline
    real(C_DOUBLE) :: yv, lr, li

    yv = real(iy, C_DOUBLE)
    lr = real(mod(iline, 17_C_INT), C_DOUBLE)
    li = real(mod(iline, 13_C_INT), C_DOUBLE)
    autotune_exact_value = cmplx(1.0_C_DOUBLE + 0.031_C_DOUBLE*yv - 0.00021_C_DOUBLE*yv*yv + &
                                 0.000003_C_DOUBLE*yv*yv*yv + 0.0017_C_DOUBLE*lr, &
                                 -0.25_C_DOUBLE + 0.017_C_DOUBLE*yv + 0.00013_C_DOUBLE*yv*yv + &
                                 0.0023_C_DOUBLE*li, kind=C_DOUBLE)
  end function autotune_exact_value

  subroutine node_ids(node)
    integer(C_INT), allocatable, intent(out) :: node(:)
    type(MPI_Comm) :: shared
    integer :: ierr, rank, nranks, local_rank
    integer(C_INT) :: root
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr); call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shared, ierr)
    call MPI_Comm_rank(shared, local_rank, ierr)
    root = merge(int(rank, C_INT), -1_C_INT, local_rank == 0)
    call MPI_Bcast(root, 1, MPI_INTEGER, 0, shared, ierr)
    allocate (node(0:nranks - 1))
    call MPI_Allgather(root, 1, MPI_INTEGER, node, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_Comm_free(shared, ierr)
  end subroutine node_ids

  subroutine abort_msg(rank, message)
    integer, intent(in) :: rank
    character(*), intent(in) :: message
    integer :: ierr
    if (rank == 0) print *, "Error: ", trim(message)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end subroutine abort_msg

  pure logical function mpi_autotune_has_pass_sequence(npy, passes)
    integer(C_INT), intent(in) :: npy, passes(:)
    integer(C_INT) :: i, product
    product = 1_C_INT
    mpi_autotune_has_pass_sequence = (npy == 1_C_INT .and. size(passes) == 0)
    if (size(passes) == 0) return
    do i = 1_C_INT, int(size(passes), C_INT)
      if (.not. any(passes(i) == ARITY)) return
      product = product*passes(i)
    end do
    mpi_autotune_has_pass_sequence = (product == npy)
  end function mpi_autotune_has_pass_sequence

  integer function tune_repeats()
    character(64) :: text
    integer :: status, length, io, value
    call get_environment_variable("CHANNEL_MPI_AUTOTUNE_REPEATS", text, length, status)
    tune_repeats = 2
    if (status == 0) then
      read (text(:length), *, iostat=io) value
      if (io == 0) tune_repeats = max(1, value)
    end if
  end function tune_repeats

  character(128) function pass_string(passes, npass)
    integer(C_INT), intent(in) :: passes(:), npass
    integer(C_INT) :: i
    character(8) :: item
    pass_string = "none"
    if (npass == 0_C_INT) return
    pass_string = ""
    do i = 1_C_INT, npass; write (item, '(I0)') passes(i); if (i > 1_C_INT) pass_string = trim(pass_string)//","; pass_string = trim(pass_string)//trim(item); end do
  end function pass_string

  character(16) function exchange_string(exchange)
    integer(C_INT), intent(in) :: exchange
    select case (exchange)
    case (YS_SCHUR_EXCHANGE_AUTO); exchange_string = "auto"
    case (YS_SCHUR_EXCHANGE_ALLTOALL); exchange_string = "alltoall"
    case (YS_SCHUR_EXCHANGE_ALLGATHER); exchange_string = "allgather"
    case default; exchange_string = "invalid"
    end select
  end function exchange_string

  character(16) function y_solver_string(y_solver)
    integer(C_INT), intent(in) :: y_solver
    select case (y_solver)
    case (Y_SOLVER_SCHUR); y_solver_string = "schur"
    case (Y_SOLVER_PIPELINED_LU); y_solver_string = "pipelined_lu"
    case default; y_solver_string = "invalid"
    end select
  end function y_solver_string

  subroutine print_config(label, npxz, npy, passes, npass, exchange)
    character(*), intent(in) :: label
    integer(C_INT), intent(in) :: npxz, npy, passes(:), npass, exchange
    write (*, '(*(g0,1x))') trim(label)//": npxz=", npxz, "npy=", npy, "passes=", trim(pass_string(passes, npass)), &
      "exchange=", trim(exchange_string(exchange)), "xz_enabled=", npxz > 1_C_INT
  end subroutine print_config
end module mpi_autotune
