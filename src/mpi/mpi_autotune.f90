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
  use y_pipeline_nccl, only: CHANNEL_COMM_BACKEND_AUTO, CHANNEL_COMM_BACKEND_MPI, CHANNEL_COMM_BACKEND_NCCL, &
                             channel_comm_available, channel_comm_backend_from_env, channel_comm_backend_name, &
                             channel_comm_set_backend_override, channel_comm_clear_backend_override
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
  integer(C_INT), save, public :: mpi_autotune_selected_comm_backend = CHANNEL_COMM_BACKEND_AUTO
  integer(C_INT), save :: mpi_autotune_comm_mode = CHANNEL_COMM_BACKEND_AUTO
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
    integer(C_INT) :: best_comm_backend
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
    mpi_autotune_comm_mode = channel_comm_backend_from_env()
    mpi_autotune_selected_comm_backend = CHANNEL_COMM_BACKEND_AUTO
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
      case ("report_only", "REPORT_ONLY", "report-only", "REPORT-ONLY", "scan", "SCAN"); mode = 3
      end select
    end if
    manual = has_npy .or. has_npxz .or. has_passes .or. has_exchange .or. has_y_solver
    found = .false.; applied = .false.; ran_scan = .false.; best_score = huge(0.0_C_DOUBLE)
    if ((mode == 1 .and. .not. manual) .or. mode == 2 .or. mode == 3) then
      ran_scan = .true.
      call node_ids(node)
      call scan(int(nranks, C_INT), nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, best_score, &
                best_npxz, best_npy, best_pass, best_npass, best_exchange, best_comm_backend, &
                best_y_solver, best_y_batches, found)
      if (rank == 0 .and. found) &
        call print_config("MPI autotune recommendation", best_npxz, best_npy, best_pass, best_npass, best_exchange, &
                          best_comm_backend, best_y_solver, best_y_batches)
      if (found .and. (mode == 1 .or. mode == 2) .and. .not. manual) then
        npxz_out = best_npxz; npy_out = best_npy; exchange = best_exchange
        mpi_autotune_selected_y_solver = best_y_solver
        mpi_autotune_selected_y_batches = best_y_batches
        mpi_autotune_selected_comm_backend = best_comm_backend
        call channel_comm_set_backend_override(best_comm_backend)
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
    if (.not. applied) call apply_comm_env_override()
    if (rank == 0 .and. applied) &
      call print_config("MPI autotune selected", npxz_out, npy_out, passes, int(size(passes), C_INT), exchange, &
                        mpi_autotune_selected_comm_backend)
    if (rank == 0 .and. .not. applied .and. .not. (ran_scan .and. mode == 3)) then
      if (ran_scan .and. mode == 2 .and. manual) &
        print *, "MPI autotune report: manual decomposition/backend settings override the recommendation."
      call print_config("MPI decomposition retained", npxz_out, npy_out, passes, int(size(passes), C_INT), exchange)
    end if
  end subroutine configure_mpi_decomposition

  subroutine scan(nranks, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, best_score, &
                  best_npxz, best_npy, best_pass, best_npass, best_exchange, best_comm_backend, &
                  best_y_solver, best_y_batches, found)
    integer(C_INT), intent(in) :: nranks, nxpp, nxd, nzd, nz, ny, nphi, node(0:)
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(out) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange, best_comm_backend
    integer(C_INT), intent(out) :: best_y_solver, best_y_batches
    logical, intent(inout) :: found
    integer(C_INT) :: npxz, path(MAXP), icomm, comm_backend
    path = 1_C_INT
    do npxz = 1_C_INT, nranks
      if (mod(nranks, npxz) == 0_C_INT) then
        if (.not. autotune_npxz_candidate(nranks, npxz, node)) cycle
        do icomm = 1_C_INT, autotune_comm_candidate_count()
          comm_backend = autotune_comm_candidate(icomm)
          call try_pipelined_candidate(nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
                                       best_score, best_npxz, best_npy, best_pass, best_npass, best_exchange, &
                                       best_comm_backend, best_y_solver, best_y_batches, found)
          call gen(nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
                   nranks/npxz, path, 0_C_INT, best_score, best_npxz, best_npy, best_pass, &
                   best_npass, best_exchange, best_comm_backend, best_y_solver, best_y_batches, found)
        end do
      end if
    end do
  end subroutine scan

  subroutine try_pipelined_candidate(nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
                                     best_score, best_npxz, best_npy, best_pass, best_npass, best_exchange, &
                                     best_comm_backend, best_y_solver, best_y_batches, found)
    integer(C_INT), intent(in) :: nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi, node(0:)
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(inout) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange, best_comm_backend
    integer(C_INT), intent(inout) :: best_y_solver, best_y_batches
    logical, intent(inout) :: found
    real(C_DOUBLE) :: xz_forward_cost, xz_back_cost, y_cost, score
    real(C_DOUBLE) :: y_error
    integer(C_INT) :: ib, batches, npy, nlines, placeholder_path(MAXP), placeholder_npass, batch_candidates(3)
    integer(C_INT), allocatable :: placeholder_passes(:)
    logical :: y_ok
    integer :: ierr, rank

    if (.not. valid_decomposition(nranks, npxz, nxpp, nzd, ny, node)) return
    npy = nranks/npxz
    if (npy <= 1_C_INT) return

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call channel_comm_set_backend_override(comm_backend)
    call time_xz_sweep(nxpp, nxd, nzd, nz, ny, nphi, overlapping, npy, xz_forward_cost, xz_back_cost)

    placeholder_path = 1_C_INT
    call ys_schur_default_pass_counts(npy, placeholder_passes)
    placeholder_npass = int(size(placeholder_passes), C_INT)
    if (placeholder_npass > MAXP) error stop "autotune placeholder pass list exceeds MAXP"
    if (placeholder_npass > 0_C_INT) placeholder_path(1:placeholder_npass) = placeholder_passes

    nlines = (nxpp/npxz)*(2_C_INT*nz + 1_C_INT)
    batch_candidates = [max(1_C_INT, npy/2_C_INT), npy, 2_C_INT*npy]
    do ib = 1_C_INT, int(size(batch_candidates), C_INT)
      batches = batch_candidates(ib)
      if (batches > nlines) cycle
      y_cost = time_y(nxpp, nzd, nz, ny, nphi, overlapping, npy, placeholder_path, placeholder_npass, &
                      YS_SCHUR_EXCHANGE_AUTO, Y_SOLVER_PIPELINED_LU, batches, y_ok, y_error)
      if (.not. y_ok) then
        if (rank == 0) write (*, '(*(g0,1x))') "MPI autotune rejected: npxz=", npxz, &
          "npy=", npy, "comm=", trim(channel_comm_backend_name(comm_backend)), &
          "passes=", trim(display_pass_string(placeholder_path, placeholder_npass, Y_SOLVER_PIPELINED_LU)), &
          "exchange= auto y_solver=", trim(y_solver_string(Y_SOLVER_PIPELINED_LU)), &
          "y_batches=", batches, "y_correctness_error=", y_error
        cycle
      end if
      score = RK_SUBSTEPS*(xz_forward_cost*real(3_C_INT + nphi, C_DOUBLE) + &
                           xz_back_cost*real(6_C_INT + 3_C_INT*nphi, C_DOUBLE) + &
                           y_cost*real(3_C_INT + nphi, C_DOUBLE))
      if (rank == 0) write (*, '(*(g0,1x))') "MPI autotune tested: npxz=", npxz, &
        "npy=", npy, "comm=", trim(channel_comm_backend_name(comm_backend)), &
        "passes=", trim(display_pass_string(placeholder_path, placeholder_npass, Y_SOLVER_PIPELINED_LU)), &
        "exchange= auto xz_forward_ms=", 1d3*xz_forward_cost, &
        "xz_back_ms=", 1d3*xz_back_cost, "y_solver=", trim(y_solver_string(Y_SOLVER_PIPELINED_LU)), &
        "y_batches=", batches, "y_ms=", 1d3*y_cost, "y_correctness_error=", y_error, &
        "score_timestep_ms=", 1d3*score
      if (.not. found .or. score < best_score) then
        found = .true.; best_score = score; best_npxz = npxz; best_npy = npy
        best_pass = placeholder_path; best_npass = placeholder_npass; best_exchange = YS_SCHUR_EXCHANGE_AUTO
        best_comm_backend = comm_backend
        best_y_solver = Y_SOLVER_PIPELINED_LU; best_y_batches = batches
      end if
    end do
    if (allocated(placeholder_passes)) deallocate (placeholder_passes)
  end subroutine try_pipelined_candidate

  recursive subroutine gen(nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
                           remaining, path, npass, best_score, best_npxz, best_npy, best_pass, &
                           best_npass, best_exchange, best_comm_backend, best_y_solver, best_y_batches, found)
    integer(C_INT), intent(in) :: nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi
    integer(C_INT), intent(in) :: node(0:), remaining, path(MAXP), npass
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(inout) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange, best_comm_backend
    integer(C_INT), intent(inout) :: best_y_solver, best_y_batches
    logical, intent(inout) :: found
    integer(C_INT) :: i, next_path(MAXP), min_arity
    if (remaining == 1_C_INT) then
      call try_candidate(nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, path, &
                         npass, YS_SCHUR_EXCHANGE_ALLTOALL, best_score, best_npxz, best_npy, &
                         best_pass, best_npass, best_exchange, best_comm_backend, best_y_solver, best_y_batches, found)
      if (npass > 0_C_INT) then
        if (path(npass) < 4_C_INT) &
          call try_candidate(nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, path, &
                             npass, YS_SCHUR_EXCHANGE_ALLGATHER, best_score, best_npxz, best_npy, &
                             best_pass, best_npass, best_exchange, best_comm_backend, best_y_solver, best_y_batches, found)
      end if
      return
    end if
    if (npass >= MAXP) return
    min_arity = min_schur_arity(nranks/npxz)
    do i = 1_C_INT, int(size(ARITY), C_INT)
      if (ARITY(i) < min_arity) cycle
      if (mod(remaining, ARITY(i)) /= 0_C_INT) cycle
      next_path = path; next_path(npass + 1_C_INT) = ARITY(i)
      if (.not. schur_pass_prefix_has_uniform_lines(nxpp, npxz, nz, next_path, npass + 1_C_INT)) cycle
      call gen(nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
               remaining/ARITY(i), next_path, npass + 1_C_INT, best_score, best_npxz, &
               best_npy, best_pass, best_npass, best_exchange, best_comm_backend, best_y_solver, best_y_batches, found)
    end do
  end subroutine gen

  subroutine try_candidate(nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, path, &
                           npass, exchange, best_score, best_npxz, best_npy, best_pass, &
                           best_npass, best_exchange, best_comm_backend, best_y_solver, best_y_batches, found)
    integer(C_INT), intent(in) :: nranks, npxz, comm_backend, nxpp, nxd, nzd, nz, ny, nphi
    integer(C_INT), intent(in) :: node(0:), path(MAXP), npass, exchange
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(inout) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange, best_comm_backend
    integer(C_INT), intent(inout) :: best_y_solver, best_y_batches
    logical, intent(inout) :: found
    integer :: ierr, rank
    real(C_DOUBLE) :: xz_forward_cost, xz_back_cost, y_cost, score
    real(C_DOUBLE) :: y_error
    logical :: y_ok
    if (.not. valid(nranks, npxz, nxpp, nzd, nz, ny, node, path, npass)) return
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call channel_comm_set_backend_override(comm_backend)
    call time_xz_sweep(nxpp, nxd, nzd, nz, ny, nphi, overlapping, nranks/npxz, xz_forward_cost, xz_back_cost)

    call score_y_backend(Y_SOLVER_SCHUR, 0_C_INT)

  contains
    subroutine score_y_backend(y_solver, y_batches)
      integer(C_INT), intent(in) :: y_solver, y_batches

      y_cost = time_y(nxpp, nzd, nz, ny, nphi, overlapping, nranks/npxz, path, npass, exchange, &
                      y_solver, y_batches, y_ok, y_error)
      if (.not. y_ok) then
        if (rank == 0) write (*, '(*(g0,1x))') "MPI autotune rejected: npxz=", npxz, &
          "npy=", nranks/npxz, "comm=", trim(channel_comm_backend_name(comm_backend)), &
          "passes=", trim(display_pass_string(path, npass, y_solver)), &
          "exchange=", trim(exchange_string(exchange)), "y_solver=", trim(y_solver_string(y_solver)), &
          "y_batches=", y_batches, "y_correctness_error=", y_error
        return
      end if
      score = RK_SUBSTEPS*(xz_forward_cost*real(3_C_INT + nphi, C_DOUBLE) + &
                           xz_back_cost*real(6_C_INT + 3_C_INT*nphi, C_DOUBLE) + &
                           y_cost*real(3_C_INT + nphi, C_DOUBLE))
      if (rank == 0) write (*, '(*(g0,1x))') "MPI autotune tested: npxz=", npxz, &
        "npy=", nranks/npxz, "comm=", trim(channel_comm_backend_name(comm_backend)), &
        "passes=", trim(display_pass_string(path, npass, y_solver)), &
        "exchange=", trim(exchange_string(exchange)), "xz_forward_ms=", 1d3*xz_forward_cost, &
        "xz_back_ms=", 1d3*xz_back_cost, "y_solver=", trim(y_solver_string(y_solver)), &
        "y_batches=", y_batches, "y_ms=", 1d3*y_cost, "y_correctness_error=", y_error, &
        "score_timestep_ms=", 1d3*score
      if (.not. found .or. score < best_score) then
        found = .true.; best_score = score; best_npxz = npxz; best_npy = nranks/npxz
        best_pass = path; best_npass = npass; best_exchange = exchange
        best_comm_backend = comm_backend
        best_y_solver = y_solver; best_y_batches = y_batches
      end if
    end subroutine score_y_backend
  end subroutine try_candidate

  logical function valid(nranks, npxz, nxpp, nzd, nz, ny, node, path, npass)
    integer(C_INT), intent(in) :: nranks, npxz, nxpp, nzd, nz, ny, node(0:), path(MAXP), npass
    integer(C_INT) :: npy
    valid = .false.; npy = nranks/npxz
    if (.not. valid_decomposition(nranks, npxz, nxpp, nzd, ny, node)) return
    if (.not. schur_pass_prefix_has_uniform_lines(nxpp, npxz, nz, path, npass)) return
    valid = .true.
  end function valid

  logical function valid_decomposition(nranks, npxz, nxpp, nzd, ny, node)
    integer(C_INT), intent(in) :: nranks, npxz, nxpp, nzd, ny, node(0:)
    integer(C_INT) :: npy, ipy, y0, yN
    valid_decomposition = .false.
    npy = nranks/npxz
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
    valid_decomposition = .true.
  end function valid_decomposition

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

  logical function autotune_npxz_candidate(nranks, npxz, node)
    integer(C_INT), intent(in) :: nranks, npxz, node(0:)
    integer(C_INT) :: node_npxz

    node_npxz = autotune_node_size(nranks, node)
    autotune_npxz_candidate = (npxz == 1_C_INT .or. npxz == node_npxz)
  end function autotune_npxz_candidate

  integer(C_INT) function autotune_node_size(nranks, node)
    integer(C_INT), intent(in) :: nranks, node(0:)
    integer(C_INT) :: i, j, count

    autotune_node_size = 1_C_INT
    do i = 0_C_INT, nranks - 1_C_INT
      count = 0_C_INT
      do j = 0_C_INT, nranks - 1_C_INT
        if (node(j) == node(i)) count = count + 1_C_INT
      end do
      autotune_node_size = max(autotune_node_size, count)
    end do
  end function autotune_node_size

  integer(C_INT) function min_schur_arity(npy)
    integer(C_INT), intent(in) :: npy
    min_schur_arity = 1_C_INT
    do while (min_schur_arity*min_schur_arity*min_schur_arity < npy)
      min_schur_arity = min_schur_arity + 1_C_INT
    end do
  end function min_schur_arity

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
        call repack_zTOx_local(VVdz(:, :, :, 1), VVdx(:, :, :, 1))
      else
        call pack_zTOx(VVdz(:, :, :, 1), sendbuf(:, 1))
        call alltoall(sendbuf(:, 1), recvbuf(:, 1), request, "zTOx autotune_xz_sweep")
        call MPI_Wait(request, status, ierr)
        call unpack_zTOx(recvbuf(:, 1), VVdx(:, :, :, 1))
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (iter > 0) forward_elapsed = forward_elapsed + MPI_Wtime() - t0

      t0 = MPI_Wtime()
      if (fft_transpose_is_local) then
        call repack_xTOz_local(VVdx(:, :, :, 1), VVdz(:, :, :, 1))
      else
        call pack_xTOz(VVdx(:, :, :, 1), sendbuf(:, 1))
        call alltoall(sendbuf(:, 1), recvbuf(:, 1), request, "xTOz autotune_xz_sweep")
        call MPI_Wait(request, status, ierr)
        call unpack_xTOz(recvbuf(:, 1), VVdz(:, :, :, 1))
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
    complex(C_DOUBLE_COMPLEX), allocatable :: dst(:, :, :), expected(:, :, :)
    integer(C_INT) :: active_n, nlines_z, nx_count
    integer(C_INT64_T) :: t0, t1, rate
    integer :: ierr, iter, bad_local, bad_global
    real(C_DOUBLE) :: local_elapsed, local_error, global_error

    call ys_prepare_assembled_workspace(ny, nz, row_start, row_end, line_start, nlines, passes, exchange, &
                                        prepare_schur=(y_solver == Y_SOLVER_SCHUR))
    active_n = row_end - row_start + 1_C_INT
    nlines_z = 2_C_INT*nz + 1_C_INT
    nx_count = nlines/nlines_z
    allocate (dst(active_n + 4_C_INT, 2_C_INT*nz + 1_C_INT, nx_count))
    allocate (expected(active_n, nlines_z, nx_count))
    call fill_y_endpoint_expected(expected, row_start, nlines)
    !$omp target enter data map(alloc: dst)
    !$omp target enter data map(to: expected)

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
      call check_y_endpoint_solution_device(dst, expected, local_error, bad_local)
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
    !$omp target exit data map(delete: expected)
    deallocate (dst)
    deallocate (expected)
    call ys_release_workspace(.true.)
  end subroutine time_y_endpoint_solve

  subroutine fill_y_endpoint_expected(expected, row_start, nlines)
    complex(C_DOUBLE_COMPLEX), intent(out) :: expected(:, :, :)
    integer(C_INT), intent(in) :: row_start, nlines
    integer(C_INT) :: active_n, nlines_z, nx_count, ix, iz_index, irow, iline, iy

    active_n = size(expected, 1, kind=C_INT)
    nlines_z = size(expected, 2, kind=C_INT)
    nx_count = nlines/nlines_z
    do ix = 1_C_INT, nx_count
      do iz_index = 1_C_INT, nlines_z
        iline = (ix - 1_C_INT)*nlines_z + iz_index
        do irow = 1_C_INT, active_n
          iy = row_start + irow - 1_C_INT
          expected(irow, iz_index, ix) = autotune_exact_value(iy, iline)
        end do
      end do
    end do
  end subroutine fill_y_endpoint_expected

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

  subroutine check_y_endpoint_solution_device(dst, expected, local_error, bad)
    complex(C_DOUBLE_COMPLEX), intent(in) :: dst(:, :, :), expected(:, :, :)
    real(C_DOUBLE), intent(out) :: local_error
    integer, intent(out) :: bad
    complex(C_DOUBLE_COMPLEX) :: got
    real(C_DOUBLE) :: got_real, got_imag, err
    integer(C_INT) :: active_n, nlines_z, nx_count, ix, iz_index, irow

    active_n = size(expected, 1, kind=C_INT)
    nlines_z = size(dst, 2, kind=C_INT)
    nx_count = size(expected, 3, kind=C_INT)
    local_error = 0.0_C_DOUBLE
    bad = 0
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(dst, expected, active_n, nlines_z, nx_count) &
    !$omp private(ix, iz_index, irow, got, got_real, got_imag, err) reduction(max:local_error, bad)
    do ix = 1_C_INT, nx_count
      do iz_index = 1_C_INT, nlines_z
        do irow = 1_C_INT, active_n
          got = dst(irow + 2_C_INT, iz_index, ix)
          got_real = real(got, C_DOUBLE)
          got_imag = aimag(got)
          if (got_real /= got_real .or. got_imag /= got_imag .or. &
              abs(got_real) > huge(0.0_C_DOUBLE) .or. abs(got_imag) > huge(0.0_C_DOUBLE)) then
            bad = 1
            local_error = huge(0.0_C_DOUBLE)
          else
            err = abs(got - expected(irow, iz_index, ix))
            local_error = max(local_error, err)
          end if
        end do
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine check_y_endpoint_solution_device

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

  pure logical function schur_pass_prefix_has_uniform_lines(nxpp, npxz, nz, passes, npass)
    integer(C_INT), intent(in) :: nxpp, npxz, nz, passes(MAXP), npass
    integer(C_INT) :: level, nlines

    schur_pass_prefix_has_uniform_lines = .false.
    if (npxz < 1_C_INT .or. mod(nxpp, npxz) /= 0_C_INT) return
    nlines = (nxpp/npxz)*(2_C_INT*nz + 1_C_INT)
    do level = 1_C_INT, npass
      if (mod(nlines, passes(level)) /= 0_C_INT) return
      nlines = nlines/passes(level)
    end do
    schur_pass_prefix_has_uniform_lines = .true.
  end function schur_pass_prefix_has_uniform_lines

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

  character(128) function display_pass_string(passes, npass, y_solver)
    integer(C_INT), intent(in) :: passes(:), npass, y_solver

    if (y_solver == Y_SOLVER_PIPELINED_LU) then
      display_pass_string = "none"
    else
      display_pass_string = pass_string(passes, npass)
    end if
  end function display_pass_string

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

  integer(C_INT) function autotune_comm_candidate_count()
    if (mpi_autotune_comm_mode == CHANNEL_COMM_BACKEND_AUTO .and. channel_comm_available()) then
      autotune_comm_candidate_count = 2_C_INT
    else
      autotune_comm_candidate_count = 1_C_INT
    end if
  end function autotune_comm_candidate_count

  integer(C_INT) function autotune_comm_candidate(index)
    integer(C_INT), intent(in) :: index

    select case (mpi_autotune_comm_mode)
    case (CHANNEL_COMM_BACKEND_MPI, CHANNEL_COMM_BACKEND_NCCL)
      autotune_comm_candidate = mpi_autotune_comm_mode
    case default
      if (index == 2_C_INT .and. channel_comm_available()) then
        autotune_comm_candidate = CHANNEL_COMM_BACKEND_NCCL
      else
        autotune_comm_candidate = CHANNEL_COMM_BACKEND_MPI
      end if
    end select
  end function autotune_comm_candidate

  subroutine apply_comm_env_override()
    select case (mpi_autotune_comm_mode)
    case (CHANNEL_COMM_BACKEND_MPI, CHANNEL_COMM_BACKEND_NCCL)
      call channel_comm_set_backend_override(mpi_autotune_comm_mode)
      mpi_autotune_selected_comm_backend = mpi_autotune_comm_mode
    case default
      call channel_comm_clear_backend_override()
      mpi_autotune_selected_comm_backend = CHANNEL_COMM_BACKEND_AUTO
    end select
  end subroutine apply_comm_env_override

  subroutine print_config(label, npxz, npy, passes, npass, exchange, comm_backend, y_solver, y_batches)
    character(*), intent(in) :: label
    integer(C_INT), intent(in) :: npxz, npy, passes(:), npass, exchange
    integer(C_INT), intent(in), optional :: comm_backend
    integer(C_INT), intent(in), optional :: y_solver, y_batches
    integer(C_INT) :: solver, batches, backend
    solver = mpi_autotune_selected_y_solver
    batches = mpi_autotune_selected_y_batches
    backend = mpi_autotune_selected_comm_backend
    if (present(comm_backend)) backend = comm_backend
    if (present(y_solver)) solver = y_solver
    if (present(y_batches)) batches = y_batches
    write (*, '(*(g0,1x))') trim(label)//": npxz=", npxz, "npy=", npy, &
      "passes=", trim(display_pass_string(passes, npass, solver)), &
      "exchange=", trim(exchange_string(exchange)), "comm=", trim(channel_comm_backend_name(backend)), &
      "xz_enabled=", npxz > 1_C_INT, &
      "y_solver=", trim(y_solver_string(solver)), "y_batches=", batches
  end subroutine print_config
end module mpi_autotune
