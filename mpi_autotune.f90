#include "header.h"

module mpi_autotune
  use, intrinsic :: iso_c_binding
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
                            ys_solve_endpoint_schur, ys_gpsv_ds, ys_gpsv_dl, &
                            ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x
  use y_schur_solver, only: ys_schur_default_pass_counts, YS_SCHUR_EXCHANGE_AUTO, &
                            YS_SCHUR_EXCHANGE_ALLTOALL, YS_SCHUR_EXCHANGE_ALLGATHER
  use mpi_f08
  implicit none
  private
  integer(C_INT), parameter :: MAXP = 16_C_INT
  integer(C_INT), parameter :: ARITY(5) = [2_C_INT, 3_C_INT, 4_C_INT, 6_C_INT, 8_C_INT]
  real(C_DOUBLE), parameter :: RK_SUBSTEPS = 3.0_C_DOUBLE
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
    integer(C_INT), allocatable :: node(:)
    logical :: has_npy, has_npxz, has_passes, has_exchange, manual, applied, found
    real(C_DOUBLE) :: best_score

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

    npy_out = max(1_C_INT, requested_npy); npxz_out = 1_C_INT
    call get_environment_variable("CHANNEL_NPY", text, length, status); has_npy = (status == 0); env_npy = 0_C_INT
    if (has_npy) read (text(:length), *, iostat=io) env_npy
    call get_environment_variable("CHANNEL_NPXZ", text, length, status); has_npxz = (status == 0); env_npxz = 0_C_INT
    if (has_npxz) read (text(:length), *, iostat=io) env_npxz
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
      if (mod(int(nranks, C_INT), npy_out) /= 0_C_INT) call abort_msg(rank, "input npy must divide rank count")
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
    manual = has_npy .or. has_npxz .or. has_passes .or. has_exchange
    found = .false.; applied = .false.; best_score = huge(0.0_C_DOUBLE)
    if ((mode == 1 .and. .not. manual) .or. mode == 2) then
      call node_ids(node)
      call scan(int(nranks, C_INT), nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, best_score, &
                best_npxz, best_npy, best_pass, best_npass, best_exchange, found)
      if (rank == 0 .and. found) &
        call print_config("MPI autotune recommendation", best_npxz, best_npy, best_pass, best_npass, best_exchange)
      if (found .and. mode == 1 .and. .not. manual) then
        npxz_out = best_npxz; npy_out = best_npy; exchange = best_exchange; applied = .true.
      else if (.not. found .and. rank == 0) then
        print *, "Warning: MPI autotune found no valid candidates; keeping configured decomposition."
      end if
      deallocate (node)
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
                  best_npxz, best_npy, best_pass, best_npass, best_exchange, found)
    integer(C_INT), intent(in) :: nranks, nxpp, nxd, nzd, nz, ny, nphi, node(0:)
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(out) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange
    logical, intent(inout) :: found
    integer(C_INT) :: npxz, path(MAXP)
    path = 1_C_INT
    do npxz = 1_C_INT, nranks
      if (mod(nranks, npxz) == 0_C_INT) &
        call gen(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
                 nranks/npxz, path, 0_C_INT, best_score, best_npxz, best_npy, best_pass, &
                 best_npass, best_exchange, found)
    end do
  end subroutine scan

  recursive subroutine gen(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
                           remaining, path, npass, best_score, best_npxz, best_npy, best_pass, &
                           best_npass, best_exchange, found)
    integer(C_INT), intent(in) :: nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, node(0:), remaining, path(MAXP), npass
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(inout) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange
    logical, intent(inout) :: found
    integer(C_INT) :: i, next_path(MAXP)
    if (remaining == 1_C_INT) then
      call try_candidate(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, path, &
                         npass, YS_SCHUR_EXCHANGE_ALLTOALL, best_score, best_npxz, best_npy, &
                         best_pass, best_npass, best_exchange, found)
      if (npass > 0_C_INT) then
        if (path(npass) < 4_C_INT) &
          call try_candidate(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, path, &
                             npass, YS_SCHUR_EXCHANGE_ALLGATHER, best_score, best_npxz, best_npy, &
                             best_pass, best_npass, best_exchange, found)
      end if
      return
    end if
    if (npass >= MAXP) return
    do i = 1_C_INT, int(size(ARITY), C_INT)
      if (mod(remaining, ARITY(i)) /= 0_C_INT) cycle
      next_path = path; next_path(npass + 1_C_INT) = ARITY(i)
      call gen(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, &
               remaining/ARITY(i), next_path, npass + 1_C_INT, best_score, best_npxz, &
               best_npy, best_pass, best_npass, best_exchange, found)
    end do
  end subroutine gen

  subroutine try_candidate(nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, overlapping, node, path, &
                           npass, exchange, best_score, best_npxz, best_npy, best_pass, &
                           best_npass, best_exchange, found)
    integer(C_INT), intent(in) :: nranks, npxz, nxpp, nxd, nzd, nz, ny, nphi, node(0:), path(MAXP), npass, exchange
    logical, intent(in) :: overlapping
    real(C_DOUBLE), intent(inout) :: best_score
    integer(C_INT), intent(inout) :: best_npxz, best_npy, best_pass(MAXP), best_npass, best_exchange
    logical, intent(inout) :: found
    integer :: ierr, rank
    real(C_DOUBLE) :: xz_forward_cost, xz_back_cost, y_cost, score
    if (.not. valid(nranks, npxz, nxpp, nzd, nz, node, path, npass)) return
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call time_xz_sweep(nxpp, nxd, nzd, nz, ny, nphi, overlapping, nranks/npxz, xz_forward_cost, xz_back_cost)
    y_cost = time_y(nxpp, nzd, nz, ny, nphi, overlapping, nranks/npxz, path, npass, exchange)
    score = RK_SUBSTEPS*(xz_forward_cost*real(3_C_INT + nphi, C_DOUBLE) + &
                         xz_back_cost*real(6_C_INT + 3_C_INT*nphi, C_DOUBLE) + &
                         y_cost*real(3_C_INT + nphi, C_DOUBLE))
    if (rank == 0) write (*, '(*(g0,1x))') "MPI autotune tested: npxz=", npxz, &
      "npy=", nranks/npxz, "passes=", trim(pass_string(path, npass)), &
      "exchange=", trim(exchange_string(exchange)), "xz_forward_ms=", 1d3*xz_forward_cost, &
      "xz_back_ms=", 1d3*xz_back_cost, &
      "y_schur_ms=", 1d3*y_cost, "score_timestep_ms=", 1d3*score
    if (.not. found .or. score < best_score) then
      found = .true.; best_score = score; best_npxz = npxz; best_npy = nranks/npxz
      best_pass = path; best_npass = npass; best_exchange = exchange
    end if
  end subroutine try_candidate

  logical function valid(nranks, npxz, nxpp, nzd, nz, node, path, npass)
    integer(C_INT), intent(in) :: nranks, npxz, nxpp, nzd, nz, node(0:), path(MAXP), npass
    integer(C_INT) :: npy, ipy, ipxz, prev, level
    valid = .false.; npy = nranks/npxz
    if (mod(nxpp, npxz) /= 0_C_INT .or. mod(nzd, npxz) /= 0_C_INT) return
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

  real(C_DOUBLE) function time_y(nxpp, nzd, nz, ny, nphi, overlapping, npy, path, npass, exchange)
    integer(C_INT), intent(in) :: nxpp, nzd, nz, ny, nphi, npy, path(MAXP), npass, exchange
    logical, intent(in) :: overlapping
    integer(C_INT), allocatable :: passes(:)
    integer(C_INT) :: nlines

    call init_MPI(nxpp, nz, ny, nzd, nphi, overlapping, npy, .true.)
    allocate (passes(npass))
    if (npass > 0_C_INT) passes = path(1:npass)
    nlines = nxB*(2_C_INT*nz + 1_C_INT)
    call time_y_endpoint_solve(ny, nz, ny0, nyN, 1_C_INT, nlines, passes, exchange, &
                               tune_repeats(), time_y)
    deallocate (passes)
    call free_MPI()
  end function time_y

  subroutine time_y_endpoint_solve(ny, nz, row_start, row_end, line_start, nlines, &
                                   passes, exchange, repeats, elapsed)
    integer(C_INT), intent(in) :: ny, nz, row_start, row_end, line_start, nlines
    integer(C_INT), intent(in) :: passes(:), exchange
    integer, intent(in) :: repeats
    real(C_DOUBLE), intent(out) :: elapsed
    complex(C_DOUBLE_COMPLEX), allocatable :: dst(:, :, :)
    integer(C_INT) :: active_n, nx_count
    integer(C_INT64_T) :: t0, t1, rate
    integer :: ierr, iter
    real(C_DOUBLE) :: local_elapsed

    call ys_prepare_assembled_workspace(ny, nz, row_start, row_end, line_start, nlines, passes, exchange)
    active_n = row_end - row_start + 1_C_INT
    nx_count = nlines/(2_C_INT*nz + 1_C_INT)
    allocate (dst(active_n + 4_C_INT, 2_C_INT*nz + 1_C_INT, nx_count))
    !$omp target enter data map(alloc: dst)

    call system_clock(count_rate=rate)
    local_elapsed = 0.0_C_DOUBLE
    do iter = 0, repeats
      call seed_y_endpoint_system(active_n, nlines)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call system_clock(t0)
      call ys_solve_endpoint_schur(dst, .true.)
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call system_clock(t1)
      if (iter > 0) local_elapsed = local_elapsed + real(t1 - t0, C_DOUBLE)/real(rate, C_DOUBLE)
    end do

    local_elapsed = local_elapsed/real(max(1, repeats), C_DOUBLE)
    call MPI_Allreduce(local_elapsed, elapsed, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    !$omp target exit data map(delete: dst)
    deallocate (dst)
    call ys_release_workspace()
  end subroutine time_y_endpoint_solve

  subroutine seed_y_endpoint_system(active_n, nlines)
    integer(C_INT), intent(in) :: active_n, nlines
    integer(C_INT) :: irow, iline, p

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, active_n, nlines) &
    !$omp private(irow, iline, p)
    do irow = 1_C_INT, active_n
      do iline = 1_C_INT, nlines
        p = (irow - 1_C_INT)*nlines + iline
        ys_gpsv_ds(p) = (0.0d0, 0.0d0)
        ys_gpsv_dl(p) = (0.0d0, 0.0d0)
        ys_gpsv_d(p) = (1.25d0, 0.0d0)
        ys_gpsv_du(p) = (0.0d0, 0.0d0)
        ys_gpsv_dw(p) = (0.0d0, 0.0d0)
        if (irow > 2_C_INT) ys_gpsv_ds(p) = (-0.015d0, 0.0d0)
        if (irow > 1_C_INT) ys_gpsv_dl(p) = (-0.08d0, 0.0d0)
        if (irow < active_n) ys_gpsv_du(p) = (-0.08d0, 0.0d0)
        if (irow < active_n - 1_C_INT) ys_gpsv_dw(p) = (-0.015d0, 0.0d0)
        ys_gpsv_x(p) = cmplx(1.0_C_DOUBLE + 0.001_C_DOUBLE*real(mod(iline, 17_C_INT), C_DOUBLE), &
                             0.002_C_DOUBLE*real(mod(irow, 13_C_INT), C_DOUBLE), kind=C_DOUBLE)
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine seed_y_endpoint_system

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

  subroutine print_config(label, npxz, npy, passes, npass, exchange)
    character(*), intent(in) :: label
    integer(C_INT), intent(in) :: npxz, npy, passes(:), npass, exchange
    write (*, '(*(g0,1x))') trim(label)//": npxz=", npxz, "npy=", npy, "passes=", trim(pass_string(passes, npass)), &
      "exchange=", trim(exchange_string(exchange)), "xz_enabled=", npxz > 1_C_INT
  end subroutine print_config
end module mpi_autotune
