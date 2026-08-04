!============================================!
!                                            !
!     Direct Numerical Simulation (DNS)      !
!       of a turbulent channel flow          !
!                                            !
!============================================!
!
! This program has been written following the
! KISS (Keep it Simple and Stupid) philosophy
!
! Author: Dr.-Ing. Davide Gatti
!

#include "build_options.h"

MODULE driver
  USE roctx, ONLY: roctxPush, roctxPop
  USE env_options, ONLY: env_flag, env_int

CONTAINS
  !==========================================================
  SUBROUTINE initialize(config_file, restart_file, solveNS)
    use config, only: ini_config, read_ini_file
    USE case_setup
    USE mpi_transpose, only: init_MPI
    USE channel_transforms, only: transform_to_physical, compute_cfl
    USE mpi_autotune, only: configure_mpi_decomposition
    USE convvelo, only: init_convvelo_runtime, get_convvelo_memory_estimate, get_convvelo_workspace_estimate, configure_convvelo
    USE ffts, only: get_fft_memory_estimate
#ifdef HAVE_CUDA
    USE ffts, only: init_cufft, free_fft, acquire_fft_workspace, release_fft_workspace, get_fft_workspace_bytes_for_dims
#elif defined(HAVE_HIP)
    USE ffts, only: init_hipfft, free_fft, acquire_fft_workspace, release_fft_workspace, get_fft_workspace_bytes_for_dims
#else
    USE ffts, only: init_fft, free_fft
#endif
USE pressure_output, only: init_pressure_output, free_pressure_output, get_pressure_memory_estimate, get_pressure_workspace_estimate
    USE y_line_solvers, only: ys_get_gpusparse_buffer_bytes
    use byte_workspace, only: workspace_reserve
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    use omp_lib
#endif
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: config_file, restart_file
    LOGICAL, OPTIONAL, INTENT(IN) :: solveNS
    type(ini_config) :: cfg
    REAL(C_DOUBLE) :: deltat_from_dnsin
    real(C_DOUBLE) :: total_mib
    integer(C_INT64_T) :: solver_floats, fft_floats, pressure_floats, convvelo_floats, external_floats, persistent_floats
    integer(C_SIZE_T) :: solver_workspace_bytes, fft_workspace_bytes, pressure_workspace_bytes, convvelo_workspace_bytes
    integer(C_SIZE_T) :: workspace_peak_bytes, sparse_external_bytes
    integer :: iPhi
    integer(C_INT) :: npxz_tuned
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    integer :: num_dev, dev, local_rank, num_dev_min, num_dev_max
#endif
    logical :: run_solver, exit_after_autotune
    complex(C_DOUBLE_COMPLEX), allocatable :: zero_mode(:)

    run_solver = .true.
    if (present(solveNS)) run_solver = solveNS

    ! Init MPI
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    num_dev = omp_get_num_devices()
    call MPI_Allreduce(num_dev, num_dev_min, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(num_dev, num_dev_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (num_dev_min < 1 .or. num_dev_min /= num_dev_max) then
      if (num_dev < 1 .or. num_dev == num_dev_min) then
        print *, 'ERROR: rank', iproc, 'sees', num_dev, 'OpenMP target devices'
      end if
      if (iproc == 0) then
        print *, 'ERROR: inconsistent OpenMP target device count across ranks.'
        print *, '       Minimum devices =', num_dev_min, ' maximum devices =', num_dev_max
      end if
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    local_rank = mpi_local_rank_from_env()
    if (local_rank >= 0) then
      dev = mod(local_rank, num_dev)
    else
      dev = mod(iproc, num_dev)
    end if

    call omp_set_default_device(dev)

    print *, 'Rank', iproc, 'of', nproc, 'local rank', local_rank, 'using device', dev, 'out of', num_dev

    !$omp target
    print *, 'Hello from GPU on rank', iproc, 'device', dev
    !$omp end target
#endif

    call read_ini_file(config_file, cfg)
    CALL read_dnsin(cfg)

    ! Choose the rank decomposition and the y-Schur pass hierarchy.  This may
    ! run a timing scan, so it is driven from here rather than from read_dnsin.
    if (allocated(schur_pass_counts)) deallocate (schur_pass_counts)
    call configure_mpi_decomposition(nx + 1, nxd, nzd, nz, ny, nPhi, overlapping, 1_C_INT, &
                                     npy, npxz_tuned, schur_pass_counts, schur_exchange_mode)

    exit_after_autotune = .false.
    call env_flag("CHANNEL_EXIT_AFTER_MPI_AUTOTUNE", exit_after_autotune)
    if (exit_after_autotune) then
      if (iproc == 0) print *, "CHANNEL_EXIT_AFTER_MPI_AUTOTUNE set; exiting after MPI autotune/configuration."
      CALL MPI_FINALIZE(ierr)
      stop
    end if
    call configure_convvelo(cfg)
    deltat_from_dnsin = deltat
    CALL init_MPI(nx + 1, nz, ny, nzd, nPhi, overlapping, npy)
    call get_solver_memory_estimate(run_solver, solver_floats)
    call get_fft_memory_estimate(nxd, nxB, nzd, nzB, nPhi, overlapping, fft_floats)
    call get_pressure_memory_estimate(pressure_floats)
    call get_convvelo_memory_estimate(convvelo_floats)
    call get_mpi_buffer_memory_estimate(external_floats)
    call get_solver_workspace_estimate(run_solver, solver_workspace_bytes)
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    call get_fft_workspace_bytes_for_dims(nxd, nxB, nzd, nzB, nPhi, overlapping, fft_workspace_bytes)
    fft_floats = 0_C_INT64_T
#else
    fft_workspace_bytes = 0_C_SIZE_T
#endif
    call get_pressure_workspace_estimate(pressure_workspace_bytes)
    call get_convvelo_workspace_estimate(convvelo_workspace_bytes)
    call ys_get_gpusparse_buffer_bytes(sparse_external_bytes)
    workspace_peak_bytes = max(solver_workspace_bytes, fft_workspace_bytes)
    workspace_peak_bytes = max(workspace_peak_bytes, pressure_workspace_bytes)
    workspace_peak_bytes = max(workspace_peak_bytes, convvelo_workspace_bytes)
    persistent_floats = solver_floats + fft_floats + pressure_floats + convvelo_floats
    if (has_terminal) then
      write (*, *) "Estimated memory per rank before allocation:"
      call print_memory_line("Persistent solver", solver_floats)
      call print_memory_line("Persistent FFT", fft_floats)
      call print_memory_line("Persistent pressure", pressure_floats)
      call print_memory_line("Persistent convvelo", convvelo_floats)
      call print_memory_line("Known external MPI", external_floats)
      write (*, '(A,": device=",F12.3,A)') "  Known external sparse", &
        real(sparse_external_bytes, C_DOUBLE)/(1024.0d0*1024.0d0), " MiB"
      write (*, '(A,": device=",F12.3,A)') "  Workspace peak", &
        real(workspace_peak_bytes, C_DOUBLE)/(1024.0d0*1024.0d0), " MiB"
      total_mib = floats_to_mib(persistent_floats + external_floats) + &
                  real(sparse_external_bytes + workspace_peak_bytes, C_DOUBLE)/(1024.0d0*1024.0d0)
      write (*, '(A,F12.3,A)') "  Total peak : device=", total_mib, " MiB"
      write (*, '(A)') "  Note       : sparse-library gpsv buffers are updated after their first size query."
      write (*, *) " "
    end if
    CALL init_memory(run_solver)

    ! Init various subroutines
#ifdef HAVE_CUDA
    CALL init_cufft(nxd, nxB, nzd, nzB, nPhi, overlapping)
#elif defined(HAVE_HIP)
    CALL init_hipfft(nxd, nxB, nzd, nzB, nPhi, overlapping)
#elif defined(HAVE_FFTW)
    CALL init_fft(nxd, nxB, nzd, nzB, nPhi, overlapping)
#endif
    CALL setup_derivatives()
    CALL setup_boundary_conditions()
    CALL read_restart_file(restart_file, V)
    if (.not. run_solver) then
      !$omp target update to(V)
    end if
    CALL init_pressure_output()
    call init_convvelo_runtime()
    call workspace_reserve(workspace_peak_bytes)

    ! Field number (for output)
    ifield = FLOOR((time + 0.5*deltat)/dt_field)
    time0 = time

    ! Reset desired dt
    IF (cflmax == 0) THEN
      deltat = deltat_from_dnsin
    END IF

    IF (has_terminal) THEN
      ! Output DNS.in
      WRITE (*, *) " "
      WRITE (*, *) "!====================================================!"
      WRITE (*, *) "!                     D   N   S                      !"
      WRITE (*, *) "!====================================================!"
      WRITE (*, *) " "
      WRITE (*, "(A,I5,A,I5,A,I5)") "   nx =", nx, "   ny =", ny, "   nz =", nz
      WRITE (*, "(A,I5,A,I5)") "   nxd =", nxd, "  nzd =", nzd
      WRITE (*, "(A,F11.6,A,F11.6,A,F8.6)") "   alfa0 =", alfa0, "       beta0 =", beta0, "   ni =", ni
      WRITE (*, "(A,F11.6,A,F11.6)") "   meanpx =", meanpx, "      meanpz =", meanpz
      WRITE (*, "(A,F11.6,A,F11.6)") "   meanflowx =", meanflowx, "   meanflowz =", meanflowz
      WRITE (*, "(A,I6,A,L1)") "   nsteps =", nstep, "   time_from_restart =", time_from_restart
      WRITE (*, *) "NUM_SCALARS", nPhi, "PRANDTL NUMBERS:"
      do iPhi = 1, nPhi
        write (*, '(F10.4)') 1/pra(iPhi)
      end do
      WRITE (*, *) " "

      print *, "Overlapping communication and computation:", overlapping
    END IF

    allocate (zero_mode(-1:ny + 1))

    if (run_solver) then
      ! Compute CFL
      if (deltat == 0.0) deltat = 1.0
      !$omp target update to(V)
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
      call acquire_fft_workspace("init_fft")
#endif
      CALL transform_to_physical()
      call compute_cfl()
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
      call release_fft_workspace("init_fft")
#endif
      ! cfl is still the rank-local maximum here (outstats reduces it), so
      ! print it from one rank only rather than once per rank.
      IF (has_terminal) print *, "CFL", deltat, cfl
      ! Compute flow rate
      IF (has_average) THEN
        call gather_full_y_line(ny, V(:, 0, 0, 1), zero_mode); fr(1) = yintegr(zero_mode, y); 
        call gather_full_y_line(ny, V(:, 0, 0, 3), zero_mode); fr(2) = yintegr(zero_mode, y); 
        DO iPhi = 1, nPhi
          call gather_full_y_line(ny, V(:, 0, 0, 3 + iPhi), zero_mode); fr(3 + iPhi) = yintegr(zero_mode, y); 
        END DO
      END IF
      CALL outstats()
    end if
  END SUBROUTINE initialize

  ! Node-local rank as published by whichever launcher started us; -1 when none
  ! of them did, in which case the caller falls back to the global rank.
  integer function mpi_local_rank_from_env()
    use, intrinsic :: iso_c_binding, only: C_INT
    implicit none
    character(len=*), parameter :: LAUNCHER_VARS(5) = [character(len=26) :: &
                                                       "OMPI_COMM_WORLD_LOCAL_RANK", &
                                                       "MPI_LOCALRANKID", &
                                                       "MV2_COMM_WORLD_LOCAL_RANK", &
                                                       "SLURM_LOCALID", &
                                                       "PMI_LOCAL_RANK"]
    integer(C_INT) :: parsed
    integer :: i

    do i = 1, size(LAUNCHER_VARS)
      if (env_int(trim(LAUNCHER_VARS(i)), parsed)) then
        mpi_local_rank_from_env = int(parsed)
        return
      end if
    end do
    mpi_local_rank_from_env = -1
  end function mpi_local_rank_from_env

  !==========================================================
  SUBROUTINE timeloop()
    USE case_setup
    USE channel_transforms, only: transform_to_physical, transform_back_and_build_rhs, compute_cfl
    USE convvelo, only: advance_convvelo_runtime
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    USE ffts, only: acquire_fft_workspace, release_fft_workspace
#endif
    IMPLICIT NONE
    integer:: iPhi, ix, iz, i
#ifdef chron
    REAL timei, timee, elapsed_run_time

    elapsed_run_time = 0.0
#endif

    DO WHILE ((time < t_max - deltat/2.0) .AND. (istep < nstep))
#ifdef chron
      CALL CPU_TIME(timei)
#endif
      call roctxPush("timestep")
      ! apply boundary conditions from input file (Couette-like)
      call roctxPush("boundary_conditions")
      IF (has_average) THEN
        !$omp target
        bc0(0, 0, 1) = u0; bcn(0, 0, 1) = uN
        !$omp end target
      END IF
      !$omp target teams distribute parallel do collapse(3) private(iPhi, ix, iz)
      DO iPhi = 1, nPhi
        DO ix = nx0, nxN
          DO iz = -nz, nz
            IF (ix == 0 .and. iz == 0) THEN
              bc0(iz, ix, 5 + iPhi) = t0
              bcn(iz, ix, 5 + iPhi) = tn
            ELSE
              bc0(iz, ix, 5 + iPhi) = 0
              bcn(iz, ix, 5 + iPhi) = 0
            END IF
          END DO
        END DO
      END DO
      !$omp end target teams distribute parallel do
      call roctxPop("boundary_conditions")
      ! Increment number of steps
      istep = istep + 1

      ! Loop over sub timestep
      do i = 1, 3
        call roctxPush("rk_substep")
        time = time + 2.0/RK_rai(1, i)*deltat
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
        call acquire_fft_workspace("rk_fft")
#endif
        call roctxPush("transform_to_physical")
        CALL transform_to_physical()
        call roctxPop("transform_to_physical")

        if (i .eq. 3) THEN
          call roctxPush("compute_cfl")
          call compute_cfl
          call roctxPop("compute_cfl")
        END IF

        !only depends on data from the previous substep, updates V(:, :, :, 1:2) and oldrhs(:, :, :, 1:2)
        !can be done in parallel to FFTs
        call roctxPush("buildrhs_prepare")
        CALL buildrhs_prepare(RK_rai(:, i))
        call roctxPop("buildrhs_prepare")

        call roctxPush("transform_back_and_build_rhs")
        CALL transform_back_and_build_rhs(RK_rai(:, i))
        call roctxPop("transform_back_and_build_rhs")
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
        call release_fft_workspace("rk_fft")
#endif

        !depends on V(:, :, :, 1:2), updates V(:, :, :, 1:3)
        call roctxPush("linsolve_velocity")
        CALL linsolve(RK_rai(1, i)/deltat)
        call roctxPop("linsolve_velocity")
        call roctxPush("linsolve_scalar")
        do iPhi = 1, nPhi
          !depends on (V(:, :, :, 3+iPhi), updates V(:, :, :, 3+iPhi)
          CALL linsolve_scalar(RK_rai(1, i)/deltat, iPhi)
        end do
        call roctxPop("linsolve_scalar")
        call roctxPop("rk_substep")
      end do

      call roctxPush("convvelo_runtime")
      call advance_convvelo_runtime()
      call roctxPop("convvelo_runtime")

      ! Write runtime file
      call roctxPush("outstats")
      CALL outstats()
      call roctxPop("outstats")

#ifdef chron
      CALL CPU_TIME(timee)
      elapsed_run_time = elapsed_run_time + (timee - timei)
      if (has_terminal) then
        write (*, '(A,I0,A,I0,A,F12.6,A,F12.6)') "STEP ", istep, "/", nstep, " TIME PER TIMESTEP ", timee - timei, &
          " ELAPSED RUN TIME ", elapsed_run_time
      end if
#endif
      call roctxPop("timestep")
    END DO
  END SUBROUTINE timeloop

  SUBROUTINE finalize()
    USE case_setup, only: disable_restart_write, has_terminal, time, save_restart_file, V, &
#ifdef HAVE_FFTW
                       VVdz, VVdx, rVVdx, &
#endif
                       free_memory
    USE convvelo, only: finalize_convvelo_runtime
    USE pressure_output, only: free_pressure_output
    USE byte_workspace, only: workspace_finalize
    USE y_schur_solver, only: ys_schur_finalize_contexts
    USE y_line_solvers, only: ys_finalize_nccl_contexts
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    USE ffts, only: free_fft
#elif defined(HAVE_FFTW)
    USE ffts, only: free_fft
#endif
    USE mpi_transpose, only: free_MPI, finalize_xcomm_nccl_contexts
    ! Without this, MPI_Finalize below resolves to the implicit external
    ! mpi_finalize_ -- the old binding, whose ierror argument is mandatory --
    ! and the zero-argument call makes it write through a garbage pointer.
    USE mpi_f08, only: MPI_Finalize
    IMPLICIT NONE
    if (disable_restart_write) then
      IF (has_terminal) WRITE (*, *) "End of time/iterations loop: restart write disabled for benchmark profiling at time ", time
    else
      IF (has_terminal) WRITE (*, *) "End of time/iterations loop: writing restart file at time ", time
      CALL save_restart_file("Dati.cart.out", V)
    end if

    IF (has_terminal) CLOSE (102)
    call finalize_convvelo_runtime()
    CALL free_pressure_output()
    ! Realease memory
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    CALL free_fft()
#elif defined(HAVE_FFTW)
    CALL free_fft()
#endif
    call ys_schur_finalize_contexts()
    call ys_finalize_nccl_contexts()
    call finalize_xcomm_nccl_contexts()
    CALL free_MPI()
    CALL free_memory(.TRUE.)
    call workspace_finalize()
    CALL MPI_Finalize()
  END SUBROUTINE finalize

  subroutine print_memory_line(label, n_floats)
    use, intrinsic :: iso_c_binding, only: C_INT64_T, C_DOUBLE
    implicit none
    character(len=*), intent(in) :: label
    integer(C_INT64_T), intent(in) :: n_floats
    real(C_DOUBLE) :: total_mib

    total_mib = floats_to_mib(n_floats)
    write (*, '(A,": device=",F12.3,A)') "  "//trim(label), total_mib, " MiB"
  end subroutine print_memory_line

  real(C_DOUBLE) function floats_to_mib(n_floats)
    use, intrinsic :: iso_c_binding, only: C_INT64_T, C_DOUBLE
    implicit none
    integer(C_INT64_T), intent(in) :: n_floats

    floats_to_mib = 8.0d0*real(n_floats, C_DOUBLE)/(1024.0d0*1024.0d0)
  end function floats_to_mib

END MODULE driver
