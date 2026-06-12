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

#include "header.h"

MODULE driver

CONTAINS
  !==========================================================
  SUBROUTINE initialize(config_file, restart_file, solveNS)
    use config, only: ini_config, read_ini_file
    USE dnsdata
    USE convvelo, only: init_convvelo_runtime, get_convvelo_memory_estimate, configure_convvelo
    USE ffts, only: get_fft_memory_estimate
#ifdef HAVE_CUDA
    USE ffts, only: init_cufft
#elif defined(HAVE_HIP)
    USE ffts, only: init_hipfft
#else
    USE ffts, only: init_fft, free_fft
#endif
    USE pressure_output, only: init_pressure_output, free_pressure_output, get_pressure_memory_estimate
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    use omp_lib
#endif
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: config_file, restart_file
    LOGICAL, OPTIONAL, INTENT(IN) :: solveNS
    type(ini_config) :: cfg
    REAL(C_DOUBLE) :: deltat_from_dnsin
    real(C_DOUBLE) :: total_mib
    integer(C_INT64_T) :: solver_floats, fft_floats, pressure_floats, convvelo_floats
    integer :: iy, iPhi, num_dev, dev
    logical :: run_solver
    complex(C_DOUBLE_COMPLEX), allocatable :: zero_mode(:)

    run_solver = .true.
    if (present(solveNS)) run_solver = solveNS

    ! Init MPI
#ifdef HAVE_MPI
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
#else
    iproc = 0
    nproc = 1
#endif

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    num_dev = omp_get_num_devices()
    dev = mod(iproc, num_dev)

    call omp_set_default_device(dev)

    print *, 'Rank', iproc, 'of', nproc, 'using device', dev, 'out of', num_dev

    !$omp target
    print *, 'Hello from GPU on rank', iproc, 'device', dev
    !$omp end target
#endif

    call read_ini_file(config_file, cfg)
    CALL read_dnsin(cfg)
    call configure_convvelo(cfg)
    deltat_from_dnsin = deltat
    CALL init_MPI(nx + 1, nz, ny, nzd, nPhi, overlapping, npy)
    call get_solver_memory_estimate(run_solver, solver_floats)
    call get_fft_memory_estimate(nxd, nxB, nzd, nzB, nPhi, overlapping, fft_floats)
    call get_pressure_memory_estimate(pressure_floats)
    call get_convvelo_memory_estimate(convvelo_floats)
    if (has_terminal) then
      write (*, *) "Estimated memory per rank before allocation:"
      call print_memory_line("Solver", solver_floats)
      call print_memory_line("FFT", fft_floats)
      call print_memory_line("Pressure", pressure_floats)
      call print_memory_line("Convvelo", convvelo_floats)
      total_mib = floats_to_mib(solver_floats + fft_floats + pressure_floats + convvelo_floats)
      write (*, '(A,F12.3,A)') "  Total      : device=", total_mib, " MiB"
      write (*, *) " "
    end if
    CALL init_memory(run_solver)

    ! Init various subroutines
#ifdef HAVE_CUDA
    CALL init_cufft(nxd, nxB, nzd, nzB, nPhi, overlapping)
#elif defined(HAVE_HIP)
    CALL init_hipfft(nxd, nxB, nzd, nzB, nPhi, overlapping)
#elif defined(HAVE_FFTW)
    CALL init_fft(VVdz, VVdx, rVVdx, nxd, nxB, nzd, nzB, nPhi, overlapping)
#endif
    CALL setup_derivatives()
    CALL setup_boundary_conditions()
    CALL read_restart_file(restart_file, V)
    if (.not. run_solver) then
      !$omp target update to(V)
    end if
    CALL init_pressure_output()
    call init_convvelo_runtime()

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
      CALL transform_to_physical()
      call compute_cfl()
      print *, "CFL", deltat, cfl
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

  !==========================================================
  SUBROUTINE timeloop()
    USE dnsdata
    USE convvelo, only: advance_convvelo_runtime
    IMPLICIT NONE
    integer:: iPhi, ix, iz, i, ic
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
        bc0(0, 0, 1) = u0; bcn(0, 0, 1) = uN
      END IF
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
      !$omp target update to(bc0, bcn)
      call roctxPop("boundary_conditions")
      ! Increment number of steps
      istep = istep + 1

      ! Loop over sub timestep
      do i = 1, 3
        call roctxPush("rk_substep")
        time = time + 2.0/RK_rai(1, i)*deltat
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
    USE dnsdata, only: disable_restart_write, has_terminal, time, save_restart_file, V, &
#ifdef HAVE_FFTW
                       VVdz, VVdx, rVVdx, &
#endif
                       free_memory
    USE convvelo, only: finalize_convvelo_runtime
    USE pressure_output, only: free_pressure_output
#ifdef HAVE_FFTW
    USE ffts, only: free_fft
#endif
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
#ifdef HAVE_FFTW
    CALL free_fft(VVdz, VVdx, rVVdx)
#endif
    CALL free_memory(.TRUE.)
#ifdef HAVE_MPI
    CALL MPI_Finalize()
#endif
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
