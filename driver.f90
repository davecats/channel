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
  SUBROUTINE initialize(config_file, restart_file)
    USE dnsdata
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    use omp_lib
#endif
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: config_file, restart_file
    REAL(C_DOUBLE) :: deltat_from_dnsin
    integer :: iy, iPhi, num_dev, dev

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

    CALL read_dnsin(config_file)
    deltat_from_dnsin = deltat
    CALL init_MPI(nx + 1, nz, ny, nxd + 1, nzd, nPhi)
    CALL init_memory(.TRUE.)

    ! Init various subroutines
#ifdef HAVE_CUDA
    CALL init_cufft(nxd, nxB, ny, nzd, nzB, nPhi)
#elif defined(HAVE_HIP)
    CALL init_hipfft(nxd, nxB, ny, nzd, nzB, nPhi)
#elif defined(HAVE_FFTW)
    CALL init_fft(VVdz, VVdx, rVVdx, nxd, nxB, ny, nzd, nzB, nPhi)
#endif
    CALL setup_derivatives()
    CALL setup_boundary_conditions()
    CALL read_restart_file(restart_file, V)

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
      WRITE (*, *) " "
    END IF

    ! Compute CFL
    if (deltat == 0.0) deltat = 1.0
    !$omp target update to(V)
    CALL transform_to_physical()
    call compute_cfl()
    print *, "CFL", deltat, cfl
    ! Compute flow rate
    IF (has_average) THEN
      fr(1) = yintegr(V(:, 0, 0, 1), y); fr(2) = yintegr(V(:, 0, 0, 3), y); 
      DO iPhi = 1, nPhi
        fr(3 + iPhi) = yintegr(V(:, 0, 0, 3 + iPhi), y)
      END DO
    END IF
    CALL outstats()
  END SUBROUTINE initialize

  !==========================================================
  SUBROUTINE timeloop()
    USE dnsdata
    IMPLICIT NONE
    integer:: iPhi, ix, iz, i, ic
#ifdef chron
    REAL timei, timee
#endif

    DO WHILE ((time < t_max - deltat/2.0) .AND. (istep < nstep))
#ifdef chron
      CALL CPU_TIME(timei)
#endif
      ! apply boundary conditions from input file (Couette-like)
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
      ! Increment number of steps
      istep = istep + 1

      ! Loop over sub timestep
      do i = 1, 3
        time = time + 2.0/RK_rai(1, i)*deltat
        CALL transform_to_physical()

        if (i .eq. 3) THEN
          call compute_cfl
        END IF

        !only depends on data from the previous substep, updates V(:, :, :, 1:2) and oldrhs(:, :, :, 1:2)
        !can be done in parallel to FFTs
        CALL buildrhs_prepare(RK_rai(:, i))

        CALL transform_back_and_build_rhs(RK_rai(:, i))

        !depends on V(:, :, :, 1:2), updates V(:, :, :, 1:3)
        CALL linsolve(RK_rai(1, i)/deltat)
        do iPhi = 1, nPhi
          !depends on (V(:, :, :, 3+iPhi), updates V(:, :, :, 3+iPhi)
          CALL linsolve_scalar(RK_rai(1, i)/deltat, iPhi)
        end do
      end do

      ! Write runtime file
      CALL outstats()

#ifdef chron
      CALL CPU_TIME(timee)
      IF (has_terminal) WRITE (*, *) "TIME PER TIMESTEP", timee - timei
#endif
    END DO
  END SUBROUTINE timeloop

  SUBROUTINE finalize()
    USE dnsdata
    IMPLICIT NONE
    CHARACTER(len=40) :: end_filename

    IF (has_terminal) WRITE (*, *) "End of time/iterations loop: writing restart file at time ", time
    end_filename = "Dati.cart.out"; CALL save_restart_file(end_filename, V)

    IF (has_terminal) CLOSE (102)
    ! Realease memory
#ifdef HAVE_FFTW
    CALL free_fft(VVdz, VVdx, rVVdx)
#endif
    CALL free_memory(.TRUE.)
#ifdef HAVE_MPI
    CALL MPI_Finalize()
#endif
  END SUBROUTINE finalize
END MODULE driver
