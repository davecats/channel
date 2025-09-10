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
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: config_file, restart_file
    REAL(C_DOUBLE) :: deltat_from_dnsin
    integer :: iy

    ! Init MPI
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

    CALL read_dnsin(config_file)
    deltat_from_dnsin = deltat
    CALL init_MPI(nx + 1, nz, ny, nxd + 1, nzd)
    CALL init_memory(.TRUE.)

    ! Init various subroutines
    CALL init_fft(VVdz, VVdx, rVVdx, nxd, nxB, ny, nzd, nzB)
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
    DO iy = ny0, nyN
      IF (deltat == 0) deltat = 1.0; 
      CALL convolutions(iy, .TRUE.)
    END DO
    ! Compute flow rate
    IF (has_average) THEN
      fr(1) = yintegr(dreal(V(:, 0, 0, 1)), y); fr(2) = yintegr(dreal(V(:, 0, 0, 3)), y); 
    END IF
    CALL outstats()
  END SUBROUTINE initialize

  !==========================================================
  SUBROUTINE timeloop()
    USE dnsdata
    IMPLICIT NONE
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
      ! Increment number of steps
      istep = istep + 1

      ! Solve (RK3 - Step1)
      time = time + 2.0/RK1_rai(1)*deltat
      CALL buildrhs(RK1_rai, .FALSE.)
      CALL linsolve(RK1_rai(1)/deltat)

      ! Solve (RK3 - Step2)
      time = time + 2.0/RK2_rai(1)*deltat
      CALL buildrhs(RK2_rai, .FALSE.)
      CALL linsolve(RK2_rai(1)/deltat)

      ! Solve (RK3 - Step3)
      time = time + 2.0/RK3_rai(1)*deltat
      CALL buildrhs(RK3_rai, .TRUE.)
      CALL linsolve(RK3_rai(1)/deltat)

      ! Write runtime file
      CALL outstats()

#ifdef chron
      CALL CPU_TIME(timee)
      IF (has_terminal) WRITE (*, *) timee - timei
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
    CALL free_fft(VVdz, VVdx, rVVdx)
    CALL free_memory(.TRUE.)
    CALL MPI_Finalize()
  END SUBROUTINE finalize
END MODULE driver
