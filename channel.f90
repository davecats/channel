!============================================!
!                                            !
!     Direct Numerical Simulation (DNS)      !
!       of a turbulent channel flow          !
!         (OpenMP GPU-accelerated)           !
!                                            !
!============================================!
!
! This program has been written following the
! KISS (Keep it Simple and Stupid) philosophy
!
! Author: Dr.-Ing. Davide Gatti
! GPU Adaptation: Gemini
!

#include "header.h"

PROGRAM channel

  USE dnsdata
#ifdef crhon
  REAL timei,timee
#endif
  REAL(C_DOUBLE) :: deltat_from_dnsin
  CHARACTER(len = 40) :: end_filename

  ! Init MPI
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  CALL read_dnsin()
  deltat_from_dnsin = deltat
  CALL init_MPI(nx+1,nz,ny,nxd+1,nzd)
  CALL init_memory(.TRUE.)

  ! Init various subroutines
  ! Use cuFFT for GPU
  CALL init_fft(nxd,nxB,nzd,nzB)
  CALL setup_derivatives()
  CALL setup_boundary_conditions()
  fname="Dati.cart.out"
  CALL read_restart_file(fname,V)

  ! Field number (for output)
  ifield=FLOOR((time+0.5*deltat)/dt_field)
  time0=time

  ! Reset desired dt
  IF (cflmax == 0) THEN
    deltat = deltat_from_dnsin
  END IF

IF (has_terminal) THEN
  ! Output DNS.in
  WRITE(*,*) " "
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) "!              D   N   S  (GPU Version)               !"
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) " "
  WRITE(*,"(A,I5,A,I5,A,I5)") "   nx =",nx,"   ny =",ny,"   nz =",nz
  WRITE(*,"(A,I5,A,I5)") "   nxd =",nxd,"  nzd =",nzd
  WRITE(*,"(A,F11.6,A,F11.6,A,F8.6)") "   alfa0 =",alfa0,"       beta0 =",beta0,"   ni =",ni
  WRITE(*,"(A,F11.6,A,F11.6)") "   meanpx =",meanpx,"      meanpz =",meanpz
  WRITE(*,"(A,F11.6,A,F11.6)") "   meanflowx =",meanflowx, "   meanflowz =", meanflowz
  WRITE(*,"(A,I6,A,L1)"   ) "   nsteps =",nstep, "   time_from_restart =", time_from_restart
  WRITE(*,*) " "
END IF

!$OMP TARGET DATA MAP(TOFROM: V, der) &
!$OMP               MAP(ALLOC: memrhs, oldrhs, bc0, bcn) &
!$OMP               MAP(ALLOC: D0mat, etamat, eta00mat, D2vmat) &
!$OMP               MAP(TO: y, dy, izd, ialfa, ibeta, k2) &
!$OMP               MAP(TO: d040,d140,d240,d14m1,d24m1,d04n,d14n,d24n,d14np1,d24np1) &
!$OMP               MAP(TO: v0bc,v0m1bc,vnbc,vnp1bc,eta0bc,eta0m1bc,etanbc,etanp1bc) &
!$OMP               MAP(TO: RK1_rai, RK2_rai, RK3_rai) &
!$OMP               MAP(ALLOC: VVdz, VVdx, rVVdx)

  ! Manually update device with initial data
!$OMP TARGET UPDATE TO(V, der, y, dy, izd, ialfa, ibeta, k2)
!$OMP TARGET UPDATE TO(d040,d140,d240,d14m1,d24m1,d04n,d14n,d24n,d14np1,d24np1)
!$OMP TARGET UPDATE TO(v0bc,v0m1bc,vnbc,vnp1bc,eta0bc,eta0m1bc,etanbc,etanp1bc)
!$OMP TARGET UPDATE TO(RK1_rai, RK2_rai, RK3_rai)

  ! Compute CFL
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
  DO iy=ny0,nyN
    IF (deltat==0) deltat=1.0;
    CALL convolutions(iy,1,.TRUE.)
  END DO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

!$OMP TARGET UPDATE FROM(cfl)

  ! Compute flow rate
  IF (has_average) THEN
!$OMP TARGET MAP(TOFROM: fr)
    fr(1)=yintegr(dreal(V(:,0,0,1))); fr(2)=yintegr(dreal(V(:,0,0,3)));
!$OMP END TARGET
  END IF

  CALL outstats()
  ! Time loop
  timeloop: DO WHILE ((time<t_max-deltat/2.0) .AND. (istep<nstep))
#ifdef chron
    CALL CPU_TIME(timei)
#endif
    ! apply boundary conditions from dns.in (Couette-like)
    IF (has_average) THEN
!$OMP TARGET
      bc0(0,0,1)=u0; bcn(0,0,1)=uN
!$OMP END TARGET
    END IF
    ! Increment number of steps
    istep=istep+1

    ! Solve (RK3 - Step1)
    time=time+2.0/RK1_rai(1)*deltat
!$OMP TARGET NOWAIT
    CALL buildrhs(RK1_rai,.FALSE. )
!$OMP END TARGET
!$OMP TARGET NOWAIT
    CALL linsolve(RK1_rai(1)/deltat)
!$OMP END TARGET

    ! Solve (RK3 - Step2)
    time=time+2.0/RK2_rai(1)*deltat
!$OMP TARGET NOWAIT
    CALL buildrhs(RK2_rai,.FALSE.)
!$OMP END TARGET
!$OMP TARGET NOWAIT
    CALL linsolve(RK2_rai(1)/deltat)
!$OMP END TARGET

    ! Solve (RK3 - Step3)
    time=time+2.0/RK3_rai(1)*deltat
!$OMP TARGET NOWAIT
    CALL buildrhs(RK3_rai,.TRUE.)
!$OMP END TARGET
!$OMP TARGET NOWAIT
    CALL linsolve(RK3_rai(1)/deltat)
!$OMP END TARGET
!$OMP TASKWAIT

    ! Write runtime file
    CALL outstats()

#ifdef chron
    CALL CPU_TIME(timee)
    IF (has_terminal) WRITE(*,*) timee-timei
#endif
  END DO timeloop

!$OMP END TARGET DATA

  IF (has_terminal) WRITE(*,*) "End of time/iterations loop: writing restart file at time ", time
  end_filename="Dati.cart.out";  CALL save_restart_file(end_filename,V)

  IF (has_terminal) CLOSE(102)
  ! Realease memory
  CALL free_fft()
  CALL free_memory(.TRUE.)
  CALL MPI_Finalize()


END PROGRAM channel

