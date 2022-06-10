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

PROGRAM channel

  USE dnsdata
#ifdef runtimestats
  USE runtime
#endif
#ifdef crhon
  REAL timei,timee
#endif
  REAL(C_DOUBLE) :: frl(1:3), deltat_from_dnsin
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
#ifdef ibm
  CALL init_fft(VVdz,VVdx,rVVdx,Fdz,Fdx,rFdx,nxd,nxB,nzd,nzB)
#else
  CALL init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB)
#endif
  CALL setup_derivatives()
  CALL setup_boundary_conditions()
  fname="Dati.cart.out";       CALL read_restart_file(fname,V)
  ! move cursor to desired record
  IF (time_from_restart .AND. rtd_exists) THEN
    CALL get_record(time)
  ELSE IF (.NOT. time_from_restart) THEN
    CALL read_dnsin()
  END IF
#ifdef ibm
  fname="ibm.bin";             CALL read_body_file(fname,InBody(:,:,:,:))
  fname="dUint.cart.out";      CALL read_body_file(fname,dUint(:,:,:,:,0))
#endif

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
  WRITE(*,*) "!                     D   N   S                      !"
  WRITE(*,*) "!====================================================!"
  WRITE(*,*) " "
  WRITE(*,"(A,I5,A,I5,A,I5)") "   nx =",nx,"   ny =",ny,"   nz =",nz
  WRITE(*,"(A,I5,A,I5)") "   nxd =",nxd,"  nzd =",nzd
  WRITE(*,"(A,F11.6,A,F11.6,A,F8.6)") "   alfa0 =",alfa0,"       beta0 =",beta0,"   ni =",ni
  WRITE(*,"(A,F11.6,A,F11.6)") "   meanpx =",meanpx,"      meanpz =",meanpz
  WRITE(*,"(A,F11.6,A,F11.6)") "   meanflowx =",meanflowx, "   meanflowz =", meanflowz
  WRITE(*,"(A,I6,A,L1)"   ) "   nsteps =",nstep, "   time_from_restart =", time_from_restart
#ifdef ibm
  WRITE(*,"(A)") "   Immersed Boundary Method (IBM) active."
#endif
  WRITE(*,*) " "
END IF
#ifdef runtimestats
  call runtime_setup()
#endif
#ifdef bodyforce
  call config_body_force()
#endif

  ! Compute CFL
  DO iy=ny0,nyN
    IF (deltat==0) deltat=1.0; 
    CALL convolutions(iy,1,.TRUE.,.FALSE.,RK1_rai(1))
  END DO
  ! Compute flow rate flow rate
  IF (has_average) THEN
    frl(1)=yintegr(dreal(V(:,0,0,1))); frl(2)=yintegr(dreal(V(:,0,0,3))); 
    CALL MPI_Allreduce(frl,fr,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Y)
    IF (CPI) THEN
      SELECT CASE (CPI_type)
        CASE (0)
          meanpx = (1-gamma)*6*ni/fr(1)
        CASE (1)
          meanpx = (1.5d0/gamma)*fr(1)*ni 
        CASE DEFAULT
          WRITE(*,*) "Wrong selection of CPI_Type"
          STOP
      END SELECT
    END IF
  END IF
  CALL outstats()
  ! Time loop 
  timeloop: DO WHILE ((time<t_max-deltat/2.0) .AND. (istep<nstep))
#ifdef chron
    CALL CPU_TIME(timei)
#endif
    IF (has_average) THEN ! apply boundary conditions from dns.in (Couette-like)
      bc0(0,0)%u=u0; bcn(0,0)%u=uN
    END IF
    ! Increment number of steps
    istep=istep+1
    ! Solve (RK3 - Step1)
    time=time+2.0/RK1_rai(1)*deltat
#ifdef bodyforce
    CALL set_body_force(); 
#endif
    CALL buildrhs(RK1_rai,.FALSE. )
#ifdef nonblockingY
    CALL MPI_Barrier(MPI_COMM_Y)
#endif
    CALL linsolve(RK1_rai(1)/deltat)
#ifdef nonblockingY
    CALL vetaTOuvw(); CALL computeflowrate(RK1_rai(1)/deltat)
#endif
    ! Solve (RK3 - Step2)
    time=time+2.0/RK2_rai(1)*deltat
#ifdef bodyforce
    CALL set_body_force(); 
#endif
    CALL buildrhs(RK2_rai,.FALSE.)  
#ifdef nonblockingY
    CALL MPI_Barrier(MPI_COMM_Y)
#endif
    CALL linsolve(RK2_rai(1)/deltat)
#ifdef nonblockingY
    CALL vetaTOuvw(); CALL computeflowrate(RK2_rai(1)/deltat)
#endif
    ! Solve (RK3 - Step3)
    time=time+2.0/RK3_rai(1)*deltat
#ifdef bodyforce
    CALL set_body_force(); 
#endif
    CALL buildrhs(RK3_rai,.TRUE.)
#ifdef nonblockingY
    CALL MPI_Barrier(MPI_COMM_Y)
#endif
    CALL linsolve(RK3_rai(1)/deltat)
#ifdef nonblockingY
    CALL vetaTOuvw(); CALL computeflowrate(RK3_rai(1)/deltat)
#endif
    ! Write runtime file
    CALL outstats()
#ifdef runtimestats
    IF (rtstats_savenow) THEN
      IF (has_terminal) print *, "Saving runtime stats at time", time
      CALL runtime_save()
      rtstats_savenow = .FALSE.
    END IF
#endif
#ifdef chron
    CALL CPU_TIME(timee)
    IF (has_terminal) WRITE(*,*) timee-timei
#endif
  END DO timeloop
  IF (has_terminal) WRITE(*,*) "End of time/iterations loop: writing restart file at time ", time
  end_filename="Dati.cart.out";  CALL save_restart_file(end_filename,V)
#ifdef ibm
  IF (has_terminal) WRITE(*,*) "End of time/iterations loop: writing dUint.cart.out at time ", time
  filename="dUint.cart.out"; CALL save_body_file(filename,dUint(:,:,:,:,0))
#endif
  IF (has_terminal) CLOSE(102)
  ! Realease memory
#ifdef runtimestats
  call runtime_finalise()
#endif
  CALL free_fft()
  CALL free_memory(.TRUE.) 
  CALL MPI_Finalize()


END PROGRAM channel
