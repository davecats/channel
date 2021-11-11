!============================================!
!                                            !
!    Data Structures, Definitions and I/O    !
!                  for the                   !
!     Direct Numerical Simulation (DNS)      !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: Dr.-Ing. Davide Gatti
! Date  : 28/Jul/2015
!

#include "header.h"

MODULE dnsdata

  USE, intrinsic :: iso_c_binding
  USE rbmat
  USE mpi_transpose
  USE ffts

  IMPLICIT NONE

  !Simulation parameters
  real(C_DOUBLE) :: PI=3.1415926535897932384626433832795028841971
  integer(C_INT) :: nx,ny,nz,nxd,nzd
  real(C_DOUBLE) :: alfa0,beta0,ni,a,ymin,ymax,deltat,cflmax,time,time0=0,dt_field,dt_save,t_max,gamma
  real(C_DOUBLE) :: u0,uN
  logical :: CPI
  integer(C_INT) :: CPI_type
  real(C_DOUBLE) :: meanpx,meanpz,meanflowx,meanflowz
  integer(C_INT), allocatable :: izd(:)
  complex(C_DOUBLE_COMPLEX), allocatable :: ialfa(:),ibeta(:)
  real(C_DOUBLE), allocatable :: k2(:,:)
  logical :: time_from_restart
  !Grid
  integer(C_INT), private :: iy
  real(C_DOUBLE), allocatable :: y(:),dy(:)
  real(C_DOUBLE) :: dx,dz,factor
  !Derivatives
  TYPE(Di), allocatable :: der(:)
  real(C_DOUBLE), dimension(-2:2) :: d040,d140,d240,d14m1,d24m1,d04n,d14n,d24n,d14np1,d24np1
  real(C_DOUBLE), allocatable :: D0mat(:,:), etamat(:,:,:), eta00mat(:,:), D2vmat(:,:,:)
  !Fourier-transformable arrays (allocated in ffts.f90)
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:) :: VVdx, VVdz
  real(C_DOUBLE), pointer, dimension(:,:,:,:) :: rVVdx
  !Solution
  TYPE(RHSTYPE),  allocatable :: memrhs(:,:,:), oldrhs(:,:,:)
  complex(C_DOUBLE_COMPLEX), allocatable :: V(:,:,:,:)
#ifdef bodyforce
  complex(C_DOUBLE_COMPLEX), allocatable :: F(:,:,:,:)
#ifdef ibm
  real(C_DOUBLE) :: ibmp=5000, ibmi=10000
  real(C_DOUBLE), allocatable :: dUint(:,:,:,:,:) 
  real(C_DOUBLE), allocatable :: Ut(:,:,:,:) 
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:) :: Fdx, Fdz   ! Fourier-transformable (ffts.f90)
  real(C_DOUBLE), pointer, dimension(:,:,:,:) :: rFdx                  ! Fourier-transformable (ffts.f90)
  real(C_DOUBLE), allocatable, dimension(:,:,:,:) :: InBody          
#endif
#endif
  !Boundary conditions
  real(C_DOUBLE), dimension(-2:2) :: v0bc,v0m1bc,vnbc,vnp1bc,eta0bc,eta0m1bc,etanbc,etanp1bc
  TYPE(BCOND),    allocatable :: bc0(:,:), bcn(:,:)
  !Mean pressure correction
  real(C_DOUBLE), private :: corrpx=0.d0, corrpz=0.d0
  !ODE Library
  real(C_DOUBLE) :: RK1_rai(1:3)=(/ 120.0d0/32.0d0, 2.0d0, 0.0d0 /), &
                    RK2_rai(1:3)=(/ 120.0d0/8.0d0,  50.0d0/8.0d0,  34.0d0/8.0d0 /), &
                    RK3_rai(1:3)=(/ 120.0d0/20.0d0, 90.0d0/20.0d0, 50.0d0/20.0d0 /)
  !Outstats
  real(C_DOUBLE) :: cfl=0.0d0
  integer(C_SIZE_T) :: istep,nstep,ifield
  real(C_DOUBLE) :: fr(1:3)
  logical :: prev_was_close = .FALSE.
  !Restart file
  character(len=40) :: fname
  !Convection velocity calculation
#ifdef convvel
  complex(C_DOUBLE_COMPLEX), allocatable :: Voldz(:,:,:,:)
  real(C_DOUBLE), allocatable :: uconv(:,:,:,:)
  integer(C_LONG) :: convvel_cnt=-1
  logical, save :: compute_convvel=.FALSE.
#endif
logical::rtd_exists ! flag to check existence of Runtimedata


  CONTAINS

  !--------------------------------------------------------------!
  !---------------------- Read input files ----------------------!
  SUBROUTINE read_dnsin(in_in_parent)
    logical, optional, intent(in) :: in_in_parent
    logical :: i
    if (present(in_in_parent)) then
      if (in_in_parent) then
        OPEN(15, file='../dns.in')
      else
        OPEN(15, file='dns.in')
      end if
    else
      OPEN(15, file='dns.in')
    end if
    READ(15, *) nx, ny, nz; READ(15, *) alfa0, beta0; nxd=3*(nx+1)/2;nzd=3*nz
#ifdef useFFTfit
    i=fftFIT(nxd); DO WHILE (.NOT. i); nxd=nxd+1; i=fftFIT(nxd); END DO
    i=fftFIT(nzd); DO WHILE (.NOT. i); nzd=nzd+1; i=fftFIT(nzd); END DO
#endif
    READ(15, *) ni; READ(15, *) a, ymin, ymax; ni=1/ni
    READ(15, *) CPI, CPI_type, gamma
    READ(15, *) meanpx, meanpz; READ(15, *) meanflowx, meanflowz
    READ(15, *) u0,uN
    READ(15, *) deltat, cflmax, time
    READ(15, *) dt_field, dt_save, t_max, time_from_restart
    READ(15, *) nstep
    READ(15, *) npy
    CLOSE(15)
    dx=PI/(alfa0*nxd); dz=2.0d0*PI/(beta0*nzd);  factor=1.0d0/(2.0d0*nxd*nzd)
  END SUBROUTINE read_dnsin

  !--------------------------------------------------------------!
  !---------------- Allocate memory for solution ----------------!
  SUBROUTINE init_memory(solveNS)
    INTEGER(C_INT) :: ix,iz
    logical, intent(IN) :: solveNS
    ALLOCATE(V(ny0-2:nyN+2,-nz:nz,nx0:nxN,1:3)); V=0
#ifdef bodyforce
    ALLOCATE(F(ny0-2:nyN+2,-nz:nz,nx0:nxN,1:3)); F=0
#ifdef ibm
    ALLOCATE(dUint(1:2*nxd,1:nzB,ny0-2:nyN+2,1:3,-1:0)); dUint=0
    ALLOCATE(Ut(1:2*nxd,1:nzB,ny0-2:nyN+2,1:3)); Ut=0
    ALLOCATE(InBody(1:2*nxd,1:nzB,ny0-2:nyN+2,1:1)); 
#endif
#endif
#ifdef convvel
    ALLOCATE(Voldz(1:nzd,1:nxB,ny0-2:nyN+2,1:3), uconv(1:nzd,1:nxB,ny0-2:nyN+2,1:3)); Voldz=0; uconv=0
#endif
    IF (solveNS) ALLOCATE(memrhs(0:2,-nz:nz,nx0:nxN),oldrhs(MAX(1,ny0-2):MIN(ny-1,nyN+2),-nz:nz,nx0:nxN),bc0(-nz:nz,nx0:nxN),bcn(-nz:nz,nx0:nxN))
#define newrhs(iy,iz,ix) memrhs(MOD(iy+1000,3),iz,ix)
#define imod(iy) MOD(iy+1000,5)
    ALLOCATE(der(MAX(1,ny0-2):MIN(ny-1,nyN+2)),d0mat(ny0:nyN+2,-2:2),etamat(ny0:nyN+2,-2:2,-nz:nz),eta00mat(ny0:nyN+2,-2:2),D2vmat(ny0:nyN+2,-2:2,-nz:nz))
    ALLOCATE(y(MAX(-1,ny0-4):MIN(ny+1,nyN+4)),dy(MAX(1,ny0-2):MIN(ny-1,nyN+2)))
    ALLOCATE(izd(-nz:nz),ialfa(nx0:nxN),ibeta(-nz:nz),k2(-nz:nz,nx0:nxN)) 
#ifdef halfchannel
    FORALL (iy=MAX(-1,ny0-4):MIN(ny+1,nyN+4)) y(iy)=ymin+(ymax-ymin)*(tanh(a*(real(iy)/real(ny)-1))/tanh(a)+1)
#else
    FORALL (iy=MAX(-1,ny0-4):MIN(ny+1,nyN+4)) y(iy)=ymin+0.5d0*(ymax-ymin)*(tanh(a*(2.0d0*real(iy)/real(ny)-1))/tanh(a)+1)
#endif
    FORALL (iy=MAX(1,ny0-2):MIN(ny-1,nyN+2)) dy(iy)=0.5d0*(y(iy+1)-y(iy-1))
    izd=(/(merge(iz,nzd+iz,iz>=0),iz=-nz,nz)/);     ialfa=(/(dcmplx(0.0d0,ix*alfa0),ix=nx0,nxN)/);
    ibeta=(/(dcmplx(0.0d0,iz*beta0),iz=-nz,nz)/); 
    FORALL  (iz=-nz:nz,ix=nx0:nxN) k2(iz,ix)=(alfa0*ix)**2.0d0+(beta0*iz)**2.0d0
    INQUIRE(FILE="Runtimedata", EXIST=rtd_exists)
    IF (solveNS .AND. has_terminal) THEN
      IF (time_from_restart .AND. rtd_exists) THEN
        WRITE(*,*) 'Found existing Runtimedata...'
        OPEN(UNIT=101,FILE='Runtimedata',ACTION='readwrite')
      ELSE
        IF (time_from_restart) THEN
          WRITE(*,*) 'Runtimedata not found...'
#ifdef warnings_are_fatal
          WRITE(*,*) 'Stopping!'
          STOP
#endif
          CALL wanna_continue()
        END IF
        WRITE(*,*) 'Creating new Runtimedata.'
        OPEN(UNIT=101,FILE='Runtimedata',ACTION='write')
      END IF
    END IF
  END SUBROUTINE init_memory

  !-------------------------------------------------------------------------------------!
  !--------------- Move cursor to correct instant in time in Runtimedata ---------------!
  SUBROUTINE get_record(threshold)
  IMPLICIT NONE
    REAL(C_DOUBLE), INTENT(IN) :: threshold
    REAL(C_DOUBLE) :: selectime,a,b,c,d,e,f,g,h,i,curr_dt
    INTEGER :: negative_if_eof = 0 
    LOGICAL :: threshold_reached
    IF (has_terminal) THEN
      DO WHILE (.NOT. threshold_reached .AND. negative_if_eof >= 0)
        READ(101,*,IOSTAT=negative_if_eof) selectime,a,b,c,d,e,f,g,h,i,curr_dt
        threshold_reached = ABS(selectime - threshold) < (0.5*curr_dt)
      END DO
      IF (negative_if_eof >= 0) THEN
        BACKSPACE(101)
        PRINT *, 'In Runtimedata: starting from time', selectime
        deltat = curr_dt ! this is only executed by machine WITH TERMINAL!
      ELSE
        PRINT *, ''
        PRINT *, '###############'
        PRINT *, '#   WARNING   #'
        PRINT *, '###############'
        PRINT *, 'No instant of time matching restart file has been found in Runtimedata.'
        PRINT *, ''
        WRITE(101,*) ''
#ifdef warnings_are_fatal
        WRITE(*,*) 'Stopping!'
        STOP
#endif
        CALL wanna_continue()
        PRINT *, 'Skipping one line and appending.'
      END IF
    END IF

    ! Broadcast deltat of machine with terminal to the remaining ones.
    ! This is not always needed, but it is if the deltat = curr_dt
    ! instruction above has been executed.
    CALL MPI_bcast(deltat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD) ! iproc == 0 has terminal

  END SUBROUTINE get_record

  !--------------------------------------------------------------!
  !--------------- Deallocate memory for solution ---------------!
  SUBROUTINE free_memory(solveNS)
    LOGICAL, intent(IN) :: solveNS
    DEALLOCATE(V,der,d0mat,etamat,D2vmat,y,dy)
    IF (solveNS) THEN
      DEALLOCATE(memrhs,oldrhs,bc0,bcn)
      IF (has_terminal) CLOSE(UNIT=101)
    END IF
  END SUBROUTINE free_memory

  !--------------------------------------------------------------!
  !--------------- Set-up the compact derivatives ---------------!
  SUBROUTINE setup_derivatives()
    real(C_DOUBLE)    :: M(0:4,0:4), t(0:4)
    integer(C_INT)    :: iy,i,j
    TYPE(MPI_REQUEST) :: Rs 
    DO iy=MAX(1,ny0-2),MIN(ny-1,nyN+2)
      FORALL (i=0:4, j=0:4) M(i,j)=(y(iy-2+j)-y(iy))**(4.0d0-i); CALL LUdecomp(M)
      t=0; t(0)=24
      der(iy)%d4(-2:2)=M.bs.t
      FORALL (i=0:4, j=0:4) M(i,j)=(5.0d0-i)*(6.0d0-i)*(7.0d0-i)*(8.0d0-i)*(y(iy-2+j)-y(iy))**(4.0d0-i); CALL LUdecomp(M)
      FORALL (i=0:4) t(i)=sum( der(iy)%d4(-2:2)*(y(iy-2:iy+2)-y(iy))**(8.0d0-i) )
      der(iy)%d0(-2:2)=M.bs.t
      FORALL (i=0:4, j=0:4) M(i,j)=(y(iy-2+j)-y(iy))**(4.0d0-i); CALL LUdecomp(M)
      t=0; FORALL (i=0:2) t(i)=sum( der(iy)%d0(-2:2)*(4.0d0-i)*(3.0d0-i)*(y(iy-2:iy+2)-y(iy))**(2.0d0-i) )
      der(iy)%d2(-2:2)=M.bs.t
      t=0; FORALL (i=0:3) t(i)=sum( der(iy)%d0(-2:2)*(4.0d0-i)*(y(iy-2:iy+2)-y(iy))**(3.0d0-i) )
      der(iy)%d1(-2:2)=M.bs.t
    END DO
    IF (first) THEN
      FORALL (i=0:4, j=0:4) M(i,j)=(y(-1+j)-y(0))**(4.0d0-i); CALL LUdecomp(M)
      t=0; t(3)=1.0; d140(-2:2)=M.bs.t
      t=0; t(2)=2.0; d240(-2:2)=M.bs.t
      FORALL (i=0:4, j=0:4) M(i,j)=(y(-1+j)-y(-1))**(4.0d0-i); CALL LUdecomp(M)
      t=0; t(3)=1.0; d14m1(-2:2)=M.bs.t
      t=0; t(2)=2.0; d24m1(-2:2)=M.bs.t
      d040=0; d040(-1)=1
    END IF
    IF (last) THEN 
      FORALL (i=0:4, j=0:4) M(i,j)=(y(ny-3+j)-y(ny))**(4.0d0-i); CALL LUdecomp(M)
      t=0; t(3)=1; d14n(-2:2)=M.bs.t
      t=0; t(2)=2; d24n(-2:2)=M.bs.t
      FORALL (i=0:4, j=0:4) M(i,j)=(y(ny-3+j)-y(ny+1))**(4.0d0-i); CALL LUdecomp(M)
      t=0; t(3)=1; d14np1(-2:2)=M.bs.t
      t=0; t(2)=2; d24np1(-2:2)=M.bs.t
      d04n=0; d04n(1)=1;
    END IF
    FORALL (iy=ny0:nyN) D0mat(iy,-2:2)=der(iy)%d0(-2:2); 
#ifdef nonblockingY
    CALL LU5decompStep(D0mat,Rs,0)
    IF ( .NOT. first ) CALL MPI_Wait(Rs,MPI_STATUS_IGNORE)
#else
    CALL LU5decompStep(D0mat)
#endif
  END SUBROUTINE setup_derivatives

  !--------------------------------------------------------------!
  !--------------- Set-up the boundary conditions ---------------!
  SUBROUTINE setup_boundary_conditions()
    INTEGER :: i
    IF (first) THEN 
      v0bc=d040; v0m1bc=d140; eta0bc=d040
      eta0m1bc=der(1)%d4
      v0bc(-1:2)=v0bc(-1:2)-v0bc(-2)*v0m1bc(-1:2)/v0m1bc(-2)
      eta0bc(-1:2)=eta0bc(-1:2)-eta0bc(-2)*eta0m1bc(-1:2)/eta0m1bc(-2)
    END IF
    IF (last) THEN 
#ifdef halfchannel
      vnbc=d04n; vnp1bc=d24n; etanbc=d14n
#else
      vnbc=d04n; vnp1bc=d14n; etanbc=d04n
#endif
      etanp1bc=der(ny-1)%d4
      vnbc(-2:1)=vnbc(-2:1)-vnbc(2)*vnp1bc(-2:1)/vnp1bc(2)
      etanbc(-2:1)=etanbc(-2:1)-etanbc(2)*etanp1bc(-2:1)/etanp1bc(2)
    END IF 
  END SUBROUTINE setup_boundary_conditions

  !--------------------------------------------------------------!
  !---------------- integral in the y-direction -----------------!
  PURE FUNCTION yintegr(f) result(II)
    real(C_DOUBLE), intent(in) :: f(ny0-2:nyN+2)
    real(C_DOUBLE) :: II, yp1, ym1, a1, a2, a3
    integer(C_INT) :: iy
    II=0.0d0
    DO iy=(ny0/2)*2+1,nyN,2
      yp1=y(iy+1)-y(iy); ym1=y(iy-1)-y(iy)
      a1=-1.0d0/3.0d0*ym1+1.0d0/6.0d0*yp1+1.0d0/6.0d0*yp1*yp1/ym1
      a3=+1.0d0/3.0d0*yp1-1.0d0/6.0d0*ym1-1.0d0/6.0d0*ym1*ym1/yp1
      a2=yp1-ym1-a1-a3
      II=II+a1*f(iy-1)+a2*f(iy)+a3*f(iy+1)
    END DO
  END FUNCTION yintegr

#define rD0(f,g,k) sum(dcmplx(der(iy)%d0(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d0(-2:2)*dreal(f(iy-2:iy+2,iz,ix,k))))
#define rD1(f,g,k) sum(dcmplx(der(iy)%d1(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d1(-2:2)*dreal(f(iy-2:iy+2,iz,ix,k))))
#define rD2(f,g,k) sum(dcmplx(der(iy)%d2(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d2(-2:2)*dreal(f(iy-2:iy+2,iz,ix,k))))
#define rD4(f,g,k) sum(dcmplx(der(iy)%d4(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d4(-2:2)*dreal(f(iy-2:iy+2,iz,ix,k))))
#define D0(f,g) sum(dcmplx(der(iy)%d0(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d0(-2:2)*dimag(f(iy-2:iy+2,iz,ix,g))))
#define D1(f,g) sum(dcmplx(der(iy)%d1(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d1(-2:2)*dimag(f(iy-2:iy+2,iz,ix,g))))
#define D2(f,g) sum(dcmplx(der(iy)%d2(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d2(-2:2)*dimag(f(iy-2:iy+2,iz,ix,g))))
#define D4(f,g) sum(dcmplx(der(iy)%d4(-2:2)*dreal(f(iy-2:iy+2,iz,ix,g)) ,der(iy)%d4(-2:2)*dimag(f(iy-2:iy+2,iz,ix,g))))
  !--------------------------------------------------------------!
  !---COMPLEX----- derivative in the y-direction ----------------!
#ifdef nonblockingY
  SUBROUTINE COMPLEXderiv(f0,f1,Rs,itag)
#else
  SUBROUTINE COMPLEXderiv(f0,f1)
#endif
  ! WARNING: TO GET THE PROPER DERIVATIVE, A CALL TO LeftLU5divStep2
  ! IS NEEDED AFTER CALLING THIS FUNCTION!
    complex(C_DOUBLE_COMPLEX), intent(in)  :: f0(ny0-2:nyN+2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: f1(ny0-2:nyN+2)
#ifdef nonblockingY
    TYPE(MPI_REQUEST), intent(out) :: Rs
    integer(C_INT), intent(in) :: itag
#endif
    IF (first) THEN 
      f1(0)=sum(d140(-2:2)*f0(-1:3))
      f1(-1)=sum(d14m1(-2:2)*f0(-1:3))
    END IF
    IF (last) THEN 
      f1(ny)=sum(d14n(-2:2)*f0(ny-3:ny+1))
      f1(ny+1)=sum(d14np1(-2:2)*f0(ny-3:ny+1))
    END IF
    DO CONCURRENT (iy=ny0:nyN)
      f1(iy)=sum(der(iy)%d1(-2:2)*f0(iy-2:iy+2))
    END DO
    IF (first) THEN
      f1(1)=f1(1)-(der(1)%d0(-1)*f1(0)+der(1)%d0(-2)*f1(-1))
      f1(2)=f1(2)-der(2)%d0(-2)*f1(0)
    END IF
    IF (last) THEN
      f1(ny-1)=f1(ny-1)-(der(ny-1)%d0(1)*f1(ny)+der(ny-1)%d0(2)*f1(ny+1))
      f1(ny-2)=f1(ny-2)-der(ny-2)%d0(2)*f1(ny)
    END IF
#ifdef nonblockingY
    CALL LeftLU5divStep1(f1,D0mat,f1,Rs,itag)
#else
    CALL LeftLU5divStep1(f1,D0mat,f1)
#endif
  END SUBROUTINE COMPLEXderiv

  !--------------------------------------------------------------!
  !---COMPLEX----- second derivative in the y-direction ---------!
#ifdef nonblockingY
  SUBROUTINE COMPLEXderiv2(f0,f1,Rs,itag)
#else
  SUBROUTINE COMPLEXderiv2(f0,f1)
#endif
  ! WARNING: TO GET THE PROPER DERIVATIVE, A CALL TO LeftLU5divStep2
  ! IS NEEDED AFTER CALLING THIS FUNCTION!
    complex(C_DOUBLE_COMPLEX), intent(in)  :: f0(ny0-2:nyN+2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: f1(ny0-2:nyN+2)
#ifdef nonblockingY
    TYPE(MPI_REQUEST), intent(out) :: Rs
    integer(C_INT), intent(in) :: itag
#endif
    IF (first) THEN 
      f1(0)=sum(d240(-2:2)*f0(-1:3))
      f1(-1)=sum(d24m1(-2:2)*f0(-1:3))
    END IF
    IF (last) THEN 
      f1(ny)=sum(d24n(-2:2)*f0(ny-3:ny+1))
      f1(ny+1)=sum(d24np1(-2:2)*f0(ny-3:ny+1))
    END IF
    DO CONCURRENT (iy=ny0:nyN)
      f1(iy)=sum(der(iy)%d2(-2:2)*f0(iy-2:iy+2))
    END DO
    IF (first) THEN
      f1(1)=f1(1)-(der(1)%d0(-1)*f1(0)+der(1)%d0(-2)*f1(-1))
      f1(2)=f1(2)-der(2)%d0(-2)*f1(0)
    END IF
    IF (last) THEN
      f1(ny-1)=f1(ny-1)-(der(ny-1)%d0(1)*f1(ny)+der(ny-1)%d0(2)*f1(ny+1))
      f1(ny-2)=f1(ny-2)-der(ny-2)%d0(2)*f1(ny)
    END IF
#ifdef nonblockingY
    CALL LeftLU5divStep1(f1,D0mat,f1,Rs,itag)
#else
    CALL LeftLU5divStep1(f1,D0mat,f1)
#endif
  END SUBROUTINE COMPLEXderiv2

  !-------------------------------------------------------------!
  !-----REAL----- derivative in the y-direction ----------------!
  SUBROUTINE REALderiv(f0,f1)
  ! WARNING: THIS ROUTINE ALREADY INCLUDES CALL TO LeftLU5divStep2
  ! NO NEED TO CALL IT!
  ! --------------------------------------------------------------
  ! TODO FIXME: this actually calls the complex routine!

    real(C_DOUBLE), intent(in)  :: f0(ny0-2:nyN+2)
    real(C_DOUBLE), intent(out) :: f1(ny0-2:nyN+2)
    complex(C_DOUBLE_COMPLEX) :: temp_in(ny0-2:nyN+2), temp_out(ny0-2:nyN+2)

    temp_in = dcmplx(f0)
    call COMPLEXderiv(temp_in, temp_out)
    call LeftLU5divStep2(D0mat, temp_out)
    f1 = dreal(temp_out)

  END SUBROUTINE REALderiv

  !--------------------------------------------------------------!
  !-----REAL----- second derivative in the y-direction ----------!
  SUBROUTINE REALderiv2(f0,f1)
  ! WARNING: THIS ROUTINE ALREADY INCLUDES CALL TO LeftLU5divStep2
  ! NO NEED TO CALL IT!
  ! --------------------------------------------------------------
  ! TODO FIXME: this actually calls the complex routine!

    real(C_DOUBLE), intent(in)  :: f0(ny0-2:nyN+2)
    real(C_DOUBLE), intent(out) :: f1(ny0-2:nyN+2)
    complex(C_DOUBLE_COMPLEX) :: temp_in(ny0-2:nyN+2), temp_out(ny0-2:nyN+2)

    temp_in = dcmplx(f0)
    call COMPLEXderiv2(temp_in, temp_out)
    call LeftLU5divStep2(D0mat, temp_out)
    f1 = dreal(temp_out)

  END SUBROUTINE REALderiv2

  !--------------------------------------------------------------!
  !----------------- apply the boundary conditions --------------!
  PURE SUBROUTINE applybc_0(EQ,bc0,bc0m1)
    real(C_DOUBLE), intent(inout) :: EQ(ny0:nyN+2,-2:2)
    real(C_DOUBLE), intent(in) :: bc0(-2:2),bc0m1(-2:2)
    EQ(1,-1:2)=EQ(1,-1:2)-EQ(1,-2)*bc0m1(-1:2)/bc0m1(-2)
    EQ(1, 0:2)=EQ(1, 0:2)-EQ(1,-1)*bc0(0:2)/bc0(-1)
    EQ(2,-1:1)=EQ(2,-1:1)-EQ(2,-2)*bc0(0:2)/bc0(-1)
  END SUBROUTINE applybc_0

  PURE SUBROUTINE applybc_n(EQ,bcn,bcnp1)
    real(C_DOUBLE), intent(inout) :: EQ(ny0:nyN+2,-2:2)
    real(C_DOUBLE), intent(in) :: bcn(-2:2),bcnp1(-2:2)
    EQ(ny-1,-2:1)=EQ(ny-1,-2:1)-EQ(ny-1,2)*bcnp1(-2:1)/bcnp1(2)
    EQ(ny-1,-2:0)=EQ(ny-1,-2:0)-EQ(ny-1,1)*bcn(-2:0)/bcn(1)
    EQ(ny-2,-1:1)=EQ(ny-2,-1:1)-EQ(ny-2,2)*bcn(-2:0)/bcn(1)
  END SUBROUTINE applybc_n

  !--------------------------------------------------------------!
  !------------------- solve the linear system  -----------------!
#define OS(iy,j) (ni*(der(iy)%d4(j)-2.0d0*k2(iz,ix)*der(iy)%d2(j)+k2(iz,ix)*k2(iz,ix)*der(iy)%d0(j)))
#define SQ(iy,j) (ni*(der(iy)%d2(j)-k2(iz,ix)*der(iy)%d0(j)))
#ifdef nonblockingY
#include "linsolve_nonblocking.inc"
#else
#include "linsolve_blocking.inc"
#endif

  !--------------------------------------------------------------!
  !------------------------ convolutions ------------------------!
#define timescheme(rhs,old,unkn,impl,expl) rhs=ODE(1)*(unkn)/deltat+(impl)+ODE(2)*(expl)-ODE(3)*(old); old=expl                   
  SUBROUTINE convolutions(iy,i,compute_cfl,in_timeloop,ODE)
     real(C_DOUBLE), intent(in) :: ODE(1:3)
     integer(C_INT), intent(in) :: iy,i
     logical, intent(in) :: compute_cfl, in_timeloop
     integer(C_INT) :: iV,ix,iz
#ifdef nonblockingXZ
     TYPE(MPI_REQUEST) :: Rs(1:6)
#endif
#ifdef ibm
     real(C_DOUBLE) :: newcoef,oldcoef
     real(C_DOUBLE) :: dU(1:3)
#endif
#ifdef convvel
     integer(C_INT) :: ix0
     real(C_DOUBLE) :: cu
     complex(C_DOUBLE_COMPLEX) :: dtu,ust
#endif
     VVdz(1:nz+1,1:nxB,1:3,i)=V(iy,0:nz,nx0:nxN,1:3);         VVdz(nz+2:nzd-nz,1:nxB,1:3,i)=0;
     VVdz(nzd+1-nz:nzd,1:nxB,1:3,i)=V(iy,-nz:-1,nx0:nxN,1:3); 
#ifdef ibm
     newcoef=ODE(2)/ODE(1); oldcoef=-ODE(3)/ODE(1)
#endif
     DO iV=1,3
       CALL IFT(VVdz(1:nzd,1:nxB,iV,i))  
#ifdef nonblockingXZ
       CALL MPI_IAlltoall(VVdz(:,:,iV,i), 1, Mdz, VVdx(:,:,iV,i), 1, Mdx, MPI_COMM_X, Rs(iV))
#endif
#ifdef convvel
       IF (compute_convvel) THEN
         IF (convvel_cnt>-1) THEN
           DO CONCURRENT (iz=1:nzd, ix=1:nxB) 
              ix0=ix-1+nx0; 
#define phase(n) ATAN2(dimag(n),dreal(n))
              IF (ix0>0) THEN 
                dtu = (VVdz(iz,ix,iV,i) - Voldz(iz,ix,iy,iV))/deltat
                ust = 0.5*(VVdz(iz,ix,iV,i) + Voldz(iz,ix,iy,iV))
                cu = dimag(dconjg(ust)*dtu)/(ix0*alfa0*ust*dconjg(ust))
                uconv(iz,ix,iy,iV)=  uconv(iz,ix,iy,iV) + cu
              END IF
           END DO
         END IF
         Voldz(1:nzd,1:nxB,iy,iV)=VVdz(1:nzd,1:nxB,iV,i)
       END IF
#endif
#ifndef nonblockingXZ
       CALL MPI_Alltoall(VVdz(:,:,iV,i), 1, Mdz, VVdx(:,:,iV,i), 1, Mdx, MPI_COMM_X)
       VVdx(nx+2:nxd+1,1:nzB,iV,i)=0;    CALL RFT(VVdx(1:nxd+1,1:nzB,iV,i),rVVdx(1:2*nxd+2,1:nzB,iV,i))
#endif
     END DO
#ifdef nonblockingXZ
     DO iV=1,3
       CALL MPI_wait(Rs(iV),MPI_STATUS_IGNORE); 
       VVdx(nx+2:nxd+1,1:nzB,iV,i)=0; 
       CALL RFT(VVdx(1:nxd+1,1:nzB,iV,i),rVVdx(1:2*nxd+2,1:nzB,iV,i));
     END DO
#endif
#ifdef convvel
     IF (iy==nyN+2 .AND. compute_convvel) THEN
       convvel_cnt=convvel_cnt+1
       compute_convvel=.FALSE.
     END IF
#endif
     VVdx(nx+2:nxd+1,1:nzB,4:6,i)=0
     IF (compute_cfl .and. iy>=1 .and. iy<=ny-1) THEN
           cfl=max(cfl,(maxval(abs(rVVdx(1:2*nxd,1:nzB,1,i))/dx     + &
                               abs(rVVdx(1:2*nxd,1:nzB,2,i))/dy(iy) + &
                               abs(rVVdx(1:2*nxd,1:nzB,3,i))/dz)))
     END IF
#ifdef ibm
     IF (in_timeloop) THEN
       rFdx=0; Fdz=0
       DO CONCURRENT (ix=1:2*nxd, iz=1:nzB, iV=1:3)
         ! Velocity difference
         dU(iV)=(Ut(ix,iz,iy,iV) - rVVdx(ix,iz,iV,i))*InBody(ix,iz,iy,1)
         ! Integral velocity difference
         dUint(ix,iz,iy,iV, 0)=dUint(ix,iz,iy,iV,0)+deltat*(dUint(ix,iz,iy,iV,-1)*oldcoef+dU(iV)*newcoef)
         dUint(ix,iz,iy,iV,-1)=dU(iV)
         !timescheme(dUint(ix,iz,iy,iV,0),dUint(ix,iz,iy,iV,-1),dUint(ix,iz,iy,iV,0),0,dU(iV));
         !dUint(ix,iz,iy,iV,0)=dUint(ix,iz,iy,iV,0)*deltat/ODE(1)
         ! Define body force
         rFdx(ix,iz,iV,1)=(dU(iV)*ibmp + dUint(ix,iz,iy,iV,0)*ibmi)
       END DO
       DO iV=1,3
         CALL HFT(rFdx(1:2*nxd+2,1:nzB,iV,1), Fdx(1:nxd+1,1:nzB,iV,1))
         CALL MPI_Alltoall(Fdx(:,:,iV,1), 1, Mdx, Fdz(:,:,iV,1), 1, Mdz, MPI_COMM_X)
         CALL FFT(Fdz(1:nzd,1:nxB,iV,1))
       END DO
       F(iy,0:nz,nx0:nxN,1:3)=Fdz(1:nz+1,1:nxB,1:3,1)*factor
       F(iy,-nz:-1,nx0:nxN,1:3)=Fdz(nzd+1-nz:nzd,1:nxB,1:3,1)*factor
     END IF
#endif
     rVVdx(1:2*nxd,1:nzB,4,i)  = rVVdx(1:2*nxd,1:nzB,1,i)  * rVVdx(1:2*nxd,1:nzB,2,i)*factor
     rVVdx(1:2*nxd,1:nzB,5,i)  = rVVdx(1:2*nxd,1:nzB,2,i)  * rVVdx(1:2*nxd,1:nzB,3,i)*factor
     rVVdx(1:2*nxd,1:nzB,6,i)  = rVVdx(1:2*nxd,1:nzB,1,i)  * rVVdx(1:2*nxd,1:nzB,3,i)*factor
     rVVdx(1:2*nxd,1:nzB,1:3,i)= rVVdx(1:2*nxd,1:nzB,1:3,i)* rVVdx(1:2*nxd,1:nzB,1:3,i)*factor
     DO iV=1,6
       CALL HFT(rVVdx(1:2*nxd+2,1:nzB,iV,i),VVdx(1:nxd+1,1:nzB,iV,i)); 
#ifndef nonblockingXZ
       CALL MPI_Alltoall(VVdx(:,:,iV,i), 1, Mdx, VVdz(:,:,iV,i), 1, Mdz, MPI_COMM_X)
       CALL FFT(VVdz(1:nzd,1:nxB,iV,i));
#else 
       CALL MPI_IAlltoall(VVdx(:,:,iV,i), 1, Mdx, VVdz(:,:,iV,i), 1, Mdz, MPI_COMM_X, Rs(iV))
#endif
     END DO
#ifdef nonblockingXZ
     DO iV=1,6
       CALL MPI_Wait(Rs(iV),MPI_STATUS_IGNORE)
       CALL FFT(VVdz(1:nzd,1:nxB,iV,i));
     END DO
#endif
  END SUBROUTINE convolutions


  !--------------------------------------------------------------!
  !-------------------------- buildRHS --------------------------!
  ! (u,v,w) = (1,2,3)
  ! (uu,vv,ww,uv,vw,uw) = (1,2,3,4,5,6)
#define DD(f,k) ( der(iy)%f(-2)*VVdz(izd(iz)+1,ix+1-nx0,k,im2)+der(iy)%f(-1)*VVdz(izd(iz)+1,ix+1-nx0,k,im1)+der(iy)%f(0)*VVdz(izd(iz)+1,ix+1-nx0,k,i0)+ \
                  der(iy)%f(1 )*VVdz(izd(iz)+1,ix+1-nx0,k,i1 )+der(iy)%f(2 )*VVdz(izd(iz)+1,ix+1-nx0,k,i2 ) )
  SUBROUTINE buildrhs(ODE,compute_cfl)
    logical, intent(in) :: compute_cfl
    real(C_DOUBLE), intent(in) :: ODE(1:3)
    integer(C_INT) :: iy,iz,ix,i,im2,im1,i0,i1,i2
    complex(C_DOUBLE_COMPLEX) :: rhsu,rhsv,rhsw,DD0_6,DD1_6,expl
#ifdef bodyforce
    IF (first) THEN 
        iy=1; F(-1:0,:,:,:)=0
        DO CONCURRENT (iz=-nz:nz, ix=nx0:nxN, i=1:3)
          F(-1,iz,ix,i)=-D4(F,i)/der(iy)%d4(-2) 
        END DO
    END IF 
    IF (last) THEN
        iy=ny-1; F(ny:ny+1,:,:,:)=0
        DO CONCURRENT (iz=-nz:nz, ix=nx0:nxN, i=1:3)
          F(ny+1,iz,ix,i)=-D4(F,i)/der(iy)%d4(2)
        END DO
    END IF
#endif
    DO iy=ny0-4,nyN+2 !-3,ny+1
      IF (iy<=nyN) THEN !(iy<=ny-1) THEN
      CALL convolutions(iy+2,imod(iy+2)+1,compute_cfl,.TRUE.,ODE)
      IF (iy>=ny0) THEN !(iy>=1) THEN
        im2=imod(iy-2)+1; im1=imod(iy-1)+1; i0=imod(iy)+1; i1=imod(iy+1)+1; i2=imod(iy+2)+1;
        DO iz=-nz,nz 
        DO ix=nx0,nxN
            DD0_6=DD(d0,6); DD1_6=DD(d1,6);
            rhsu=-ialfa(ix)*DD(d0,1)-DD(d1,4)-ibeta(iz)*DD0_6
            rhsv=-ialfa(ix)*DD(d0,4)-DD(d1,2)-ibeta(iz)*DD(d0,5)
            rhsw=-ialfa(ix)*DD0_6-DD(d1,5)-ibeta(iz)*DD(d0,3)
            expl=(ialfa(ix)*(ialfa(ix)*DD(d1,1)+DD(d2,4)+ibeta(iz)*DD1_6)+&
                  ibeta(iz)*(ialfa(ix)*DD1_6+DD(d2,5)+ibeta(iz)*DD(d1,3))-k2(iz,ix)*rhsv &
#ifdef bodyforce
                 - k2(iz,ix)*D0(F,2)-ialfa(ix)*D1(F,1)-ibeta(iz)*D1(F,3) &
#endif
                 )
            timescheme(newrhs(iy,iz,ix)%D2v, oldrhs(iy,iz,ix)%D2v, D2(V,2)-k2(iz,ix)*D0(V,2),sum(OS(iy,-2:2)*V(iy-2:iy+2,iz,ix,2)),expl); !(D2v)
            IF (ix==0 .AND. iz==0) THEN
              expl=(dcmplx(dreal(rhsu)+meanpx,dreal(rhsw)+meanpz) &
#ifdef bodyforce
                   +rD0(F,1,3) &
#endif
                   )
              timescheme(newrhs(iy,0,0)%eta,oldrhs(iy,0,0)%eta,rD0(V,1,3),ni*rD2(V,1,3),expl)!(Ubar, Wbar)
            ELSE
              expl=(ibeta(iz)*rhsu-ialfa(ix)*rhsw &
#ifdef bodyforce
                   +ibeta(iz)*D0(F,1)-ialfa(ix)*D0(F,3) &
#endif
                   )
              timescheme(newrhs(iy,iz,ix)%eta, oldrhs(iy,iz,ix)%eta,ibeta(iz)*D0(V,1)-ialfa(ix)*D0(V,3),sum(SQ(iy,-2:2)*[ibeta(iz)*V(iy-2:iy+2,iz,ix,1)-ialfa(ix)*V(iy-2:iy+2,iz,ix,3)]),expl) !(eta)
            END IF
        END DO
        END DO
      END IF
      END IF
      IF (iy-2>=ny0) THEN !(iy-2>=1) THEN
        DO CONCURRENT (ix=nx0:nxN, iz=-nz:nz) 
          V(iy-2,iz,ix,1) = newrhs(iy-2,iz,ix)%eta; V(iy-2,iz,ix,2) = newrhs(iy-2,iz,ix)%d2v; 
        END DO
      END IF      
    END DO
  END SUBROUTINE buildrhs

  !--------------------------------------------------------------!
  !-------------------- read_restart_file -----------------------! 
  SUBROUTINE read_restart_file(filename,R)
    complex(C_DOUBLE_COMPLEX), intent(INOUT) :: R(ny0-2:nyN+2,-nz:nz,nx0:nxN,1:3)
    character(len=40), intent(IN) :: filename
    integer(C_SIZE_T) :: ix,iy,iz,io
    INTEGER(MPI_OFFSET_KIND) :: disp = 3*C_INT + 7*C_DOUBLE
    TYPE(MPI_File) :: fh
    real(C_DOUBLE) :: rn(1:3)

    OPEN(UNIT=100,FILE=TRIM(filename),access="stream",status="old",action="read",iostat=io)
    IF (io==0) THEN
      if (has_terminal) print *, "Reading from file "//filename
      READ(100,POS=1) nx,ny,nz,alfa0,beta0,ni,a,ymin,ymax,time
      CLOSE(100)
      call MPI_file_open(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh)
      call MPI_file_set_view(fh, disp, MPI_DOUBLE_COMPLEX, vel_read_type, 'native', MPI_INFO_NULL)
      call MPI_file_read_all(fh, R, 1, vel_field_type, MPI_STATUS_IGNORE)
      call MPI_file_close(fh)
    ELSE
      CLOSE(100)
      R=0
      IF (has_terminal) WRITE(*,*) "Generating initial field..."
      DO iy=ny0-2,nyN+2; DO ix=nx0,nxN; DO iz=-nz,nz
          CALL RANDOM_NUMBER(rn)
          R(iy,iz,ix,1) = 0.01*EXP(dcmplx(0,rn(1)-0.5));  R(iy,iz,ix,2) = 0.01*EXP(dcmplx(0,rn(2)-0.5));  R(iy,iz,ix,3) = 0.01*EXP(dcmplx(0,rn(3)-0.5));
      END DO;        END DO;        END DO
      IF (has_average) THEN
        DO CONCURRENT (iy=ny0-2:nyN+2)
          R(iy,0,0,1)=3*0.5*y(iy)*(2-y(iy))  !+ 0.01*SIN(8*y(iy)*2*PI)/ni
          !V(iy,0,0,1)=y(iy)*(2-y(iy))*3.d0/2.d0 + 0.001*SIN(8*y(iy)*2*PI);
          !V(iy,0,0,1)=y(iy)-1
        END DO
      END IF
    END IF
  END SUBROUTINE read_restart_file

#ifdef ibm
  !--------------------------------------------------------------!
  !---------------------- read_body_file ------------------------! 
  SUBROUTINE read_body_file(filename,R)
    real(C_DOUBLE), intent(INOUT) :: R(1:,1:,ny0-2:,1:)
    character(len=40), intent(IN) :: filename
    integer(C_SIZE_T) :: iV,nV,ix,iy,iz,io,nxB_t,nx_t,nz_t,ny_t,iproc_t,br=8,bc=16,b1=1,b7=7,b3=3
    integer(C_SIZE_T) :: pos
    nV=SIZE(R,4)
    OPEN(UNIT=100,FILE=TRIM(filename),access="stream",status="old",action="read",iostat=io)
    nx_t=2*nxd; ny_t=ny+3; nz_t=nzd; nV=SIZE(R,4)
    IF (io==0) THEN
      IF (has_terminal) WRITE(*,*) "Reading "//TRIM(ADJUSTL(filename))//" ..."
       DO iV=1,nV
        DO iy=ny0-2,nyN+2
          DO iz=nz0,nzN
              pos=b1+br*(INT(iV-1,C_SIZE_T))*nx_t*nz_t*ny_t + &
                     br*(INT(iy+1,C_SIZE_T))*nx_t*nz_t + &
                     br*(INT(iz,  C_SIZE_T))*nx_t
              READ(100,POS=pos) R(:,iz-nz0+1,iy,iV)
          END DO
        END DO
       END DO
      CLOSE(100)
    ELSE
      IF (has_terminal) WRITE(*,*) TRIM(ADJUSTL(filename)) // " not loaded: file does not exist"
    END IF
  END SUBROUTINE read_body_file

  !--------------------------------------------------------------!
  !---------------------- save_body_file ------------------------! 
  SUBROUTINE save_body_file(filename,R)
    real(C_DOUBLE), intent(IN) :: R(1:,1:,ny0-2:,1:)
    character(len=40), intent(IN) :: filename
    integer(C_SIZE_T) :: iV,nV,ix,iy,iz,io,i,nxB_t,nx_t,nz_t,ny_t,iproc_t,br=8,bc=16,b1=1,b7=7,b3=3
    integer(C_SIZE_T) :: pos
    nV=SIZE(R,4)
    DO i=0,nproc-1
      IF (i==iproc) THEN
      OPEN(UNIT=100,FILE=TRIM(filename),access="stream",action="write")
      nx_t=2*nxd; ny_t=ny+3; nz_t=nzd
        IF (has_terminal) WRITE(*,*) "Writing "//TRIM(ADJUSTL(filename))//" ..."
        DO iV=1,nV
          DO iy=ny0-2,nyN+2
            DO iz=nz0,nzN
                pos=b1+br*(INT(iV-1,C_SIZE_T))*nx_t*nz_t*ny_t + &
                       br*(INT(iy+1,C_SIZE_T))*nx_t*nz_t + &
                       br*(INT(iz,  C_SIZE_T))*nx_t
                WRITE(100,POS=pos) R(:,iz-nz0+1,iy,iV)
            END DO
          END DO
        END DO
        CLOSE(100)
      END IF
      CALL MPI_Barrier(MPI_COMM_WORLD)
    END DO
  END SUBROUTINE save_body_file

!  !--------------------------------------------------------------!
!  !------------------ define where the body is ------------------!
!  PURE FUNCTION InBody(ix,iz,iy)
!     integer(C_INT), intent(in) :: ix,iz,iy
!     integer(C_INT) :: InBody
!     InBody = MERGE(1,0,( MOD(iz-1,nzd/8) < MAX(0,16-NINT(16*MERGE(ymax-y(iy),y(iy),y(iy)>1)/0.2)) ) ) 
!  END FUNCTION InBody
#endif

#ifdef convvel
  !--------------------------------------------------------------!
  !-------------------- save_convvel_file -----------------------! 
  SUBROUTINE save_convvel_file(filename,R)
    real(C_DOUBLE), intent(IN) :: R(1:,1:,ny0-2:,1:)
    character(len=40), intent(IN) :: filename
    integer(C_SIZE_T) :: iV,nV,ix,iy,iz,io,i,nxB_t,nx_t,nz_t,ny_t,iproc_t,br=8,bc=16,b1=1,b7=7,b3=3
    integer(C_SIZE_T) :: pos
    nV=SIZE(R,4)
    DO i=0,nproc-1
      IF (i==iproc) THEN
        OPEN(UNIT=100,FILE=TRIM(filename),access="stream",action="write")
        nx_t=nx+1; ny_t=ny+3; nz_t=nzd
        DO iV=1,nV
          DO iy=miny,maxy
            DO ix=nx0,nxN
              pos=b1+br*(INT(iV-1,C_SIZE_T))*nx_t*nz_t*ny_t + &
                     br*(INT(iy+1,C_SIZE_T))*nx_t*nz_t + &
                     br*(INT(ix,C_SIZE_T))*nz_t
              WRITE(100,POS=pos) R(:,ix-nx0+1,iy,iV)
            END DO
          END DO
         END DO
         CLOSE(100)
      END IF
      CALL MPI_Barrier(MPI_COMM_WORLD)
    END DO
  END SUBROUTINE save_convvel_file
#endif

  !--------------------------------------------------------------!
  !-------------------- save_restart_file -----------------------!
  SUBROUTINE save_restart_file(filename,R)
    complex(C_DOUBLE_COMPLEX), intent(in) :: R(ny0-2:nyN+2,-nz:nz,nx0:nxN,1:3)
    character(len=40), intent(in) :: filename
    ! mpi stuff
    TYPE(MPI_File) :: fh
    INTEGER(MPI_OFFSET_KIND) :: disp 
    TYPE(MPI_Status) :: status
    
    ! open file
    CALL MPI_File_open(MPI_COMM_WORLD, TRIM(filename), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh) 

    ! write header
    IF (has_terminal) THEN ! only one process does this
      CALL MPI_file_write(fh, [nx,ny,nz], 3, MPI_INTEGER, status)
      CALL MPI_file_write(fh, [alfa0,beta0,ni,a,ymin,ymax,time], 7, MPI_DOUBLE_PRECISION, status)
    END IF

    ! set view to subarray
    disp = 3*C_INT + 7*C_DOUBLE ! offset to skip header
    CALL MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, writeview_type, 'native', MPI_INFO_NULL)

    ! finally write field
    CALL MPI_File_write_all(fh, R, 1, owned2write_type, status)

    ! close file
    call MPI_File_close(fh)

  END SUBROUTINE save_restart_file


  !--------------------------------------------------------------!
  !------------------------- outstats ---------------------------!
  SUBROUTINE outstats()
   real(C_DOUBLE) :: runtime_global,dudy(1:2,1:2)   !cfl
   character(len=40) :: istring, filename
   TYPE(MPI_REQUEST) :: Rs,Rr
   TYPE(MPI_STATUS) :: S
#ifdef convvel
   compute_convvel=.TRUE.
#endif
   CALL MPI_Allreduce(cfl,runtime_global,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD); cfl=0;
   IF (cflmax>0)  deltat=cflmax/runtime_global;
   IF (ipx==0) THEN 
     IF (ipy==npy-1) THEN 
       dudy(1,2)=-sum(d14n(-2:2)*dreal(V(ny-3:ny+1,0,0,1))); dudy(2,2)=-sum(d14n(-2:2)*dreal(V(ny-3:ny+1,0,0,3)))
       CALL MPI_ISend(dudy(1:2,2),2,MPI_DOUBLE_PRECISION,0,TAG_DUDY,MPI_COMM_Y,Rs)
     END IF
     IF (ipy==0) THEN
       dudy(1,1)=sum(d140(-2:2)*dreal(V(-1:3,0,0,1)));       dudy(2,1)=sum(d140(-2:2)*dreal(V(-1:3,0,0,3)))
       CALL MPI_IRecv(dudy(1:2,2),2,MPI_DOUBLE_PRECISION,npy-1,TAG_DUDY,MPI_COMM_Y,Rr)
       CALL MPI_Wait(Rr,S);
     END IF
     IF (ipy==npy-1) CALL MPI_Wait(Rs,S)
   END IF
   IF (has_terminal) THEN
     WRITE(*,"(F10.4,3X,4(F11.6,3X),4(F9.4,3X),2(F9.6,3X))") &
           time,dudy(1,1),dudy(1,2),dudy(2,1),dudy(2,2),fr(1)+corrpx*fr(3),meanpx+corrpx,fr(2)+corrpz*fr(3),meanpz+corrpz,runtime_global*deltat,deltat
     WRITE(101,*) time,dudy(1,1),dudy(1,2),dudy(2,1),dudy(2,2),fr(1)+corrpx*fr(3),meanpx+corrpx,fr(2)+corrpz*fr(3),meanpz+corrpz,runtime_global*deltat,deltat
     FLUSH(101)
   END IF
   runtime_global=0
   IF (dt_save > 0) THEN ! save restart file
     IF ( ((FLOOR((time+0.5*deltat)/dt_save) > FLOOR((time-0.5*deltat)/dt_save)) .AND. (istep>1))) THEN
       IF (has_terminal) WRITE(*,*) "Writing Dati.cart.out at time ", time
       filename="Dati.cart.out";  CALL save_restart_file(filename,V)
     END IF
#ifdef ibm
     IF ( ((FLOOR((time+0.5*deltat)/dt_save) > FLOOR((time-0.5*deltat)/dt_save)) .AND. (istep>1))) THEN
       IF (has_terminal) WRITE(*,*) "Writing dUint.cart.out at time ", time
       filename="dUint.cart.out"; CALL save_body_file(filename,dUint(:,:,:,:,0))
     END IF
#endif
   END IF ! end of saving restart file
   ! Save i-th field
   IF ( time+deltat >= (ifield+1)*dt_field ) THEN ! fast evaluation
     IF ( prev_was_close .OR. ((FLOOR((time+0.5*deltat)/dt_field) > FLOOR((time-0.5*deltat)/dt_field)) .AND. (time>time0)) ) THEN ! accurate evaluation
       ifield=ifield+1; WRITE(istring,*) ifield
       IF (has_terminal) WRITE(*,*) "Writing Dati.cart."//TRIM(ADJUSTL(istring))//".out at time ", time
       filename="Dati.cart."//TRIM(ADJUSTL(istring))//".out"; CALL save_restart_file(filename,V)
#ifdef bodyforce
       IF (has_terminal) WRITE(*,*) "Writing Force.cart."//TRIM(ADJUSTL(istring))//".out at time ", time
       filename="Force.cart."//TRIM(ADJUSTL(istring))//".out"; CALL save_restart_file(filename,F)
#endif
#ifdef convvel
       IF (has_terminal) WRITE(*,*) "Writing Convvel.cart."//TRIM(ADJUSTL(istring))//".out at time ", time
       uconv=uconv/convvel_cnt
       filename="Convvel.cart."//TRIM(ADJUSTL(istring))//".out"; CALL save_convvel_file(filename,uconv)
       uconv=0; convvel_cnt=0
#endif
       prev_was_close = .FALSE.
     ELSE
       prev_was_close = .TRUE.
     END IF
   END IF
  END SUBROUTINE outstats


  !--------------------------------------------------------------!
  !------------------------- filedmap ---------------------------!
  ! access a field on disk without loading it in memory
  function fieldmap(unit, iy_zero, iz_zero, ix_zero, ic_zero) result(element)
        integer, intent(in) :: unit, ix_zero, iy_zero, iz_zero, ic_zero
        integer(C_SIZE_T) :: ix, iy, iz, ic
        complex(C_DOUBLE_COMPLEX) :: element
        integer(C_SIZE_T) :: position, el_idx ! position
        integer :: base = 3*C_INT + 7*C_DOUBLE ! base address
        
        ! calculate indices starting from 1
        ix = ix_zero + 1
        iz = iz_zero + nz + 1
        iy = iy_zero + 2
        ic = ic_zero

        ! calculate position
        el_idx = (ic-1_C_SIZE_T)*(nx+1_C_SIZE_T)*(2_C_SIZE_T*nz+1_C_SIZE_T)*(ny+3_C_SIZE_T) + (ix-1_C_SIZE_T)*(2_C_SIZE_T*nz+1_C_SIZE_T)*(ny+3_C_SIZE_T) + (iz-1_C_SIZE_T)*(ny+3_C_SIZE_T) + iy
        position = base + 1_C_SIZE_T + (el_idx - 1) * 2*C_DOUBLE_COMPLEX

        ! read element
        read(unit,pos=position) element

    end function


!--------------------------------------------------------------!
!---------------------- wanna_continue ------------------------!
! asks user whether to continue or not; possibly stops execution
    SUBROUTINE wanna_continue()
      character :: yn 

      IF (has_terminal) THEN
        WRITE(*,"(a,$)") "Do you want to continue? (y/n) "
        READ(*,*) yn
        IF (yn == 'n') STOP
      END IF
    END SUBROUTINE

END MODULE dnsdata
