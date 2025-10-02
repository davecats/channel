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

  ! V(iy,iz,ix,i) ,                                 i={1:u,2:v,3:w, 4 to 3+nPhi: scalar}
  ! oldrhs(iy,iz,ix,i),                         i={1:eta,2:d2v, 3 to 2+nPhi: scalar}
  ! bc0(iz,ix,i), bcn(iz,ix,i),                 i={1:u, 2:v, 3:w, 4:vy, 5:eta, 6 to 5+nPhi: scalar}
  ! der(iy,i,j),                                 i={0:d0, 1:d1, 2:d2, 3:d4}, j={-2:2}

  USE, intrinsic :: iso_c_binding
  USE rbmat
  USE mpi_transpose
  USE ffts

  IMPLICIT NONE

  !Simulation parameters
  real(C_DOUBLE) :: PI = 3.1415926535897932384626433832795028841971
  integer(C_INT) :: ny, nz, nxd, nPhi
  !$omp declare target(ny)
  real(C_DOUBLE) :: alfa0, beta0, ni, a, ymin, ymax, deltat, cflmax, time, time0 = 0, dt_field, dt_save, t_max, gamma
  real(C_DOUBLE) :: u0, uN, t0, tN
  real(C_DOUBLE) :: meanpx, meanpz, meanflowx, meanflowz, meantx, meantb
  integer(C_INT), allocatable :: izd(:)
  complex(C_DOUBLE_COMPLEX), allocatable :: ialfa(:), ibeta(:)
  real(C_DOUBLE), allocatable :: k2(:, :)
  logical :: time_from_restart
  !Grid
  integer(C_INT), private :: iy
  real(C_DOUBLE), allocatable :: y(:), dy(:)
  real(C_DOUBLE) :: dx, dz, factor
  !Derivatives
  real(C_DOUBLE), allocatable :: der(:, :, :)
  real(C_DOUBLE), dimension(-2:2) :: d040, d140, d240, d14m1, d24m1, d04n, d14n, d24n, d14np1, d24np1
  !$omp declare target(d14np1, d14n, d14m1, d140)
  real(C_DOUBLE), allocatable :: D0mat(:, :), linsolve_mat(:, :, :, :), eta00mat(:, :)
#if !(defined(HAVE_CUDA) || defined(HAVE_HIP))
  !Fourier-transformable arrays (allocated in ffts.f90)
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:, :, :, :) :: VVdx, VVdz
  real(C_DOUBLE), pointer, dimension(:, :, :, :) :: rVVdx
#endif
  !Solution
  complex(C_DOUBLE_COMPLEX), allocatable :: memrhs(:, :, :, :), oldrhs(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: V(:, :, :, :)
#ifdef bodyforce
  complex(C_DOUBLE_COMPLEX), allocatable :: F(:, :, :, :)
#endif
  !Boundary conditions
  real(C_DOUBLE), dimension(-2:2) :: v0bc, v0m1bc, vnbc, vnp1bc, eta0bc, eta0m1bc, etanbc, etanp1bc,phi0bc,phi0m1bc,phinbc,phinp1bc
  complex(C_DOUBLE_COMPLEX), allocatable :: bc0(:, :, :), bcn(:, :, :)
  !Mean pressure correction
  real(C_DOUBLE), private :: corrpx = 0.d0, corrpz = 0.d0
  complex(C_DOUBLE_COMPLEX), allocatable :: ucor(:), tcor(:, :)
  real(C_DOUBLE), dimension(:), allocatable :: corrtx, pra
  !ODE Library
  real(C_DOUBLE), dimension(3, 3) :: RK_rai = reshape( &
                                     (/120.0d0/32.0d0, 2.0d0, 0.0d0, &
                                       120.0d0/8.0d0, 50.0d0/8.0d0, 34.0d0/8.0d0, &
                                       120.0d0/20.0d0, 90.0d0/20.0d0, 50.0d0/20.0d0/), &
                                     shape=(/3, 3/))
  !Outstats
  real(C_DOUBLE) :: cfl = 0.0d0
  integer(C_SIZE_T) :: istep, nstep, ifield
  real(C_DOUBLE), dimension(:), allocatable :: fr
  !Restart file
  character(len=40) :: fname

CONTAINS

  !--------------------------------------------------------------!
  !---------------------- Read input files ----------------------!
  SUBROUTINE read_dnsin(filename)
    IMPLICIT NONE
    logical :: i
    integer :: iPhi
    CHARACTER(len=*), INTENT(IN) :: filename
    OPEN (15, file=filename)
    READ (15, *) nx, ny, nz; READ (15, *) alfa0, beta0; nxd = 3*(nx + 1)/2; nzd = 3*nz
    !$omp target update to(ny)
#ifdef useFFTfit
    i = fftFIT(nxd); DO WHILE (.NOT. i); nxd = nxd + 1; i = fftFIT(nxd); END DO
    i = fftFIT(nzd); DO WHILE (.NOT. i); nzd = nzd + 1; i = fftFIT(nzd); END DO
#endif
    READ (15, *) ni; 
    READ (15, *) a, ymin, ymax; ni = 1/ni
    READ (15, *) meanpx, meanpz
    READ (15, *) meanflowx, meanflowz
    READ (15, *) meantx, meantb
    READ (15, *) u0, uN, t0, tN
    READ (15, *) nPhi
    allocate (pra(nPhi))
    READ (15, *) pra(1:nPhi)
    DO iPhi = 1, nPhi
      pra(iPhi) = 1/pra(iPhi)
    END DO
    !$omp target enter data map(to: pra)
    READ (15, *) deltat, cflmax, time
    READ (15, *) dt_field, dt_save, t_max, time_from_restart
    READ (15, *) nstep
    CLOSE (15)
    dx = PI/(alfa0*nxd); dz = 2.0d0*PI/(beta0*nzd); factor = 1.0d0/(2.0d0*nxd*nzd)
  END SUBROUTINE read_dnsin

  !--------------------------------------------------------------!
  !---------------- Allocate memory for solution ----------------!
  SUBROUTINE init_memory(solveNS)
    IMPLICIT NONE
    INTEGER(C_INT) :: ix, iz
    logical, intent(IN) :: solveNS
    ALLOCATE (V(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3 + nPhi)); V = 0
    !$omp target enter data map(to: V)
#ifdef bodyforce
    ALLOCATE (F(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3)); F = 0
    !$omp target enter data map(to: F)
#endif
    IF (solveNS) then
      ALLOCATE(memrhs(1:ny-1,-nz:nz,nx0:nxN,1:2+nPhi),oldrhs(1:ny-1,-nz:nz,nx0:nxN,1:2+nPhi),bc0(-nz:nz,nx0:nxN,1:5+nPhi),bcn(-nz:nz,nx0:nxN,1:5+nPhi)); 
      memrhs = 0.0; oldrhs = 0.0; bc0 = 0.0; bcn = 0.0
      !$omp target enter data map(to: memrhs, oldrhs, bc0, bcn)
    END IF
#define newrhs(iy,iz,ix,i) memrhs(iy,iz,ix,i)
#define imod(iy) MOD(iy+1000,5)
    ALLOCATE (der(1:ny - 1, 0:3, -2:2), d0mat(ny0:nyN + 2, -2:2), linsolve_mat(ny0:nyN + 2, -2:2, -nz:nz, nx0:nxN))
    ALLOCATE (eta00mat(ny0:nyN + 2, -2:2))
    der = 0.0; d0mat = 0.0; linsolve_mat = 0.0; eta00mat = 0.0
    allocate (corrtx(1:nPhi)); corrtx = 0.0d0
    allocate (ucor(ny0 - 2:nyN + 2), tcor(ny0 - 2:nyN + 2, 1:nPhi))
    ALLOCATE (y(-1:ny + 1), dy(1:ny - 1))
    y = 0.0; dy = 0.0
    ALLOCATE (izd(-nz:nz), ialfa(nx0:nxN), ibeta(-nz:nz), k2(-nz:nz, nx0:nxN))
#ifdef halfchannel
    FORALL (iy=-1:ny + 1) y(iy) = ymin + (ymax - ymin)*(tanh(a*(real(iy)/real(ny) - 1))/tanh(a) + 1)
#else
    FORALL (iy=-1:ny + 1) y(iy) = ymin + 0.5d0*(ymax - ymin)*(tanh(a*(2.0d0*real(iy)/real(ny) - 1))/tanh(a) + 1)
#endif
    FORALL (iy=1:ny - 1) dy(iy) = 0.5d0*(y(iy + 1) - y(iy - 1))
    izd = (/(merge(iz, nzd + iz, iz >= 0), iz=-nz, nz)/); ialfa = (/(dcmplx(0.0d0, ix*alfa0), ix=nx0, nxN)/); 
    ibeta = (/(dcmplx(0.0d0, iz*beta0), iz=-nz, nz)/); 
    FORALL (iz=-nz:nz, ix=nx0:nxN) k2(iz, ix) = (alfa0*ix)**2.0d0 + (beta0*iz)**2.0d0
    !$omp target enter data map(to: izd, ialfa, ibeta, k2, y, iy, rk_rai, linsolve_mat, ucor, tcor)
    OPEN (UNIT=195, FILE='Runtimedata.phi', ACTION='write')
    IF (solveNS .AND. has_terminal) OPEN (UNIT=121, FILE='Runtimedata', ACTION='write')

    allocate (fr(3 + 2*nPhi)); fr = 0.0
  END SUBROUTINE init_memory

  !--------------------------------------------------------------!
  !--------------- Deallocate memory for solution ---------------!
  SUBROUTINE free_memory(solveNS)
    IMPLICIT NONE
    LOGICAL, intent(IN) :: solveNS
    !$omp target exit data map(delete: d240, d24m1, d04n, d24n, d24np1, D0mat)
    !$omp target exit data map(delete: V)
    !$omp target exit data map(delete: memrhs, oldrhs, bc0, bcn, linsolve_mat)
    !$omp target exit data map(delete: izd, ialfa, ibeta, k2, ucor, tcor)
#ifdef bodyforce
    !$omp target exit data map(delete: F)
#endif
    DEALLOCATE (V, der, d0mat, linsolve_mat, y, dy)
    IF (solveNS) THEN
      DEALLOCATE (memrhs, oldrhs, bc0, bcn)
      IF (has_terminal) CLOSE (UNIT=121)
    END IF
  END SUBROUTINE free_memory

  !--------------------------------------------------------------!
  !--------------- Set-up the compact derivatives ---------------!
  SUBROUTINE setup_derivatives()
    IMPLICIT NONE
    real(C_DOUBLE)    :: M(0:4, 0:4), t(0:4)
    integer(C_INT)    :: iy, i, j
    DO iy = 1, ny - 1
      FORALL (i=0:4, j=0:4) M(i, j) = (y(iy - 2 + j) - y(iy))**(4.0d0 - i); CALL LUdecomp(M)
      t = 0; t(0) = 24
      der(iy, 3, -2:2) = M.bs.t
      FORALL (i=0:4, j=0:4) M(i, j) = (5.0d0 - i)*(6.0d0 - i)*(7.0d0 - i)*(8.0d0 - i)*(y(iy - 2 + j) - y(iy))**(4.0d0 - i); CALL LUdecomp(M)
      FORALL (i=0:4) t(i) = sum(der(iy, 3, -2:2)*(y(iy - 2:iy + 2) - y(iy))**(8.0d0 - i))
      der(iy, 0, -2:2) = M.bs.t
      FORALL (i=0:4, j=0:4) M(i, j) = (y(iy - 2 + j) - y(iy))**(4.0d0 - i); CALL LUdecomp(M)
      t = 0; FORALL (i=0:2) t(i) = sum(der(iy, 0, -2:2)*(4.0d0 - i)*(3.0d0 - i)*(y(iy - 2:iy + 2) - y(iy))**(2.0d0 - i))
      der(iy, 2, -2:2) = M.bs.t
      t = 0; FORALL (i=0:3) t(i) = sum(der(iy, 0, -2:2)*(4.0d0 - i)*(y(iy - 2:iy + 2) - y(iy))**(3.0d0 - i))
      der(iy, 1, -2:2) = M.bs.t
    END DO
    FORALL (i=0:4, j=0:4) M(i, j) = (y(-1 + j) - y(0))**(4.0d0 - i); CALL LUdecomp(M)
    t = 0; t(3) = 1.0; d140(-2:2) = M.bs.t
    t = 0; t(2) = 2.0; d240(-2:2) = M.bs.t
    FORALL (i=0:4, j=0:4) M(i, j) = (y(-1 + j) - y(-1))**(4.0d0 - i); CALL LUdecomp(M)
    t = 0; t(3) = 1.0; d14m1(-2:2) = M.bs.t
    t = 0; t(2) = 2.0; d24m1(-2:2) = M.bs.t
    d040 = 0; d040(-1) = 1
    FORALL (i=0:4, j=0:4) M(i, j) = (y(ny - 3 + j) - y(ny))**(4.0d0 - i); CALL LUdecomp(M)
    t = 0; t(3) = 1; d14n(-2:2) = M.bs.t
    t = 0; t(2) = 2; d24n(-2:2) = M.bs.t
    FORALL (i=0:4, j=0:4) M(i, j) = (y(ny - 3 + j) - y(ny + 1))**(4.0d0 - i); CALL LUdecomp(M)
    t = 0; t(3) = 1; d14np1(-2:2) = M.bs.t
    t = 0; t(2) = 2; d24np1(-2:2) = M.bs.t
    d04n = 0; d04n(1) = 1; 
    FORALL (iy=1:ny - 1) D0mat(iy, -2:2) = der(iy, 0, -2:2); 
    CALL LU5decomp(D0mat)
    !$omp target update to(d14np1, d14n, d14m1, d140)
    !$omp target enter data map(to: d240, d24m1, d04n, d14n, d24n, d24np1, D0mat, der)
  END SUBROUTINE setup_derivatives

  !--------------------------------------------------------------!
  !--------------- Set-up the boundary conditions ---------------!
  SUBROUTINE setup_boundary_conditions()
    IMPLICIT NONE
    integer :: ix, iz
    ! Bottom wall
    v0bc = d040; v0m1bc = d140; eta0bc = d040
    eta0m1bc = der(1, 3, :)
    phi0bc = d040; phi0m1bc = der(1, 3, :)        ! Dirichlet
#ifdef phiNeumann
    phi0bc = d140; phi0m1bc = der(1, 3, :)        ! Neumann
#endif
    v0bc(-1:2) = v0bc(-1:2) - v0bc(-2)*v0m1bc(-1:2)/v0m1bc(-2)
    eta0bc(-1:2) = eta0bc(-1:2) - eta0bc(-2)*eta0m1bc(-1:2)/eta0m1bc(-2)
    phi0bc(-1:2) = phi0bc(-1:2) - phi0bc(-2)*phi0m1bc(-1:2)/phi0m1bc(-2)
    ! Top wall
#ifdef halfchannel
    vnbc = d04n; vnp1bc = d24n; etanbc = d14n
#else
    vnbc = d04n; vnp1bc = d14n; etanbc = d04n
#endif
    etanp1bc = der(ny - 1, 3, :)
    phinbc = d04n; phinp1bc = der(ny - 1, 3, :) ! Dirichlet
#ifdef phiNeumann
    phinbc = d14n; phinp1bc = d04n
#endif
    vnbc(-2:1) = vnbc(-2:1) - vnbc(2)*vnp1bc(-2:1)/vnp1bc(2)
    etanbc(-2:1) = etanbc(-2:1) - etanbc(2)*etanp1bc(-2:1)/etanp1bc(2)
    phinbc(-2:1) = phinbc(-2:1) - phinbc(2)*phinp1bc(-2:1)/phinp1bc(2)
    !$omp target enter data map(to: v0bc, v0m1bc, vnbc, vnp1bc, eta0bc, eta0m1bc, etanbc, etanp1bc, phinbc, phi0bc, phi0m1bc, phinp1bc)

    !precompute bc0 and bcn
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(bc0, bcn, ialfa, ibeta, v0bc, v0m1bc, vnbc, vnp1bc, nx0, nxN, nz) &
    !$omp private(ix, iz)
    DO ix = nx0, nxN
      DO iz = -nz, nz
        IF (ix == 0 .AND. iz == 0) THEN
          bc0(iz, ix, 1) = 0
          bc0(iz, ix, 4) = 0
          bc0(iz, ix, 5) = dcmplx(dreal(bc0(iz, ix, 1)) - dimag(bc0(iz, ix, 3)), dimag(bc0(iz, ix, 1)) + dreal(bc0(iz, ix, 3)))
        ELSE
          bc0(iz, ix, 4) = -ialfa(ix)*bc0(iz, ix, 1) - ibeta(iz)*bc0(iz, ix, 3); bc0(iz, ix, 5) = ibeta(iz)*bc0(iz, ix, 1) - ialfa(ix)*bc0(iz, ix, 3)
        END IF
        bc0(iz, ix, 2) = bc0(iz, ix, 2) - v0bc(-2)*bc0(iz, ix, 4)/v0m1bc(-2)

        IF (ix == 0 .AND. iz == 0) THEN
          bcn(iz, ix, 2) = 0
          bcn(iz, ix, 4) = 0
          bcn(iz, ix, 5) = dcmplx(dreal(bcn(iz, ix, 1)) - dimag(bcn(iz, ix, 3)), dimag(bcn(iz, ix, 1)) + dreal(bcn(iz, ix, 3)))
        ELSE
          bcn(iz, ix, 4) = -ialfa(ix)*bcn(iz, ix, 1) - ibeta(iz)*bcn(iz, ix, 3); bcn(iz, ix, 5) = ibeta(iz)*bcn(iz, ix, 1) - ialfa(ix)*bcn(iz, ix, 3)
        END IF
        bcn(iz, ix, 2) = bcn(iz, ix, 2) - vnbc(2)*bcn(iz, ix, 4)/vnp1bc(2)
      END DO
    END DO
  END SUBROUTINE setup_boundary_conditions

  !--------------------------------------------------------------!
  !---------------- integral in the y-direction -----------------!
#ifdef HAVE_CUDA
  !$omp declare target(yintegr)
#endif
  PURE FUNCTION yintegr(f, y) result(II)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: f(ny0 - 2:nyN + 2)
    real(C_DOUBLE), intent(in) :: y(-1:ny + 1)
    real(C_DOUBLE) :: II, yp1, ym1, a1, a2, a3
    integer(C_INT) :: iy
    II = 0.0d0
    DO iy = (ny0/2)*2 + 1, nyN, 2
      yp1 = y(iy + 1) - y(iy); ym1 = y(iy - 1) - y(iy)
      a1 = -1.0d0/3.0d0*ym1 + 1.0d0/6.0d0*yp1 + 1.0d0/6.0d0*yp1*yp1/ym1
      a3 = +1.0d0/3.0d0*yp1 - 1.0d0/6.0d0*ym1 - 1.0d0/6.0d0*ym1*ym1/yp1
      a2 = yp1 - ym1 - a1 - a3
      II = II + a1*dreal(f(iy - 1)) + a2*dreal(f(iy)) + a3*dreal(f(iy + 1))
    END DO
  END FUNCTION yintegr

#define RD_MASTER(f,g,k,ord) ( \
  der(iy, ord, -2)*dcmplx(dreal(f(iy - 2, iz, ix, g)), dreal(f(iy - 2, iz, ix, k))) + \
  der(iy, ord, -1)*dcmplx(dreal(f(iy - 1, iz, ix, g)), dreal(f(iy - 1, iz, ix, k))) + \
  der(iy, ord, 0)*dcmplx(dreal(f(iy, iz, ix, g)), dreal(f(iy, iz, ix, k))) + \
  der(iy, ord, 1)*dcmplx(dreal(f(iy + 1, iz, ix, g)), dreal(f(iy + 1, iz, ix, k))) + \
  der(iy, ord, 2)*dcmplx(dreal(f(iy + 2, iz, ix, g)), dreal(f(iy + 2, iz, ix, k))))

#define rD0(f,g,k) RD_MASTER(f,g,k,0)
#define rD1(f,g,k) RD_MASTER(f,g,k,1)
#define rD2(f,g,k) RD_MASTER(f,g,k,2)
#define rD4(f,g,k) RD_MASTER(f,g,k,3)

#define D_MASTER(f,g,ord) ( \
  der(iy, ord, -2)*f(iy - 2, iz, ix, g) + \
  der(iy, ord, -1)*f(iy - 1, iz, ix, g) + \
  der(iy, ord, 0)*f(iy, iz, ix, g) + \
  der(iy, ord, 1)*f(iy + 1, iz, ix, g) + \
  der(iy, ord, 2)*f(iy + 2, iz, ix, g))

#define D0(f,g) D_MASTER(f,g,0)
#define D1(f,g) D_MASTER(f,g,1)
#define D2(f,g) D_MASTER(f,g,2)
#define D4(f,g) D_MASTER(f,g,3)
  !--------------------------------------------------------------!
  !---COMPLEX----- derivative in the y-direction ----------------!
#ifdef HAVE_CUDA
  !$omp declare target(COMPLEXderiv)
#endif
  SUBROUTINE COMPLEXderiv(f0, f1, der, D0mat)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in)  :: f0(-1:ny + 1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: f1(-1:ny + 1)
    real(C_DOUBLE), DIMENSION(:, :), intent(in) :: der(ny0:nyN, 0:3, -2:2)
    real(C_DOUBLE), DIMENSION(:, :), intent(in) :: D0mat(ny0:nyN + 2, -2:2)
    integer(C_INT) :: iy
    f1(0) = sum(d140(-2:2)*f0(-1:3))
    f1(-1) = sum(d14m1(-2:2)*f0(-1:3))
    f1(ny) = sum(d14n(-2:2)*f0(ny - 3:ny + 1))
    f1(ny + 1) = sum(d14np1(-2:2)*f0(ny - 3:ny + 1))
    DO iy = ny0, nyN
      f1(iy) = sum(der(iy, 1, -2:2)*f0(iy - 2:iy + 2))
    END DO
    f1(1) = f1(1) - (der(1, 0, -1)*f1(0) + der(1, 0, -2)*f1(-1))
    f1(2) = f1(2) - der(2, 0, -2)*f1(0)
    f1(ny - 1) = f1(ny - 1) - (der(ny - 1, 0, 1)*f1(ny) + der(ny - 1, 0, 2)*f1(ny + 1))
    f1(ny - 2) = f1(ny - 2) - der(ny - 2, 0, 2)*f1(ny)
    CALL LeftLU5div(f1, D0mat, f1)
  END SUBROUTINE COMPLEXderiv

  !--------------------------------------------------------------!
  !----------------- apply the boundary conditions --------------!
#ifdef HAVE_CUDA
  !$omp declare target(applybc_0)
#endif
  PURE SUBROUTINE applybc_0(EQ, bc0, bc0m1)
    real(C_DOUBLE), intent(inout) :: EQ(ny0:nyN + 2, -2:2)
    real(C_DOUBLE), intent(in) :: bc0(-2:2), bc0m1(-2:2)
    EQ(1, -1:2) = EQ(1, -1:2) - EQ(1, -2)*bc0m1(-1:2)/bc0m1(-2)
    EQ(1, 0:2) = EQ(1, 0:2) - EQ(1, -1)*bc0(0:2)/bc0(-1)
    EQ(2, -1:1) = EQ(2, -1:1) - EQ(2, -2)*bc0(0:2)/bc0(-1)
  END SUBROUTINE applybc_0

#ifdef HAVE_CUDA
  !$omp declare target(applybc_n)
#endif
  PURE SUBROUTINE applybc_n(EQ, bcn, bcnp1)
    real(C_DOUBLE), intent(inout) :: EQ(ny0:nyN + 2, -2:2)
    real(C_DOUBLE), intent(in) :: bcn(-2:2), bcnp1(-2:2)
    EQ(ny - 1, -2:1) = EQ(ny - 1, -2:1) - EQ(ny - 1, 2)*bcnp1(-2:1)/bcnp1(2)
    EQ(ny - 1, -2:0) = EQ(ny - 1, -2:0) - EQ(ny - 1, 1)*bcn(-2:0)/bcn(1)
    EQ(ny - 2, -1:1) = EQ(ny - 2, -1:1) - EQ(ny - 2, 2)*bcn(-2:0)/bcn(1)
  END SUBROUTINE applybc_n

! Orr-Sommerfeld and Squire opearators
#define OS(iy,j) (ni*(der(iy,3,j)-2.0d0*k2(iz,ix)*der(iy,2,j)+k2(iz,ix)*k2(iz,ix)*der(iy,0,j)))
#define SQ(iy,j) (ni*(der(iy,2,j)-k2(iz,ix)*der(iy,0,j)))

#ifdef HAVE_CUDA
  !$omp declare target(LU5decomp)
#endif
  !---- in-place LU Decomposition of a banded matrix ---!
  !-----------------------------------------------------!
  SUBROUTINE LU5decomp(A)
    IMPLICIT NONE
    real(C_DOUBLE), intent(inout) :: A(0:, -2:)
    integer(C_INT) :: HI1, HI2
    real(C_DOUBLE) :: piv
    INTEGER :: i, k, j
    HI1 = SIZE(A, 1) - 1; HI2 = SIZE(A, 2) - 3; 
    A(HI1 - 2, 1:2) = 0; A(HI1 - 3, 2) = 0
    DO i = HI1 - HI2, 0, -1
      DO k = HI2, 1, -1
        piv = A(i, k)
        DO j = -1, -2, -1
          A(i, j + k) = A(i, j + k) - piv*A(i + k, j)
        END DO
      END DO
      piv = 1.0d0/A(i, 0); A(i, 0) = piv; A(i, -2:-1) = A(i, -2:-1)*piv
    END DO
    A(0, -2:-1) = 0; A(1, -2) = 0
  END SUBROUTINE LU5decomp

  !- Left LU division of a banded matrix -!
  !---------------------------------------!
#ifdef HAVE_CUDA
  !$omp declare target(LeftLU5div)
#endif
  SUBROUTINE LeftLU5div(x, A, b)
    complex(C_DOUBLE_COMPLEX), intent(out) :: x(-2:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: b(-2:)
    real(C_DOUBLE), intent(in)  :: A(0:, -2:)
    integer(C_INT) :: HI1, HI2, i
    HI1 = SIZE(A, 1) - 1
    HI2 = SIZE(A, 2) - 3

    ! initialise x with rhs

    DO i = LBOUND(x, 1), UBOUND(x, 1)
      x(i) = b(i)
    END DO

    ! backward substitution
    DO i = HI1 - HI2, 0, -1
      x(i) = x(i) - (A(i, 1)*x(i + 1) + A(i, 2)*x(i + 2))
      x(i) = x(i)*A(i, 0)
    END DO

    ! forward substitution
    DO i = 0, HI1
      x(i) = x(i) - (A(i, -2)*x(i - 2) + A(i, -1)*x(i - 1))
    END DO

  END SUBROUTINE LeftLU5div

  !--------------------------------------------------------------!
  !------------------- solve the linear system  -----------------!
  SUBROUTINE linsolve(lambda)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: lambda
    integer(C_INT) :: ix, iz, i, j, iPhi
    complex(C_DOUBLE_COMPLEX) :: temp

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(V, linsolve_mat, bc0, bcn, der, lambda, ni, ny, k2, ny0, nyN, nx0, nxN, nz, v0bc, v0m1bc, vnbc, vnp1bc) &
    !$omp private(ix, iz, iy)
    DO ix = nx0, nxN
      DO iz = -nz, nz
        DO iy = ny0, nyN
          linsolve_mat(iy, -2:2, iz, ix) = lambda*(der(iy, 2, -2:2) - k2(iz, ix)*der(iy, 0, -2:2)) - OS(iy, -2:2)
        END DO

        CALL applybc_0(linsolve_mat(:, :, iz, ix), v0bc, v0m1bc)
        V(1, iz, ix, 2) = V(1, iz, ix, 2) - linsolve_mat(1, -2, iz, ix)*bc0(iz, ix, 4)/v0m1bc(-2) - linsolve_mat(1, -1, iz, ix)*bc0(iz, ix, 2)/v0bc(-1)
        V(2, iz, ix, 2) = V(2, iz, ix, 2) - linsolve_mat(2, -2, iz, ix)*bc0(iz, ix, 2)/v0bc(-1)

        CALL applybc_n(linsolve_mat(:, :, iz, ix), vnbc, vnp1bc)
        V(ny - 1, iz, ix, 2) = V(ny - 1, iz, ix, 2) - linsolve_mat(ny - 1, 2, iz, ix)*bcn(iz, ix, 4)/vnp1bc(2) - linsolve_mat(ny - 1, 1, iz, ix)*bcn(iz, ix, 2)/vnbc(1)
        V(ny - 2, iz, ix, 2) = V(ny - 2, iz, ix, 2) - linsolve_mat(ny - 2, 2, iz, ix)*bcn(iz, ix, 2)/vnbc(1)

        CALL LU5decomp(linsolve_mat(:, :, iz, ix)); 
        CALL LeftLU5div(V(:, iz, ix, 2), linsolve_mat(:, :, iz, ix), V(:, iz, ix, 2))

        V(0, iz, ix, 2) = (bc0(iz, ix, 2) - sum(V(1:3, iz, ix, 2)*v0bc(0:2)))/v0bc(-1)
        V(-1, iz, ix, 2) = (bc0(iz, ix, 4) - sum(V(0:3, iz, ix, 2)*v0m1bc(-1:2)))/v0m1bc(-2)
        V(ny, iz, ix, 2) = (bcn(iz, ix, 2) - sum(V(ny - 3:ny - 1, iz, ix, 2)*vnbc(-2:0)))/vnbc(1)
        V(ny + 1, iz, ix, 2) = (bcn(iz, ix, 4) - sum(V(ny - 3:ny, iz, ix, 2)*vnp1bc(-2:1)))/vnp1bc(2)
      END DO
    END DO

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(V, linsolve_mat, bc0, bcn, der, k2, ny, ni, lambda, ny0, nyN, nx0, nxN, nz, eta0bc, eta0m1bc, etanbc, etanp1bc) &
    !$omp private(ix, iz, iy)
    DO ix = nx0, nxN
      DO iz = -nz, nz
        DO iy = ny0, nyN
          linsolve_mat(iy, -2:2, iz, ix) = lambda*der(iy, 0, -2:2) - SQ(iy, -2:2)
        END DO

        CALL applybc_0(linsolve_mat(:, :, iz, ix), eta0bc, eta0m1bc)
        V(1, iz, ix, 1) = V(1, iz, ix, 1) - linsolve_mat(1, -1, iz, ix)*bc0(iz, ix, 5)/eta0bc(-1)
        V(2, iz, ix, 1) = V(2, iz, ix, 1) - linsolve_mat(2, -2, iz, ix)*bc0(iz, ix, 5)/eta0bc(-1)

        CALL applybc_n(linsolve_mat(:, :, iz, ix), etanbc, etanp1bc)
        V(ny - 1, iz, ix, 1) = V(ny - 1, iz, ix, 1) - linsolve_mat(ny - 1, 1, iz, ix)*bcn(iz, ix, 5)/etanbc(1)
        V(ny - 2, iz, ix, 1) = V(ny - 2, iz, ix, 1) - linsolve_mat(ny - 2, 2, iz, ix)*bcn(iz, ix, 5)/etanbc(1)

        CALL LU5decomp(linsolve_mat(:, :, iz, ix))
        CALL LeftLU5div(V(:, iz, ix, 1), linsolve_mat(:, :, iz, ix), V(:, iz, ix, 1))

        V(0, iz, ix, 1) = (bc0(iz, ix, 5) - sum(V(1:3, iz, ix, 1)*eta0bc(0:2)))/eta0bc(-1)
        V(-1, iz, ix, 1) = -sum(V(0:3, iz, ix, 1)*eta0m1bc(-1:2))/eta0m1bc(-2)
        V(ny, iz, ix, 1) = (bcn(iz, ix, 5) - sum(V(ny - 3:ny - 1, iz, ix, 1)*etanbc(-2:0)))/etanbc(1)
        V(ny + 1, iz, ix, 1) = -sum(V(ny - 3:ny, iz, ix, 1)*etanp1bc(-2:1))/etanp1bc(2)
      END DO
    END DO

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ny, ny0, nyN, nz, nx0, nxN) &
    !$omp shared(V, bc0, bcn, v0bc, v0m1bc, eta0bc, eta0m1bc, vnbc, vnp1bc, etanbc, etanp1bc, linsolve_mat) &
    !$omp shared(D0mat, der, k2, ialfa, ibeta, y, fr, ucor, meanflowx, meanflowz, corrpx, corrpz) &
    !$omp private(ix, iz, j, temp)
    DO ix = nx0, nxN
      DO iz = -nz, nz
        ! Correct flow rate
        IF (ix == 0 .AND. iz == 0) THEN
          V(:, 0, 0, 3) = dcmplx(dimag(V(:, 0, 0, 1)), 0.d0); 
          V(:, 0, 0, 1) = dcmplx(dreal(V(:, 0, 0, 1)), 0.d0); 
          ucor(ny0 - 2:ny0 - 1) = 0; ucor(ny0:nyN) = 1; ucor(nyN + 1:nyN + 2) = 0
          !linsolve_mat contains etamat
          CALL LeftLU5div(ucor, linsolve_mat(:, :, iz, ix), ucor)
          ucor(0) = -sum(ucor(1:3)*eta0bc(0:2))/eta0bc(-1)
          ucor(-1) = -sum(ucor(0:3)*eta0m1bc(-1:2))/eta0m1bc(-2)
          ucor(ny) = -sum(ucor(ny - 3:ny - 1)*etanbc(-2:0))/etanbc(1)
          ucor(ny + 1) = -sum(ucor(ny - 3:ny)*etanp1bc(-2:1))/etanp1bc(2)

          fr(1) = yintegr(V(:, 0, 0, 1), y)
          fr(2) = yintegr(V(:, 0, 0, 3), y)
          fr(3) = yintegr(ucor, y)
          IF (abs(meanflowx) > 1.0d-7) THEN
            corrpx = (meanflowx - fr(1))/fr(3)
            V(:, 0, 0, 1) = dcmplx(dreal(V(:, 0, 0, 1)) + corrpx*dreal(ucor), dimag(V(:, 0, 0, 1)))
          END IF
          IF (abs(meanflowz) > 1.0d-7) THEN
            corrpz = (meanflowz - fr(2))/fr(3)
            V(:, 0, 0, 3) = dcmplx(dreal(V(:, 0, 0, 3)) + corrpz*dreal(ucor), dimag(V(:, 0, 0, 3)))
          END IF
        ELSE
          CALL COMPLEXderiv(V(:, iz, ix, 2), V(:, iz, ix, 3), der, D0mat)
          do j = ny0 - 2, nyN + 2
            temp = (ialfa(ix)*V(j, iz, ix, 3) - ibeta(iz)*V(j, iz, ix, 1))/k2(iz, ix)
            V(j, iz, ix, 3) = (ibeta(iz)*V(j, iz, ix, 3) + ialfa(ix)*V(j, iz, ix, 1))/k2(iz, ix)
            V(j, iz, ix, 1) = temp
          end do
        END IF
      END DO
    END DO
  END SUBROUTINE linsolve

  SUBROUTINE linsolve_scalar(lambda, iPhi)
    IMPLICIT NONE
    integer(C_INT), intent(in) :: iPhi
    real(C_DOUBLE), intent(in) :: lambda
    integer(C_INT) :: ix, iz, i, j
    complex(C_DOUBLE_COMPLEX) :: temp
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(V, linsolve_mat, bc0, bcn, der, lambda, pra, ni, k2, ny, ny0, nyN, nx0, nxN, nz, phi0bc, phi0m1bc, phinbc, phinp1bc, iPhi) &
    !$omp private(ix, iz, iy)
    DO ix = nx0, nxN
      DO iz = -nz, nz
        DO iy = ny0, nyN
          linsolve_mat(iy, -2:2, iz, ix) = lambda*der(iy, 0, -2:2) - pra(iPhi)*(SQ(iy, -2:2))
        END DO

        CALL applybc_0(linsolve_mat(:, :, iz, ix), phi0bc, phi0m1bc)
        V(1, iz, ix, 3 + iPhi) = V(1, iz, ix, 3 + iPhi) - linsolve_mat(1, -1, iz, ix)*bc0(iz, ix, 5 + iPhi)/phi0bc(-1)
        V(2, iz, ix, 3 + iPhi) = V(2, iz, ix, 3 + iPhi) - linsolve_mat(2, -2, iz, ix)*bc0(iz, ix, 5 + iPhi)/phi0bc(-1)

        CALL applybc_n(linsolve_mat(:, :, iz, ix), phinbc, phinp1bc)
        V(ny - 1, iz, ix, 3 + iPhi) = V(ny - 1, iz, ix, 3 + iPhi) - linsolve_mat(ny - 1, 1, iz, ix)*bcn(iz, ix, 5 + iPhi)/phinbc(1)
        V(ny - 2, iz, ix, 3 + iPhi) = V(ny - 2, iz, ix, 3 + iPhi) - linsolve_mat(ny - 2, 2, iz, ix)*bcn(iz, ix, 5 + iPhi)/phinbc(1)

        CALL LU5decomp(linsolve_mat(:, :, iz, ix))
        CALL LeftLU5div(V(:, iz, ix, 3 + iPhi), linsolve_mat(:, :, iz, ix), V(:, iz, ix, 3 + iPhi))

        V(0, iz, ix, 3 + iPhi) = (bc0(iz, ix, 5 + iPhi) - sum(V(1:3, iz, ix, 3 + iPhi)*phi0bc(0:2)))/phi0bc(-1)
        V(-1, iz, ix, 3 + iPhi) = -sum(V(0:3, iz, ix, 3 + iPhi)*phi0m1bc(-1:2))/phi0m1bc(-2)
        V(ny, iz, ix, 3 + iPhi) = (bcn(iz, ix, 5 + iPhi) - sum(V(ny - 3:ny - 1, iz, ix, 3 + iPhi)*phinbc(-2:0)))/phinbc(1)
        V(ny + 1, iz, ix, 3 + iPhi) = -sum(V(ny - 3:ny, iz, ix, 3 + iPhi)*phinp1bc(-2:1))/phinp1bc(2)
      END DO
    END do

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ny, ny0, nyN, nz, nx0, nxN, nPhi) &
    !$omp shared(V, bc0, bcn, phi0bc, phi0m1bc, phinbc, phinp1bc,linsolve_mat, tcor, y, fr, corrtx, meantb, iPhi) &
    !$omp private(ix, iz, temp)
    DO ix = nx0, nxN
      DO iz = -nz, nz
        ! Correct flow rate
        IF (ix == 0 .AND. iz == 0) THEN
          tcor(ny0 - 2:ny0 - 1, :) = 0
          tcor(ny0:nyN, :) = 1
          tcor(nyN + 1:nyN + 2, :) = 0

          CALL LeftLU5div(tcor(:, iPhi), linsolve_mat(:, :, iz, ix), tcor(:, iPhi))
          tcor(0, iPhi) = -sum(tcor(1:3, iPhi)*phi0bc(0:2))/phi0bc(-1)
          tcor(-1, iPhi) = -sum(tcor(0:3, iPhi)*phi0m1bc(-1:2))/phi0m1bc(-2)
          tcor(ny, iPhi) = -sum(tcor(ny - 3:ny - 1, iPhi)*phinbc(-2:0))/phinbc(1)
          tcor(ny + 1, iPhi) = -sum(tcor(ny - 3:ny, iPhi)*phinp1bc(-2:1))/phinp1bc(2)
          fr(3 + iPhi) = yintegr(V(:, 0, 0, 3 + iPhi), y); 
          fr(3 + nPhi + iPhi) = yintegr(tcor(:, iPhi), y); 
          IF (abs(meantb) > 1.0d-7) THEN
            corrtx(iPhi) = (meantb - fr(3 + iPhi))/fr(3 + nPhi + iPhi)
         V(:, 0, 0, 3 + iPhi) = dcmplx(dreal(V(:, 0, 0, 3 + iPhi)) + corrtx(iPhi)*dreal(tcor(:, iPhi)), dimag(V(:, 0, 0, 3 + iPhi)))
          END IF
        END IF
      END DO
    END DO
  END SUBROUTINE linsolve_scalar

  !--------------------------------------------------------------!
  !------------------------ convolutions ------------------------!

  SUBROUTINE assemble_vvdz(m)
    IMPLICIT NONE
    integer(C_INT) :: iy
    integer(C_INT) :: i, j, k
    integer(C_INT), intent(in) :: m
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, VVdz) shared(ny, nxB, nzd, nx0, nxN, nz, m) private(i,j,k)
    DO i = 1, ny + 3
      DO j = 1, nxB
        DO k = 1, nzd
          IF (k <= nz + 1) THEN
            VVdz(k, j, i, m) = V(i - 2, k - 1, j + nx0 - 1, m)
          ELSEIF (k >= nz + 2 .AND. k <= nzd - nz) THEN
            VVdz(k, j, i, m) = 0.0
          ELSE
            VVdz(k, j, i, m) = V(i - 2, k - nzd - 1, j + nx0 - 1, m)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE assemble_vvdz

  SUBROUTINE zero_vvdx_hft(m)
    IMPLICIT NONE
    integer(C_INT) :: iy
    integer(C_INT) :: i, j, k
    integer(C_INT), intent(in) :: m
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdx) shared(nx, nxd, nzB, ny, m) private(i,j,k)
    DO i = 1, ny + 3
      DO j = 1, nzB
        DO k = nx + 2, nxd + 1
          VVdx(k, j, i, m) = 0.0
        END DO
      END DO
    END DO
  END SUBROUTINE zero_vvdx_hft

  subroutine compute_cfl()
    implicit none
    integer(C_INT) :: i, j, k
    real(C_DOUBLE) :: tmp
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp private(i,j,k,tmp) reduction(max:cfl) &
    !$omp shared(rVVdx, dx, dy, dz, ny, nxd, nzB)
    do j = 1, 2*nxd
      do k = 1, nzB
        do i = 3, ny + 1
          tmp = abs(rVVdx(j, k, i, 1))/dx + abs(rVVdx(j, k, i, 2))/dy(i - 2) + abs(rVVdx(j, k, i, 3))/dz
          cfl = max(cfl, tmp)
        end do
      end do
    end do
  END SUBROUTINE compute_cfl

  SUBROUTINE build_products(m)
    IMPLICIT NONE
    integer(C_INT), intent(in) :: m
    integer(C_INT) :: iPhi, component, first, second
    integer(C_INT) :: i, j, k
    real(C_DOUBLE) :: a, b

    if (m > 6) then ! scalar cases 7, ...
      iPhi = (m - 4)/3
      component = mod(m - 4, 3) + 1
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(rVVdx) shared(nxd, nzB, ny, factor, iPhi, iPhi, component)
      DO i = 1, ny + 3
        DO j = 1, nzB
          DO k = 1, 2*nxd
            rVVdx(k, j, i, 6 + 3*(iPhi - 1) + component) = rVVdx(k, j, i, component)*rVVdx(k, j, i, 3 + iPhi)*factor
          END DO
        END DO
      END DO
    else if (m > 3) then !cases 4, 5, 6
      first = mod(m - 1, 3) + 1
      second = mod(m, 3) + 1
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(rVVdx) shared(nxd, nzB, ny, factor, first, second, m) private(a, b)
      DO i = 1, ny + 3
        DO j = 1, nzB
          DO k = 1, 2*nxd
            a = rVVdx(k, j, i, first)
            b = rVVdx(k, j, i, second)
            rVVdx(k, j, i, m) = a*b*factor
          END DO
        END DO
      END DO
    else ! cases 1, 2, 3
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(rVVdx) shared(nxd, nzB, ny, factor, m) private(a)
      DO i = 1, ny + 3
        DO j = 1, nzB
          DO k = 1, 2*nxd
            a = rVVdx(k, j, i, m)
            rVVdx(k, j, i, m) = a*a*factor
          END DO
        END DO
      END DO
    end if
  END SUBROUTINE build_products

  SUBROUTINE transform_to_physical()
    IMPLICIT NONE
    integer(C_INT) ::  m, mm1
    type(MPI_Request), dimension(:) :: requests(3 + nPhi)
    type(MPI_Status)  :: status

    DO m = 1, 3 + nPhi + 1
      mm1 = m - 1
      if (m <= 3 + nPhi) then
        !assemble, pack and post nonblocking alltoall
        CALL assemble_vvdz(m)
        CALL IFT(VVdz(:, :, :, m), ny)
        call pack_zTOx(VVdz(:, :, :, m), sendbuf(:, m), ny)
        call alltoall(sendbuf(:, m), recvbuf(:, m), requests(m))
      end if
      if (m > 1) then
        !wait on the previous step, unpack, and fft
        CALL MPI_WAIT(requests(mm1), status, ierr)
        CALL unpack_zTOx(recvbuf(:, mm1), VVdx(:, :, :, mm1), ny)
        CALL zero_vvdx_hft(mm1)
        CALL RFT(VVdx(:, :, :, mm1), rVVdx(:, :, :, mm1), ny)
      end if
    END DO
  END SUBROUTINE transform_to_physical

  SUBROUTINE transform_back_and_build_rhs(ODE)
    IMPLICIT NONE
    integer(C_INT) ::  m, mp1
    real(C_DOUBLE), intent(in) :: ODE(1:3)
    type(MPI_Request), dimension(:) :: requests(6 + 3*nPhi)
    type(MPI_Status)  :: status

    ! Reverse pass to build rVVdx
    DO m = 6 + 3*nPhi, 0, -1
      mp1 = m + 1
      if (m >= 1) then
        ! Build, HFT, pack and post nonblocking alltoall
        CALL build_products(m)
        CALL HFT(rVVdx(:, :, :, m), VVdx(:, :, :, m), ny)
        call pack_xTOz(VVdx(:, :, :, m), sendbuf(:, m), ny)
        call alltoall(sendbuf(:, m), recvbuf(:, m), requests(m))
      end if

      if (m < 6 + 3*nPhi) then
        ! Wait on the previous step, unpack, and fft
        CALL MPI_WAIT(requests(mp1), status, ierr)
        call unpack_xTOz(recvbuf(:, mp1), VVdz(:, :, :, mp1), ny)
        CALL FFT(VVdz(:, :, :, mp1), ny)
        CALL buildrhs(ODE, mp1)
      end if
    END DO
  END SUBROUTINE transform_back_and_build_rhs

  !--------------------------------------------------------------!
  !-------------------------- buildRHS --------------------------!
  ! (u,v,w) = (1,2,3)
  ! (uu,vv,ww,uv,vw,uw) = (1,2,3,4,5,6)
#define timescheme(rhs,old,unkn,impl,expl) rhs=ODE(1)*(unkn)/deltat+(impl)+ODE(2)*(expl)-ODE(3)*(old); old=expl
#define timescheme_accum(rhs,old,expl) rhs = rhs + ODE(2)*(expl); old = old + expl
#define DD(f,k) ( der(iy,f,-2)*VVdz(izd(iz)+1,ix+1-nx0,iy,k)+ \
  der(iy, f, -1)*VVdz(izd(iz) + 1, ix + 1 - nx0, iy + 1, k) + \
  der(iy, f, 0)*VVdz(izd(iz) + 1, ix + 1 - nx0, iy + 2, k) + \
  der(iy, f, 1)*VVdz(izd(iz) + 1, ix + 1 - nx0, iy + 3, k) + \
  der(iy, f, 2)*VVdz(izd(iz) + 1, ix + 1 - nx0, iy + 4, k))
#define ACCUM_D2V(new, expl) \
  timescheme_accum(new(iy, iz, ix, 2), oldrhs(iy, iz, ix, 2), (expl))

#define ACCUM_ETA_UW(new, rhsu, rhsw) \
  ACCUM_ETA_UW2(new, ibeta(iz)*(rhsu) - ialfa(ix)*(rhsw), dcmplx(dreal((rhsu)), dreal((rhsw))))

#define ACCUM_ETA_UW2(new, expl_general, expl_mean) \
  if (ix == 0 .AND. iz == 0) then; \
    timescheme_accum(new(iy, iz, ix, 1), oldrhs(iy, iz, ix, 1), (expl_mean)); \
  else; \
    timescheme_accum(new(iy, iz, ix, 1), oldrhs(iy, iz, ix, 1), (expl_general)); \
  end if
  SUBROUTINE buildrhs_prepare(ODE)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: ODE(1:3)
    integer(C_INT) :: iy, iz, ix, i, im2, im1, i0, i1, i2, k, iPhi
    complex(C_DOUBLE_COMPLEX) :: rhsu, rhsw, rhst, expl, tmp, unkn
#ifdef bodyforce
    iy = 1; F(-1:0, :, :, :) = 0
    DO iz = -nz, nz
      DO ix = nx0, nxN
        DO i = 1, 3
          F(-1, iz, ix, i) = -D4(F, i)/der(iy, 3, -2)
        END DO
      END DO
    END DO
    iy = ny - 1; F(ny:ny + 1, :, :, :) = 0
    DO iz = -nz, nz
      DO ix = nx0, nxN
        DO i = 1, 3
          F(ny + 1, iz, ix, i) = -D4(F, i)/der(iy, 3, 2)
        END DO
      END DO
    END DO
#endif

    ! contribution known a-priori
    !$omp target teams distribute parallel do collapse(3) default(none)  &
    !$omp private(iz, ix, iy, tmp, k, unkn) &
    !$omp shared(nz, nx0, nxN, ny) shared(ialfa, ibeta) shared(k2, der) shared(memrhs, oldrhs) shared(meanpx, meanpz, ni, deltat, ode) &
    !$omp shared(vvdz, v, pra)
    DO iz = -nz, nz
      DO ix = nx0, nxN
        DO iy = 1, ny - 1
          unkn = D2(V, 2) - k2(iz, ix)*D0(V, 2)
          tmp = 0.0
          DO k = -2, 2
            tmp = tmp + OS(iy, k)*V(iy + k, iz, ix, 2)
          END DO

          !initialize D2v
          newrhs(iy, iz, ix, 2) = ODE(1)*unkn/deltat + tmp - ODE(3)*oldrhs(iy, iz, ix, 2)
          oldrhs(iy, iz, ix, 2) = 0.0
        END DO
      END DO
    END DO

    !$omp target teams distribute parallel do collapse(3) default(none)  &
    !$omp private(iz, ix, iy, tmp, k, unkn) &
    !$omp shared(nz, nx0, nxN, ny) shared(ialfa, ibeta) shared(k2, der) shared(memrhs, oldrhs) shared(meanpx, meanpz, ni, deltat, ode) &
    !$omp shared(vvdz, v, pra)
    DO iz = -nz, nz
      DO ix = nx0, nxN
        DO iy = 1, ny - 1
          IF (ix == 0 .AND. iz == 0) THEN
            unkn = rD0(V, 1, 3)
            tmp = ni*rD2(V, 1, 3)
          ELSE
            tmp = 0.0
            DO k = -2, 2
              tmp = tmp + SQ(iy, k)*(ibeta(iz)*V(iy + k, iz, ix, 1) - ialfa(ix)*V(iy + k, iz, ix, 3))
            END DO
            unkn = ibeta(iz)*D0(V, 1) - ialfa(ix)*D0(V, 3)
          END IF
          !initialize eta
          newrhs(iy, iz, ix, 1) = ODE(1)*unkn/deltat + tmp - ODE(3)*oldrhs(iy, iz, ix, 1)
          oldrhs(iy, iz, ix, 1) = 0.0
        END DO
      END DO
    END DO

    !$omp target teams distribute parallel do collapse(3) default(none)  &
    !$omp private(rhsu, rhsw, expl) private(iz, ix, iy, tmp, k, unkn) &
    !$omp shared(nz, nx0, nxN, ny) shared(ialfa, ibeta) shared(k2, der) shared(memrhs, oldrhs) shared(meanpx, meanpz, ni, deltat, ode) &
    !$omp shared(vvdz, v, pra)
    DO iz = -nz, nz
      DO ix = nx0, nxN
        DO iy = 1, ny - 1

          if (ix == 0 .AND. iz == 0) then
            timescheme_accum(newrhs(iy, iz, ix, 1), oldrhs(iy, iz, ix, 1), dcmplx(meanpx, meanpz))
          end if
#ifdef bodyforce
          ACCUM_D2V(newrhs, -k2(iz, ix)*D0(F, 2) - ialfa(ix)*D1(F, 1) - ibeta(iz)*D1(F, 3))
          ACCUM_ETA_UW2(newrhs, ibeta(iz)*D0(F, 1) - ialfa(ix)*D0(F, 3), rD0(F, 1, 3))
#endif
        end do
      end do
    end do

    !initialize phi
    DO iPhi = 1, nPhi
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(tmp, iz, ix, iy, k, unkn) &
      !$omp shared(nz, nx0, nxN, ny, pra, V, memrhs, ode, deltat, oldrhs, iphi, der, k2, ni)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = 1, ny - 1
            tmp = 0.0
            DO k = -2, 2
              tmp = tmp + SQ(iy, k)*V(iy + k, iz, ix, 3 + iPhi)
            END DO
            unkn = D0(V, 3 + iPhi)
            newrhs(iy, iz, ix, 2 + iPhi) = ODE(1)*unkn/deltat + tmp*pra(iPhi) - ODE(3)*oldrhs(iy, iz, ix, 2 + iPhi)
            oldrhs(iy, iz, ix, 2 + iPhi) = 0.0
          END DO
        END DO
      END DO
    END DO

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(memrhs, V) shared(nz, nx0, nxN, ny, nPhi) private(iy, ix, iz)
    DO iz = -nz, nz
    DO ix = nx0, nxN
    DO iy = 1, ny - 1
      V(iy, iz, ix, 1) = newrhs(iy, iz, ix, 1); 
      V(iy, iz, ix, 2) = newrhs(iy, iz, ix, 2); 
      DO iPhi = 1, nPhi
        V(iy, iz, ix, 3 + iPhi) = newrhs(iy, iz, ix, 2 + iPhi); 
      END DO
    END DO
    END DO
    END DO

  END SUBROUTINE buildrhs_prepare

  SUBROUTINE buildrhs(ODE, component)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: ODE(1:3)
    integer(C_INT), intent(in) :: component
    integer(C_INT) :: iy, iz, ix, i, iPhi
    complex(C_DOUBLE_COMPLEX) :: rhsu, rhsw, rhst, expl, tmp, unkn

    SELECT CASE (component)
    CASE (1)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = 1, ny - 1

            rhsu = -ialfa(ix)*DD(0, 1); rhsw = 0.0
            expl = ialfa(ix)*ialfa(ix)*DD(1, 1)

            ACCUM_D2V(V, expl)
            ACCUM_ETA_UW(V, rhsu, rhsw)
          END DO
        END DO
      END DO
    CASE (2)
      !contribution from VVdz(:,:,:,2)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = 1, ny - 1

            expl = DD(1, 2)*k2(iz, ix)
            ACCUM_D2V(V, expl)
          END DO
        END DO
      END DO
    CASE (3)
      !contribution from VVdz(:,:,:,3)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = 1, ny - 1
            rhsw = -ibeta(iz)*DD(0, 3); rhsu = 0.0
            expl = ibeta(iz)*ibeta(iz)*DD(1, 3)

            ACCUM_D2V(V, expl)
            ACCUM_ETA_UW(V, rhsu, rhsw)
          END DO
        END DO
      END DO
    CASE (4)
      !contribution from VVdz(:,:,:,4)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = 1, ny - 1

            rhsu = -DD(1, 4); rhsw = 0.0
            expl = ialfa(ix)*DD(2, 4) + ialfa(ix)*DD(0, 4)*k2(iz, ix)

            ACCUM_D2V(V, expl)
            ACCUM_ETA_UW(V, rhsu, rhsw)
          END DO
        END DO
      END DO
    CASE (5)
      !contribution from VVdz(:,:,:,5)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = 1, ny - 1

            ! contribution from VVdz(:,:,:,5)
            rhsw = -DD(1, 5); rhsu = 0.0
            expl = ibeta(iz)*DD(2, 5) + ibeta(iz)*DD(0, 5)*k2(iz, ix)

            ACCUM_D2V(V, expl)
            ACCUM_ETA_UW(V, rhsu, rhsw)
          END DO
        END DO
      END DO
    CASE (6)
      !contribution from VVdz(:,:,:,6)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = 1, ny - 1

            rhsu = -ibeta(iz)*DD(0, 6)
            rhsw = -ialfa(ix)*DD(0, 6)
            expl = 2*ialfa(ix)*ibeta(iz)*DD(1, 6)

            ACCUM_D2V(V, expl)
            ACCUM_ETA_UW(V, rhsu, rhsw)

          END DO
        END DO
      END DO
    CASE DEFAULT
      iPhi = (component - 4)/3
      SELECT CASE (MODULO((component - 4), 3) + 1)
      CASE (1)
        !$omp target teams distribute parallel do collapse(3) default(none)  &
        !$omp private(rhst, iz, ix, iy) &
        !$omp shared(nz, nx0, nxN, ny, ialfa, V, ode, der, oldrhs, iphi, vvdz, izd)
        DO iz = -nz, nz
          DO ix = nx0, nxN
            DO iy = 1, ny - 1
              rhst = -ialfa(ix)*DD(0, 6 + 3*(iPhi - 1) + 1)
              timescheme_accum(V(iy, iz, ix, 3 + iPhi), oldrhs(iy, iz, ix, 2 + iPhi), rhst)
            END DO
          END DO
        END DO
      CASE (2)
        !$omp target teams distribute parallel do collapse(3) default(none)  &
        !$omp private(rhst, iz, ix, iy) &
        !$omp shared(nz, nx0, nxN, ny, V, ode, der, oldrhs, iphi, vvdz, izd)
        DO iz = -nz, nz
          DO ix = nx0, nxN
            DO iy = 1, ny - 1

              rhst = -DD(1, 6 + 3*(iPhi - 1) + 2)
              timescheme_accum(V(iy, iz, ix, 3 + iPhi), oldrhs(iy, iz, ix, 2 + iPhi), rhst)
            END DO
          END DO
        END DO
      CASE (3)
        !$omp target teams distribute parallel do collapse(3) default(none)  &
        !$omp private(rhst, iz, ix, iy) &
        !$omp shared(nz, nx0, nxN, ny, ibeta,V, ode, der, oldrhs, iphi, vvdz, izd)
        DO iz = -nz, nz
          DO ix = nx0, nxN
            DO iy = 1, ny - 1
              rhst = -ibeta(iz)*DD(0, 6 + 3*(iPhi - 1) + 3)
              timescheme_accum(V(iy, iz, ix, 3 + iPhi), oldrhs(iy, iz, ix, 2 + iPhi), rhst)
            END DO
          END DO
        END DO
      END SELECT
    END SELECT

  END SUBROUTINE buildrhs

  !--------------------------------------------------------------!
  !-------------------- read_restart_file -----------------------!
  SUBROUTINE read_restart_file(filename, R)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(INOUT) :: R(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3 + nPhi)
    character(len=*), intent(IN) :: filename
    integer(C_SIZE_T) :: ix, iy, iz, io, iPhi
    integer(C_INT) :: r_nx, r_ny, r_nz
    real(C_DOUBLE) :: r_alfa0, r_beta0, r_ni, r_a, r_ymin, r_ymax
    real(C_DOUBLE) :: rn(1:3)
#ifdef HAVE_MPI
    INTEGER(MPI_OFFSET_KIND) :: disp = 3*C_INT + 7*C_DOUBLE
    TYPE(MPI_File) :: fh

    OPEN (UNIT=120, FILE=TRIM(filename), access="stream", status="old", action="read", iostat=io)
    IF (io == 0) THEN
      if (has_terminal) print *, "Reading from file "//filename
      READ (120, POS=1) r_nx, r_ny, r_nz, r_alfa0, r_beta0, r_ni, r_a, r_ymin, r_ymax, time
      call MPI_file_open(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh)
      call MPI_file_set_view(fh, disp, MPI_DOUBLE_COMPLEX, vel_read_type, 'native', MPI_INFO_NULL)
      call MPI_file_read_all(fh, R, 1, vel_field_type, MPI_STATUS_IGNORE)
      call MPI_file_close(fh)
      IF (r_nx /= nx .OR. r_ny /= ny .OR. r_nz /= nz .OR. r_alfa0 /= alfa0 .OR. r_beta0 /= beta0 .OR. r_ni /= ni .OR. r_a /= a .OR. r_ymin /= ymin .OR. r_ymax /= ymax) THEN
        IF (has_terminal) PRINT *, "ERROR: mismatch in metadata between restart file and dns.in. Stopping."
        IF (has_terminal) PRINT *, "From .out file:"
        IF (has_terminal) PRINT *, r_nx, r_ny, r_nz, r_alfa0, r_beta0, r_ni, r_a, r_ymin, r_ymax
        IF (has_terminal) PRINT *, "From dns.in:"
        IF (has_terminal) PRINT *, nx, ny, nz, alfa0, beta0, ni, a, ymin, ymax
        STOP
      END IF
    ELSE
#endif
      IF (has_terminal) PRINT *, "Restart file "//filename//" not found"
      R = 0
      IF (has_terminal) WRITE (*, *) "Generating initial field..."
      DO iy = ny0 - 2, nyN + 2; DO ix = nx0, nxN; DO iz = -nz, nz
          CALL RANDOM_NUMBER(rn)
          !R(iy,iz,ix,1) = 0.0001*EXP(dcmplx(0,rn(1)-0.5));  R(iy,iz,ix,2) = 0.0001*EXP(dcmplx(0,rn(2)-0.5));  R(iy,iz,ix,3) = 0.0001*EXP(dcmplx(0,rn(3)-0.5));
        END DO; END DO; END DO
      IF (has_average) THEN
        DO iy = ny0 - 2, nyN + 2
          R(iy, 0, 0, 1) = 3*0.5*y(iy)*(2 - y(iy)) + 0.01*SIN(8*y(iy)*2*PI)/ni
          R(iy, 0, 0, 1) = y(iy)*(2 - y(iy))*3.d0/2.d0 + 0.001*SIN(8*y(iy)*2*PI); 
          !V(iy,0,0,1)=y(iy)-1
          DO iPhi = 1, nPhi
            R(iy, 0, 0, 3 + iPhi) = 3*0.5*y(iy)*(2 - y(iy))
          END DO
        END DO
      END IF
#ifdef HAVE_MPI
    END IF
#endif
    CLOSE (120)
  END SUBROUTINE read_restart_file

  !--------------------------------------------------------------!
  !-------------------- save_restart_file -----------------------!
  SUBROUTINE save_restart_file(filename, R)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: R(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3)
    character(len=*), intent(in) :: filename
    ! mpi stuff
#ifdef HAVE_MPI
    TYPE(MPI_File) :: fh
    INTEGER(MPI_OFFSET_KIND) :: disp
    TYPE(MPI_Status) :: status

    !$omp target update from(V)
    ! open file
    CALL MPI_File_open(MPI_COMM_WORLD, TRIM(filename), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)

    ! write header
    IF (has_terminal) THEN ! only one process does this
      CALL MPI_file_write(fh, [nx, ny, nz], 3, MPI_INTEGER, status)
      CALL MPI_file_write(fh, [alfa0, beta0, ni, a, ymin, ymax, time], 7, MPI_DOUBLE_PRECISION, status)
    END IF

    ! set view to subarray
    disp = 3*C_INT + 7*C_DOUBLE ! offset to skip header
    CALL MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, writeview_type, 'native', MPI_INFO_NULL)

    ! finally write field
    CALL MPI_File_write_all(fh, R, 1, owned2write_type, status)

    ! close file
    call MPI_File_close(fh)
#endif
  END SUBROUTINE save_restart_file

  !--------------------------------------------------------------!
  !------------------------- outstats ---------------------------!
  SUBROUTINE outstats()
    IMPLICIT NONE
    real(C_DOUBLE) :: runtime_global, dudy(1:2 + nPhi, 1:2)   !cfl
    character(len=40) :: istring, filename
    integer :: iPhi

    !$omp target update from(V(ny-3:ny+1, 0, 0, 1))
    !$omp target update from(V(ny-3:ny+1, 0, 0, 2))
    !$omp target update from(V(ny-3:ny+1, 0, 0, 3))
#ifdef HAVE_MPI
    CALL MPI_Allreduce(cfl, runtime_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD); cfl = 0; 
#else
    runtime_global = cfl
#endif
    IF (cflmax > 0) deltat = cflmax/runtime_global
    IF (has_average) THEN
      dudy(1, 2) = -sum(d14n(-2:2)*dreal(V(ny - 3:ny + 1, 0, 0, 1)))
      dudy(2, 2) = -sum(d14n(-2:2)*dreal(V(ny - 3:ny + 1, 0, 0, 3)))
      DO iPhi = 1, nPhi
        dudy(2 + iPhi, 2) = sum(d14n(-2:2)*dreal(V(ny - 3:ny + 1, 0, 0, 3 + iPhi)))
      END DO
      dudy(1, 1) = sum(d140(-2:2)*dreal(V(-1:3, 0, 0, 1)))
      dudy(2, 1) = sum(d140(-2:2)*dreal(V(-1:3, 0, 0, 3)))
      DO iPhi = 1, nPhi
        dudy(2 + iPhi, 1) = sum(d140(-2:2)*dreal(V(-1:3, 0, 0, 3 + iPhi)))
      END DO
    END IF
    IF (has_terminal) THEN
      !$omp target update from(fr)
      WRITE (*, "(F10.4,3X,4(F11.6,3X),4(F9.4,3X),2(F9.6,3X))") &
           time,dudy(1,1),dudy(1,2),dudy(2,1),dudy(2,2),fr(1)+corrpx*fr(3),meanpx+corrpx,fr(2)+corrpz*fr(3),meanpz+corrpz,runtime_global*deltat,deltat
      WRITE(121,*) time,dudy(1,1),dudy(1,2),dudy(2,1),dudy(2,2),fr(1)+corrpx*fr(3),meanpx+corrpx,fr(2)+corrpz*fr(3),meanpz+corrpz,runtime_global*deltat,deltat
      WRITE (195, *) time, dudy(3:, 1), dudy(3:, 2), fr(4:3 + nPhi) + corrtx(:)*fr(3 + nPhi + 1:3 + 2*nPhi), corrtx + meantx
      FLUSH (121); FLUSH (195)
    END IF
    runtime_global = 0
    !Save Dati.cart.out
    IF (((FLOOR((time + 0.5*deltat)/dt_save) > FLOOR((time - 0.5*deltat)/dt_save)) .AND. (istep > 1)) .OR. istep == nstep) THEN
      IF (has_terminal) WRITE (*, *) "Writing Dati.cart.out at time ", time
      filename = "Dati.cart.out"; CALL save_restart_file(filename, V)
    END IF
    !Save Dati.cart.i.out
    IF ((FLOOR((time + 0.5*deltat)/dt_field) > FLOOR((time - 0.5*deltat)/dt_field)) .AND. (istep > 1)) THEN
      ifield = ifield + 1; WRITE (istring, *) ifield
      IF (has_terminal) WRITE (*, *) "Writing Dati.cart."//TRIM(ADJUSTL(istring))//".out at time ", time
      filename = "Dati.cart."//TRIM(ADJUSTL(istring))//".out"; CALL save_restart_file(filename, V)
    END IF

  END SUBROUTINE outstats

END MODULE dnsdata
