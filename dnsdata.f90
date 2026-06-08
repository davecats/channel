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
  USE config
  USE rbmat
  USE mpi_transpose
  USE roctx, only: roctxPush, roctxPop
  USE ffts

  IMPLICIT NONE

  !Simulation parameters
  real(C_DOUBLE) :: PI = 3.1415926535897932384626433832795028841971
  integer(C_INT) :: ny, nz, nxd, nPhi
  !$omp declare target(ny)
  real(C_DOUBLE) :: alfa0, beta0, ni, a, ymin, ymax, deltat, cflmax, time, time0 = 0, dt_field, dt_save, t_max, gamma
  !$omp declare target(ni)
  real(C_DOUBLE) :: u0, uN, t0, tN
  real(C_DOUBLE) :: meanpx, meanpz, meanflowx, meanflowz, meantx, meantb
  integer(C_INT), allocatable :: izd(:)
  complex(C_DOUBLE_COMPLEX), allocatable :: ialfa(:), ibeta(:)
  real(C_DOUBLE), allocatable :: k2(:, :)
  integer(C_INT) :: npy = 1
  logical :: time_from_restart
  logical :: disable_restart_write = .false.
  logical :: use_yslab_linsolve = .false.
#ifdef HAVE_CUDA
  complex(C_DOUBLE_COMPLEX), allocatable, save :: slab(:, :)
  integer(C_INT), save :: yslab_scratch_rows = -1
  integer(C_INT), save :: yslab_scratch_lines = -1
#endif
  !Grid
  integer(C_INT), private :: iy
  real(C_DOUBLE), allocatable :: y(:), dy(:)
  real(C_DOUBLE) :: dx, dz, factor
  !Derivatives
  real(C_DOUBLE), allocatable :: der(:, :, :)
  real(C_DOUBLE), dimension(-2:2) :: d040, d140, d240, d14m1, d24m1, d04n, d14n, d24n, d14np1, d24np1
  !$omp declare target(d14np1, d14n, d14m1, d140)
  real(C_DOUBLE), allocatable :: D0mat(:, :), eta00mat(:, :)
  real(C_DOUBLE), target, allocatable ::  ws(:)
  real(C_DOUBLE), pointer :: linsolve_mat(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), pointer :: memrhs(:, :, :, :)
#if !(defined(HAVE_CUDA) || defined(HAVE_HIP))
  !Fourier-transformable arrays (allocated in ffts.f90)
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:, :, :, :) :: VVdx, VVdz
  real(C_DOUBLE), pointer, dimension(:, :, :, :) :: rVVdx
#endif
  !Solution
  complex(C_DOUBLE_COMPLEX), allocatable :: oldrhs(:, :, :, :)
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
  logical :: overlapping

  public :: get_solver_memory_estimate, sync_velocity_to_device, apply_complex_derivative_current_layout

CONTAINS

  !--------------------------------------------------------------!
  !---------------------- Read input files ----------------------!
  SUBROUTINE read_dnsin(cfg)
    IMPLICIT NONE
    logical :: i
    integer :: iPhi
    type(ini_config), intent(in) :: cfg
    character(len=16) :: env_value
    integer(C_INT) :: nstep_in
    integer(C_INT) :: npy_in
    integer :: status, length
    logical :: found

    call require_integer(cfg, "mesh", "nx", nx)
    call require_integer(cfg, "mesh", "ny", ny)
    call require_integer(cfg, "mesh", "nz", nz)
    npy_in = 1
    call get_integer(cfg, "parallel", "npy", npy_in, found)
    if (.not. found) call get_integer(cfg, "mesh", "npy", npy_in, found)
    npy = npy_in
    call require_real(cfg, "mesh", "alfa0", alfa0)
    call require_real(cfg, "mesh", "beta0", beta0)
    nxd = 3*(nx + 1)/2
    nzd = 3*nz
    !$omp target update to(ny)
#ifdef useFFTfit
    i = fftFIT(nxd); DO WHILE (.NOT. i); nxd = nxd + 1; i = fftFIT(nxd); END DO
    i = fftFIT(nzd); DO WHILE (.NOT. i); nzd = nzd + 1; i = fftFIT(nzd); END DO
#endif

    call require_real(cfg, "velocity", "ni", ni)
    call require_real(cfg, "mesh", "stretching", a)
    call require_real(cfg, "mesh", "ymin", ymin)
    call require_real(cfg, "mesh", "ymax", ymax)
    ni = 1/ni
    call require_real(cfg, "velocity", "meanpx", meanpx)
    call require_real(cfg, "velocity", "meanpz", meanpz)
    call require_real(cfg, "velocity", "meanflowx", meanflowx)
    call require_real(cfg, "velocity", "meanflowz", meanflowz)
    call require_real(cfg, "velocity", "u0", u0)
    call require_real(cfg, "velocity", "un", uN)

    call require_integer(cfg, "scalars", "nphi", nPhi)
    call require_real(cfg, "scalars", "meantx", meantx)
    call require_real(cfg, "scalars", "meantb", meantb)
    call require_real(cfg, "scalars", "t0", t0)
    call require_real(cfg, "scalars", "tn", tN)
    allocate (pra(nPhi))
    if (nPhi > 0) then
      call require_real_vector(cfg, "scalars", "pr", pra(1:nPhi))
      DO iPhi = 1, nPhi
        pra(iPhi) = 1/pra(iPhi)
      END DO
    end if
    !$omp target enter data map(to: pra)

    call require_real(cfg, "timestepping", "deltat", deltat)
    call require_real(cfg, "timestepping", "cflmax", cflmax)
    call require_real(cfg, "timestepping", "time", time)
    call require_real(cfg, "timestepping", "dt_field", dt_field)
    call require_real(cfg, "timestepping", "dt_save", dt_save)
    call require_real(cfg, "timestepping", "t_max", t_max)
    call require_logical(cfg, "timestepping", "time_from_restart", time_from_restart)
    call require_integer(cfg, "timestepping", "nstep", nstep_in)
    nstep = int(nstep_in, C_SIZE_T)

    dx = PI/(alfa0*nxd); dz = 2.0d0*PI/(beta0*nzd); factor = 1.0d0/(2.0d0*nxd*nzd)
    call get_environment_variable("CHANNEL_OVERLAPPING", env_value, length, status)
    overlapping = .false.
    if (status == 0) then
      select case (adjustl(trim(env_value(:length))))
      case ("1", "true", "TRUE", "yes", "YES", "on", "ON")
        overlapping = .true.
      case ("0", "false", "FALSE", "no", "NO", "off", "OFF")
        overlapping = .false.
      case default
        print *, "Warning: invalid value for CHANNEL_OVERLAPPING:", trim(env_value(:length))
      end select
    end if
    call get_environment_variable("CHANNEL_DISABLE_RESTART_WRITE", env_value, length, status)
    disable_restart_write = .false.
    if (status == 0) then
      select case (adjustl(trim(env_value(:length))))
      case ("1", "true", "TRUE", "yes", "YES", "on", "ON")
        disable_restart_write = .true.
      case ("0", "false", "FALSE", "no", "NO", "off", "OFF")
        disable_restart_write = .false.
      case default
        print *, "Warning: invalid value for CHANNEL_DISABLE_RESTART_WRITE:", trim(env_value(:length))
      end select
    end if
    call get_environment_variable("CHANNEL_USE_YSLAB_LINSOLVE", env_value, length, status)
    use_yslab_linsolve = .false.
    if (status == 0) then
      select case (adjustl(trim(env_value(:length))))
      case ("1", "true", "TRUE", "yes", "YES", "on", "ON")
        use_yslab_linsolve = .true.
      case ("0", "false", "FALSE", "no", "NO", "off", "OFF")
        use_yslab_linsolve = .false.
      case default
        print *, "Warning: invalid value for CHANNEL_USE_YSLAB_LINSOLVE:", trim(env_value(:length))
      end select
    end if
  END SUBROUTINE read_dnsin

  !--------------------------------------------------------------!
  !---------------- Allocate memory for solution ----------------!
  SUBROUTINE init_memory(solveNS)
    use y_line_solvers, only: ys_prepare_ghost_field_workspace
    IMPLICIT NONE
    INTEGER(C_INT) :: ix, iz
#ifdef HAVE_CUDA
    integer(C_INT) :: yslab_first_line, yslab_line_count
#endif
    logical, intent(IN) :: solveNS
    ALLOCATE (V(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3 + nPhi)); V = 0
    !$omp target enter data map(to: V)
#ifdef bodyforce
    ALLOCATE (F(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3)); F = 0
    !$omp target enter data map(to: F)
#endif
    ALLOCATE (bc0(-nz:nz, nx0:nxN, 1:5 + nPhi), &
              bcn(-nz:nz, nx0:nxN, 1:5 + nPhi))
    bc0 = 0.0
    bcn = 0.0
    !$omp target enter data map(to: bc0, bcn)
    IF (solveNS) then
      ALLOCATE (memrhs(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:2 + nPhi), &
                oldrhs(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:2 + nPhi), &
                linsolve_mat(ny0:nyN + 2, -2:2, -nz:nz, nx0:nxN))
      memrhs = 0.0
      oldrhs = 0.0
      linsolve_mat = 0.0
      !$omp target enter data map(to: memrhs, linsolve_mat, oldrhs)
    END IF
#define newrhs(iy,iz,ix,i) memrhs(iy,iz,ix,i)
#define imod(iy) MOD(iy+1000,5)
    ALLOCATE (der(1:ny - 1, 0:3, -2:2), d0mat(ny0:nyN + 2, -2:2))
    ALLOCATE (eta00mat(ny0:nyN + 2, -2:2))
    der = 0.0; d0mat = 0.0; eta00mat = 0.0
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
    !$omp target enter data map(to: izd, ialfa, ibeta, k2, y, iy, rk_rai, ucor, tcor)
    IF (solveNS) OPEN (UNIT=195, FILE='Runtimedata.phi', ACTION='write')
    IF (solveNS .AND. has_terminal) OPEN (UNIT=121, FILE='Runtimedata', ACTION='write')

    allocate (fr(3 + 2*nPhi)); fr = 0.0
    call ys_prepare_ghost_field_workspace(ny, nz, nxB)
#ifdef HAVE_CUDA
    if (solveNS .and. use_yslab_linsolve) then
      call yslab_line_range(ipy, (nxN - nx0 + 1)*(2*nz + 1), yslab_first_line, yslab_line_count)
      call prepare_yslab_scratch(ny + 3, yslab_line_count)
    end if
#endif
  END SUBROUTINE init_memory

  SUBROUTINE get_solver_memory_estimate(solveNS, n_floats)
    IMPLICIT NONE
    logical, intent(in) :: solveNS
    integer(C_INT64_T), intent(out) :: n_floats
    integer(C_INT64_T) :: spectral_planes, bc_planes, linear_planes
    integer(C_INT64_T) :: sendcount64, nbufs

    spectral_planes = int(nyN - ny0 + 5, C_INT64_T)*int(2*nz + 1, C_INT64_T)*int(nxN - nx0 + 1, C_INT64_T)
    bc_planes = int(2*nz + 1, C_INT64_T)*int(nxN - nx0 + 1, C_INT64_T)*int(5 + nPhi, C_INT64_T)
    linear_planes = int(ny - 1, C_INT64_T)*int(2*nz + 1, C_INT64_T)*int(nxN - nx0 + 1, C_INT64_T)*int(2 + nPhi, C_INT64_T)
    sendcount64 = int(nxB, C_INT64_T)*int(nzB, C_INT64_T)*int(nyN - ny0 + 5, C_INT64_T)
    nbufs = int(merge(2, 1, overlapping), C_INT64_T)

    n_floats = 0_C_INT64_T
    n_floats = n_floats + 2_C_INT64_T*spectral_planes*int(3 + nPhi, C_INT64_T)
#ifdef bodyforce
    n_floats = n_floats + 2_C_INT64_T*spectral_planes*3_C_INT64_T
#endif
    n_floats = n_floats + 4_C_INT64_T*bc_planes
    if (solveNS) then
      n_floats = n_floats + 4_C_INT64_T*linear_planes
      n_floats = n_floats + 2_C_INT64_T*5_C_INT64_T*int(nyN - ny0 + 3, C_INT64_T)* &
                 int(2*nz + 1, C_INT64_T)*int(nxN - nx0 + 1, C_INT64_T)
    end if
    n_floats = n_floats + 4_C_INT64_T*sendcount64*int(nproc, C_INT64_T)*nbufs
  END SUBROUTINE get_solver_memory_estimate

  !--------------------------------------------------------------!
  !--------------- Deallocate memory for solution ---------------!
  SUBROUTINE free_memory(solveNS)
    use y_line_solvers, only: ys_release_ghost_field_workspace
    IMPLICIT NONE
    LOGICAL, intent(IN) :: solveNS
    !$omp target exit data map(delete: d240, d24m1, d04n, d24n, d24np1, D0mat)
    !$omp target exit data map(delete: V)
    !$omp target exit data map(delete: izd, ialfa, ibeta, k2, ucor, tcor)
#ifdef bodyforce
    !$omp target exit data map(delete: F)
#endif
    DEALLOCATE (V, der, d0mat, y, dy)
    !$omp target exit data map(delete: bc0, bcn)
    DEALLOCATE (bc0, bcn)
    IF (solveNS) THEN
      if (associated(memrhs) .or. associated(linsolve_mat) .or. allocated(oldrhs)) then
        !$omp target exit data map(delete: memrhs, oldrhs, linsolve_mat)
      end if
      if (associated(memrhs)) deallocate (memrhs)
      if (allocated(oldrhs)) deallocate (oldrhs)
      if (associated(linsolve_mat)) deallocate (linsolve_mat)
      CLOSE (UNIT=195)
      IF (has_terminal) CLOSE (UNIT=121)
    END IF
#ifdef HAVE_CUDA
    call release_yslab_scratch()
#endif
    call ys_release_ghost_field_workspace()
  END SUBROUTINE free_memory

  SUBROUTINE sync_velocity_to_device()
    IMPLICIT NONE

    !$omp target update to(V)
  END SUBROUTINE sync_velocity_to_device

  !--------------------------------------------------------------!
  !--------------- Set-up the compact derivatives ---------------!
  SUBROUTINE setup_derivatives()
    use y_line_solvers, only: ys_lu5decomp
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
    D0mat = 0.0d0
    FORALL (iy=max(1_C_INT, ny0):min(ny - 1, nyN)) D0mat(iy, -2:2) = der(iy, 0, -2:2)
    call ys_lu5decomp(D0mat)
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
    !$omp target update from(bc0, bcn)
  END SUBROUTINE setup_boundary_conditions

  !--------------------------------------------------------------!
  !---------------- integral in the y-direction -----------------!
#ifdef HAVE_CUDA
  !$omp declare target(yintegr)
#endif
  PURE FUNCTION yintegr(f, y) result(II)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: f(-1:ny + 1)
    real(C_DOUBLE), intent(in) :: y(-1:ny + 1)
    real(C_DOUBLE) :: II, yp1, ym1, a1, a2, a3
    integer(C_INT) :: iy
    II = 0.0d0
    DO iy = 1, ny, 2
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
  SUBROUTINE COMPLEXderiv(f0, f1, der, D0mat)
    use y_line_solvers, only: ys_solve_compact_derivative
    complex(C_DOUBLE_COMPLEX), intent(in)  :: f0(-1:ny + 1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: f1(-1:ny + 1)
    real(C_DOUBLE), DIMENSION(:, :), intent(in) :: der(ny0:nyN, 0:3, -2:2)
    real(C_DOUBLE), DIMENSION(:, :), intent(in) :: D0mat(ny0:nyN + 2, -2:2)

    call ys_solve_compact_derivative(f0, f1, der, D0mat, d140, d14m1, d14n, d14np1, ny, ny0, nyN)
  END SUBROUTINE COMPLEXderiv

  SUBROUTINE apply_complex_derivative_current_layout(src, dst, update_device)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dst(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    logical, optional, intent(in) :: update_device

#ifdef HAVE_CUDA
    if (use_yslab_linsolve) then
      call apply_complex_derivative_transposed_y(src, dst)
      return
    end if
#endif
    call assemble_compact_derivative_system(src)
    call solve_current_layout_derivative(dst)
  END SUBROUTINE apply_complex_derivative_current_layout

  subroutine assemble_compact_derivative_system(src)
    use y_line_solvers, only: ys_local_rhs, ys_local_operator, &
                              ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
                              ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: src(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer(C_INT) :: ix, iz, iy, iline, nlines_z, ix_first, ix_last, iz_first, iz_last, row_start, row_end

    nlines_z = 2*nz + 1
    ix_first = nx0
    ix_last = nxN
    iz_first = -nz
    iz_last = nz
    row_start = ny0
    row_end = nyN

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(src, der, d140, d14m1, d14n, d14np1, ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, &
    !$omp& ys_upper_ghost_rhs, ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ny, nz, nlines_z, ix_first, ix_last, &
    !$omp& iz_first, iz_last, row_start, row_end) private(ix, iz, iy, iline)
    do ix = ix_first, ix_last
      do iz = iz_first, iz_last
        iline = (ix - ix_first)*nlines_z + (iz - iz_first + 1)

        do iy = row_start, row_end
          ys_local_operator(iy, -2:2, iline) = der(iy, 0, -2:2)
          ys_local_rhs(iy, iline) = sum(der(iy, 1, -2:2)*src(iy - 2:iy + 2, iz, ix))
        end do

        if (row_start == 1) then
          ys_lower_ghost_row(:, iline) = 0.0d0
          ys_lower_boundary_row(:, iline) = 0.0d0
          ys_lower_boundary_rhs(iline) = sum(d140(-2:2)*src(-1:3, iz, ix))
          ys_lower_ghost_rhs(iline) = sum(d14m1(-2:2)*src(-1:3, iz, ix))
          ys_lower_ghost_row(-2, iline) = 1.0d0
          ys_lower_boundary_row(-1, iline) = 1.0d0
        end if
        if (row_end == ny - 1) then
          ys_upper_boundary_row(:, iline) = 0.0d0
          ys_upper_ghost_row(:, iline) = 0.0d0
          ys_upper_boundary_rhs(iline) = sum(d14n(-2:2)*src(ny - 3:ny + 1, iz, ix))
          ys_upper_ghost_rhs(iline) = sum(d14np1(-2:2)*src(ny - 3:ny + 1, iz, ix))
          ys_upper_boundary_row(1, iline) = 1.0d0
          ys_upper_ghost_row(2, iline) = 1.0d0
        end if
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine assemble_compact_derivative_system

  subroutine solve_current_layout_derivative(dst)
    use y_line_solvers, only: ys_solve_ghost_field_reduced
#ifdef HAVE_CUDA
    use y_line_solvers, only: ys_solve_ghost_field_reduced_const_operator
#endif
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)

#ifdef HAVE_CUDA
    if (npy_grid == 2) then
      call ys_solve_ghost_field_reduced_const_operator(dst, ny, nz)
    else
      call ys_solve_ghost_field_reduced(dst, ny, nz)
    end if
#else
    call ys_solve_ghost_field_reduced(dst, ny, nz)
#endif
  end subroutine solve_current_layout_derivative

  !$omp declare target(compact_boundary_rhs_value)
  subroutine compact_boundary_rhs_value(bc0_values, bcn_values, nz_value, nx0_value, iz, ix, rhs_index, value)
    implicit none
    integer(C_INT), intent(in) :: nz_value, nx0_value, iz, ix, rhs_index
    complex(C_DOUBLE_COMPLEX), intent(in) :: bc0_values(-nz_value:, nx0_value:, :)
    complex(C_DOUBLE_COMPLEX), intent(in) :: bcn_values(-nz_value:, nx0_value:, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: value

    if (rhs_index == 0_C_INT) then
      value = (0.0d0, 0.0d0)
    else if (rhs_index > 0_C_INT) then
      value = bc0_values(iz, ix, rhs_index)
    else
      value = bcn_values(iz, ix, -rhs_index)
    end if
  end subroutine compact_boundary_rhs_value

  !$omp declare target(compact_boundary_rhs_values)
  subroutine compact_boundary_rhs_values(bc0_values, bcn_values, nz_value, nx0_value, iz, ix, lower_rhs_index, lower_ghost_rhs_index, &
                                         upper_rhs_index, upper_ghost_rhs_index, lower_ghost_value, lower_boundary_value, &
                                         upper_boundary_value, upper_ghost_value)
    implicit none
    integer(C_INT), intent(in) :: nz_value, nx0_value, iz, ix
    integer(C_INT), intent(in) :: lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index
    complex(C_DOUBLE_COMPLEX), intent(in) :: bc0_values(-nz_value:, nx0_value:, :)
    complex(C_DOUBLE_COMPLEX), intent(in) :: bcn_values(-nz_value:, nx0_value:, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: lower_ghost_value, lower_boundary_value, upper_boundary_value, upper_ghost_value

    call compact_boundary_rhs_value(bc0_values, bcn_values, nz_value, nx0_value, iz, ix, lower_ghost_rhs_index, lower_ghost_value)
    call compact_boundary_rhs_value(bc0_values, bcn_values, nz_value, nx0_value, iz, ix, lower_rhs_index, lower_boundary_value)
    call compact_boundary_rhs_value(bc0_values, bcn_values, nz_value, nx0_value, iz, ix, upper_rhs_index, upper_boundary_value)
    call compact_boundary_rhs_value(bc0_values, bcn_values, nz_value, nx0_value, iz, ix, upper_ghost_rhs_index, upper_ghost_value)
  end subroutine compact_boundary_rhs_values

  !$omp declare target(compact_component_operator_coeffs)
subroutine compact_component_operator_coeffs(component_index, lambda_coeff, diffusion_coeff, ni_value, kk, der0, der2, der4, coeffs)
    implicit none
    integer(C_INT), intent(in) :: component_index
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff, ni_value, kk
    real(C_DOUBLE), intent(in) :: der0(-2:2), der2(-2:2), der4(-2:2)
    real(C_DOUBLE), intent(out) :: coeffs(-2:2)

    if (component_index == 2_C_INT) then
      coeffs = lambda_coeff*(der2 - kk*der0) - ni_value*(der4 - 2.0d0*kk*der2 + kk*kk*der0)
    else
      coeffs = lambda_coeff*der0 - diffusion_coeff*ni_value*(der2 - kk*der0)
    end if
  end subroutine compact_component_operator_coeffs

  !$omp declare target(compact_physical_boundary_equations)
  subroutine compact_physical_boundary_equations(lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, lower_boundary_value, &
                                                lower_ghost_value, upper_boundary_value, upper_ghost_value, lower_rhs0, lower_eq0, &
                                                 upper_rhsn, upper_eqn)
    implicit none
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: lower_boundary_value, lower_ghost_value, upper_boundary_value, upper_ghost_value
    complex(C_DOUBLE_COMPLEX), intent(out) :: lower_rhs0, upper_rhsn
    real(C_DOUBLE), intent(out) :: lower_eq0(-2:2), upper_eqn(-2:2)

    lower_rhs0 = lower_boundary_value - lower_ghost_value*lower_bc(-2)/lower_ghost_bc(-2)
    lower_eq0 = lower_bc - lower_ghost_bc*lower_bc(-2)/lower_ghost_bc(-2)
    lower_eq0(-2) = 0.0d0
    upper_rhsn = upper_boundary_value - upper_ghost_value*upper_bc(2)/upper_ghost_bc(2)
    upper_eqn = upper_bc - upper_ghost_bc*upper_bc(2)/upper_ghost_bc(2)
    upper_eqn(2) = 0.0d0
  end subroutine compact_physical_boundary_equations

  !$omp declare target(compact_eliminate_physical_boundary_row)
  subroutine compact_eliminate_physical_boundary_row(iy, ny_value, rhs_value, row_coeffs, lower_ghost_value, lower_rhs0, &
                                                lower_ghost_bc, lower_eq0, upper_ghost_value, upper_rhsn, upper_ghost_bc, upper_eqn)
    implicit none
    integer(C_INT), intent(in) :: iy, ny_value
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs_value
    real(C_DOUBLE), intent(inout) :: row_coeffs(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: lower_ghost_value, lower_rhs0, upper_ghost_value, upper_rhsn
    real(C_DOUBLE), intent(in) :: lower_ghost_bc(-2:2), lower_eq0(-2:2), upper_ghost_bc(-2:2), upper_eqn(-2:2)
    real(C_DOUBLE) :: fac

    if (iy == 1) then
      fac = row_coeffs(-2)/lower_ghost_bc(-2)
      rhs_value = rhs_value - lower_ghost_value*fac
      row_coeffs = row_coeffs - lower_ghost_bc*fac
      row_coeffs(-2) = 0.0d0
      fac = row_coeffs(-1)/lower_eq0(-1)
      rhs_value = rhs_value - lower_rhs0*fac
      row_coeffs = row_coeffs - lower_eq0*fac
      row_coeffs(-1) = 0.0d0
    else if (iy == 2) then
      fac = row_coeffs(-2)/lower_eq0(-1)
      rhs_value = rhs_value - lower_rhs0*fac
      row_coeffs(-2:1) = row_coeffs(-2:1) - lower_eq0(-1:2)*fac
      row_coeffs(-2) = 0.0d0
    else if (iy == ny_value - 2) then
      fac = row_coeffs(2)/upper_eqn(1)
      rhs_value = rhs_value - upper_rhsn*fac
      row_coeffs(-1:2) = row_coeffs(-1:2) - upper_eqn(-2:1)*fac
      row_coeffs(2) = 0.0d0
    else if (iy == ny_value - 1) then
      fac = row_coeffs(2)/upper_ghost_bc(2)
      rhs_value = rhs_value - upper_ghost_value*fac
      row_coeffs = row_coeffs - upper_ghost_bc*fac
      row_coeffs(2) = 0.0d0
      fac = row_coeffs(1)/upper_eqn(1)
      rhs_value = rhs_value - upper_rhsn*fac
      row_coeffs = row_coeffs - upper_eqn*fac
      row_coeffs(1) = 0.0d0
    end if
  end subroutine compact_eliminate_physical_boundary_row

  SUBROUTINE solve_compact_component_current_layout(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                                   lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                                    lambda_coeff, diffusion_coeff)
    IMPLICIT NONE
    integer(C_INT), intent(in) :: component_index, lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff

#ifdef HAVE_CUDA
    if (use_yslab_linsolve) then
      call solve_compact_component_transposed_y(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                                lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                                lambda_coeff, diffusion_coeff)
      return
    end if
    if (npy_grid == 1) then
      call solve_compact_component_local_full_y(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                                lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                                lambda_coeff, diffusion_coeff)
      return
    end if
#endif
    call assemble_compact_component_system(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                           lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                           lambda_coeff, diffusion_coeff)
    call solve_current_layout_schur(V(:, :, :, component_index))
  END SUBROUTINE solve_compact_component_current_layout

  subroutine assemble_compact_component_system(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                               lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                               lambda_coeff, diffusion_coeff)
    use y_line_solvers, only: ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, &
                            ys_upper_ghost_rhs, ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row
    implicit none
    integer(C_INT), intent(in) :: component_index, lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff
    real(C_DOUBLE) :: row_coeffs(-2:2)
    integer(C_INT) :: ix, iz, iy, iline, nlines_z, ix_first, ix_last, iz_first, iz_last, row_start, row_end

    nlines_z = 2*nz + 1
    ix_first = nx0
    ix_last = nxN
    iz_first = -nz
    iz_last = nz
    row_start = ny0
    row_end = nyN

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(V, der, k2, ni, lambda_coeff, diffusion_coeff, component_index, ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, &
    !$omp& ys_upper_boundary_rhs, ys_upper_ghost_rhs, ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, lower_bc, &
    !$omp& lower_ghost_bc, upper_bc, upper_ghost_bc, bc0, bcn, lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, nlines_z, ny, nz, &
    !$omp& ix_first, ix_last, iz_first, iz_last, row_start, row_end) private(ix, iz, iy, iline, row_coeffs)
    do ix = ix_first, ix_last
      do iz = iz_first, iz_last
        iline = (ix - ix_first)*nlines_z + (iz - iz_first + 1)

        ys_lower_ghost_row(:, iline) = lower_ghost_bc
        ys_lower_boundary_row(:, iline) = lower_bc
        ys_upper_boundary_row(:, iline) = upper_bc
        ys_upper_ghost_row(:, iline) = upper_ghost_bc

        do iy = row_start, row_end
          ys_local_rhs(iy, iline) = V(iy, iz, ix, component_index)
          call compact_component_operator_coeffs(component_index, lambda_coeff, diffusion_coeff, ni, k2(iz, ix), der(iy, 0, -2:2), &
                                                 der(iy, 2, -2:2), der(iy, 3, -2:2), row_coeffs)
          ys_local_operator(iy, -2:2, iline) = row_coeffs
        end do

        call compact_boundary_rhs_values(bc0, bcn, nz, ix_first, iz, ix, lower_rhs_index, lower_ghost_rhs_index, &
                                         upper_rhs_index, upper_ghost_rhs_index, ys_lower_ghost_rhs(iline), &
                                         ys_lower_boundary_rhs(iline), ys_upper_boundary_rhs(iline), ys_upper_ghost_rhs(iline))
      end do
    end do
    !$omp end target teams distribute parallel do
  end subroutine assemble_compact_component_system

  subroutine solve_current_layout_schur(dst)
    use y_line_solvers, only: ys_solve_ghost_field_reduced
#ifdef HAVE_CUDA
    use y_line_solvers, only: ys_solve_ghost_field_reduced_symmetric_operator
#endif
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)

#ifdef HAVE_CUDA
    if (npy_grid == 2) then
      call ys_solve_ghost_field_reduced_symmetric_operator(dst, ny, nz)
    else
      call ys_solve_ghost_field_reduced(dst, ny, nz)
    end if
#else
    call ys_solve_ghost_field_reduced(dst, ny, nz)
#endif
  end subroutine solve_current_layout_schur

#ifdef HAVE_CUDA
  subroutine prepare_yslab_scratch(nrows, nlines)
    implicit none
    integer(C_INT), intent(in) :: nrows, nlines

    if (nlines <= 0) return
    if (allocated(slab)) then
      if (yslab_scratch_rows /= nrows .or. yslab_scratch_lines /= nlines) then
        !$omp target exit data map(delete: slab)
        deallocate (slab)
        yslab_scratch_rows = -1
        yslab_scratch_lines = -1
      end if
    end if
    if (.not. allocated(slab)) then
      allocate (slab(nrows, nlines))
      !$omp target enter data map(alloc: slab)
      yslab_scratch_rows = nrows
      yslab_scratch_lines = nlines
    end if
  end subroutine prepare_yslab_scratch

  subroutine release_yslab_scratch()
    implicit none

    if (.not. allocated(slab)) return
    !$omp target exit data map(delete: slab)
    deallocate (slab)
    yslab_scratch_rows = -1
    yslab_scratch_lines = -1
  end subroutine release_yslab_scratch
#endif

#ifdef HAVE_CUDA
  subroutine solve_compact_component_local_full_y(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                                  lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                                  lambda_coeff, diffusion_coeff)
    implicit none
    integer(C_INT), intent(in) :: component_index, lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff
    integer(C_INT) :: nlines_z, nlines

    nlines_z = 2*nz + 1
    nlines = (nxN - nx0 + 1)*nlines_z
    call prepare_yslab_scratch(ny + 3, nlines)

    call roctxPush("compact full-y copy_to_slab")
    call yslab_copy_to_full(V(:, :, :, component_index), slab, ny, nz, 1_C_INT, nlines, nlines_z)
    call roctxPop("compact full-y copy_to_slab")

    call solve_full_y_lines_packed(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                   lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                   lambda_coeff, diffusion_coeff, 1_C_INT, nlines, "compact full-y gpsv")

    call roctxPush("compact full-y copy_from_slab")
    call yslab_copy_from_full(slab, V(:, :, :, component_index), ny, nz, 1_C_INT, nlines, nlines_z)
    call roctxPop("compact full-y copy_from_slab")
  end subroutine solve_compact_component_local_full_y

  subroutine solve_full_y_lines_packed(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                       lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                       lambda_coeff, diffusion_coeff, first_line, line_count, solve_label)
    use y_line_solvers, only: ys_prepare_gpsv_workspace, ys_solve_packed_gpsv, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, &
                              ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x
    implicit none
    integer(C_INT), intent(in) :: component_index, lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff
    integer(C_INT), intent(in) :: first_line, line_count
    character(len=*), intent(in) :: solve_label
    complex(C_DOUBLE_COMPLEX) :: rhs_value, lower_ghost_value, lower_boundary_value, upper_boundary_value, upper_ghost_value
    complex(C_DOUBLE_COMPLEX) :: lower_rhs0, upper_rhsn
    real(C_DOUBLE) :: lower_eq0(-2:2), upper_eqn(-2:2)
    real(C_DOUBLE) :: kk, fac, c_m2, c_m1, c_0, c_p1, c_p2
    integer(C_INT) :: nlines_z, ix_base, ilocal, iline, ix, iz, iy, p

    if (line_count <= 0) return

    nlines_z = 2*nz + 1
    ix_base = nx0
    call ys_prepare_gpsv_workspace(ny - 1, line_count)

    call roctxPush("compact full-y pack")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(slab, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, der, k2, ni, &
    !$omp& bc0, bcn, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, lambda_coeff, diffusion_coeff, component_index, &
    !$omp& lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, first_line, line_count, nlines_z, nz, ix_base, ny) &
    !$omp private(ilocal, iy, iline, ix, iz, p, kk, fac, c_m2, c_m1, c_0, c_p1, c_p2, rhs_value, lower_ghost_value, &
    !$omp& lower_boundary_value, upper_boundary_value, upper_ghost_value, lower_rhs0, upper_rhsn, lower_eq0, upper_eqn)
    do ilocal = 1, line_count
      do iy = 1, ny - 1
        iline = first_line + ilocal - 1
        ix = (iline - 1)/nlines_z + ix_base
        iz = mod(iline - 1, nlines_z) - nz

        call compact_boundary_rhs_value(bc0, bcn, nz, ix_base, iz, ix, lower_ghost_rhs_index, lower_ghost_value)
        call compact_boundary_rhs_value(bc0, bcn, nz, ix_base, iz, ix, lower_rhs_index, lower_boundary_value)
        call compact_boundary_rhs_value(bc0, bcn, nz, ix_base, iz, ix, upper_rhs_index, upper_boundary_value)
        call compact_boundary_rhs_value(bc0, bcn, nz, ix_base, iz, ix, upper_ghost_rhs_index, upper_ghost_value)

        lower_rhs0 = lower_boundary_value - lower_ghost_value*lower_bc(-2)/lower_ghost_bc(-2)
        lower_eq0 = lower_bc - lower_ghost_bc*lower_bc(-2)/lower_ghost_bc(-2)
        lower_eq0(-2) = 0.0d0
        upper_rhsn = upper_boundary_value - upper_ghost_value*upper_bc(2)/upper_ghost_bc(2)
        upper_eqn = upper_bc - upper_ghost_bc*upper_bc(2)/upper_ghost_bc(2)
        upper_eqn(2) = 0.0d0

        kk = k2(iz, ix)
        rhs_value = slab(iy + 2, ilocal)
        if (component_index == 2_C_INT) then
          c_m2 = lambda_coeff*(der(iy, 2, -2) - kk*der(iy, 0, -2)) - &
                 ni*(der(iy, 3, -2) - 2.0d0*kk*der(iy, 2, -2) + kk*kk*der(iy, 0, -2))
          c_m1 = lambda_coeff*(der(iy, 2, -1) - kk*der(iy, 0, -1)) - &
                 ni*(der(iy, 3, -1) - 2.0d0*kk*der(iy, 2, -1) + kk*kk*der(iy, 0, -1))
          c_0 = lambda_coeff*(der(iy, 2, 0) - kk*der(iy, 0, 0)) - &
                ni*(der(iy, 3, 0) - 2.0d0*kk*der(iy, 2, 0) + kk*kk*der(iy, 0, 0))
          c_p1 = lambda_coeff*(der(iy, 2, 1) - kk*der(iy, 0, 1)) - &
                 ni*(der(iy, 3, 1) - 2.0d0*kk*der(iy, 2, 1) + kk*kk*der(iy, 0, 1))
          c_p2 = lambda_coeff*(der(iy, 2, 2) - kk*der(iy, 0, 2)) - &
                 ni*(der(iy, 3, 2) - 2.0d0*kk*der(iy, 2, 2) + kk*kk*der(iy, 0, 2))
        else
          c_m2 = lambda_coeff*der(iy, 0, -2) - diffusion_coeff*ni*(der(iy, 2, -2) - kk*der(iy, 0, -2))
          c_m1 = lambda_coeff*der(iy, 0, -1) - diffusion_coeff*ni*(der(iy, 2, -1) - kk*der(iy, 0, -1))
          c_0 = lambda_coeff*der(iy, 0, 0) - diffusion_coeff*ni*(der(iy, 2, 0) - kk*der(iy, 0, 0))
          c_p1 = lambda_coeff*der(iy, 0, 1) - diffusion_coeff*ni*(der(iy, 2, 1) - kk*der(iy, 0, 1))
          c_p2 = lambda_coeff*der(iy, 0, 2) - diffusion_coeff*ni*(der(iy, 2, 2) - kk*der(iy, 0, 2))
        end if

        if (iy == 1) then
          fac = c_m2/lower_ghost_bc(-2)
          rhs_value = rhs_value - lower_ghost_value*fac
          c_m1 = c_m1 - lower_ghost_bc(-1)*fac
          c_0 = c_0 - lower_ghost_bc(0)*fac
          c_p1 = c_p1 - lower_ghost_bc(1)*fac
          c_p2 = c_p2 - lower_ghost_bc(2)*fac
          c_m2 = 0.0d0
          fac = c_m1/lower_eq0(-1)
          rhs_value = rhs_value - lower_rhs0*fac
          c_0 = c_0 - lower_eq0(0)*fac
          c_p1 = c_p1 - lower_eq0(1)*fac
          c_p2 = c_p2 - lower_eq0(2)*fac
          c_m1 = 0.0d0
        else if (iy == 2) then
          fac = c_m2/lower_eq0(-1)
          rhs_value = rhs_value - lower_rhs0*fac
          c_m1 = c_m1 - lower_eq0(0)*fac
          c_0 = c_0 - lower_eq0(1)*fac
          c_p1 = c_p1 - lower_eq0(2)*fac
          c_m2 = 0.0d0
        else if (iy == ny - 2) then
          fac = c_p2/upper_eqn(1)
          rhs_value = rhs_value - upper_rhsn*fac
          c_m1 = c_m1 - upper_eqn(-2)*fac
          c_0 = c_0 - upper_eqn(-1)*fac
          c_p1 = c_p1 - upper_eqn(0)*fac
          c_p2 = 0.0d0
        else if (iy == ny - 1) then
          fac = c_p2/upper_ghost_bc(2)
          rhs_value = rhs_value - upper_ghost_value*fac
          c_m2 = c_m2 - upper_ghost_bc(-2)*fac
          c_m1 = c_m1 - upper_ghost_bc(-1)*fac
          c_0 = c_0 - upper_ghost_bc(0)*fac
          c_p1 = c_p1 - upper_ghost_bc(1)*fac
          c_p2 = 0.0d0
          fac = c_p1/upper_eqn(1)
          rhs_value = rhs_value - upper_rhsn*fac
          c_m2 = c_m2 - upper_eqn(-2)*fac
          c_m1 = c_m1 - upper_eqn(-1)*fac
          c_0 = c_0 - upper_eqn(0)*fac
          c_p1 = 0.0d0
        end if

        p = (iy - 1)*line_count + ilocal
        ys_gpsv_x(p) = rhs_value
        ys_gpsv_ds(p) = cmplx(c_m2, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_dl(p) = cmplx(c_m1, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_d(p) = cmplx(c_0, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_du(p) = cmplx(c_p1, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_dw(p) = cmplx(c_p2, 0.0d0, kind=C_DOUBLE)
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("compact full-y pack")

    call ys_solve_packed_gpsv(ny - 1, line_count, solve_label)

    call roctxPush("compact full-y unpack")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(slab, ys_gpsv_x, der, k2, ni, bc0, bcn, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
    !$omp& lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, first_line, line_count, nlines_z, nz, ix_base, ny) &
    !$omp private(ilocal, iy, iline, ix, iz, p, lower_ghost_value, lower_boundary_value, upper_boundary_value, upper_ghost_value, &
    !$omp& lower_rhs0, upper_rhsn, lower_eq0, upper_eqn)
    do ilocal = 1, line_count
      iline = first_line + ilocal - 1
      ix = (iline - 1)/nlines_z + ix_base
      iz = mod(iline - 1, nlines_z) - nz
      do iy = 1, ny - 1
        p = (iy - 1)*line_count + ilocal
        slab(iy + 2, ilocal) = ys_gpsv_x(p)
      end do

      call compact_boundary_rhs_values(bc0, bcn, nz, ix_base, iz, ix, lower_rhs_index, lower_ghost_rhs_index, &
                                       upper_rhs_index, upper_ghost_rhs_index, lower_ghost_value, lower_boundary_value, &
                                       upper_boundary_value, upper_ghost_value)
      call compact_physical_boundary_equations(lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, lower_boundary_value, &
                                               lower_ghost_value, upper_boundary_value, upper_ghost_value, lower_rhs0, lower_eq0, &
                                               upper_rhsn, upper_eqn)

      slab(2, ilocal) = (lower_rhs0 - lower_eq0(0)*slab(3, ilocal) - lower_eq0(1)*slab(4, ilocal) - &
                         lower_eq0(2)*slab(5, ilocal))/lower_eq0(-1)
      slab(1, ilocal) = (lower_ghost_value - lower_ghost_bc(-1)*slab(2, ilocal) - lower_ghost_bc(0)*slab(3, ilocal) - &
                         lower_ghost_bc(1)*slab(4, ilocal) - lower_ghost_bc(2)*slab(5, ilocal))/lower_ghost_bc(-2)
      slab(ny + 2, ilocal) = (upper_rhsn - upper_eqn(-2)*slab(ny - 1, ilocal) - upper_eqn(-1)*slab(ny, ilocal) - &
                              upper_eqn(0)*slab(ny + 1, ilocal))/upper_eqn(1)
      slab(ny + 3, ilocal) = (upper_ghost_value - upper_ghost_bc(-2)*slab(ny - 1, ilocal) - &
                              upper_ghost_bc(-1)*slab(ny, ilocal) - upper_ghost_bc(0)*slab(ny + 1, ilocal) - &
                              upper_ghost_bc(1)*slab(ny + 2, ilocal))/upper_ghost_bc(2)
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("compact full-y unpack")
  end subroutine solve_full_y_lines_packed

  subroutine solve_compact_component_transposed_y(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                                  lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                                  lambda_coeff, diffusion_coeff)
    implicit none
    integer(C_INT), intent(in) :: component_index, lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff
    integer(C_INT) :: nlines_z, nlines, my_first_line, my_line_count

    nlines_z = 2*nz + 1
    nlines = (nxN - nx0 + 1)*nlines_z
    call yslab_line_range(ipy, nlines, my_first_line, my_line_count)
    if (my_line_count <= 0) return

    call prepare_yslab_scratch(ny + 3, my_line_count)

    call roctxPush("yslab compact transpose_to_full")
    call yslab_transpose_to_full(V(:, :, :, component_index), slab, ny, nz, nlines, nlines_z, .false.)
    call roctxPop("yslab compact transpose_to_full")

    call solve_full_y_lines_packed(component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                   lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                   lambda_coeff, diffusion_coeff, my_first_line, my_line_count, "yslab compact gpsv")

    call roctxPush("yslab compact transpose_from_full")
    call yslab_transpose_from_full(slab, V(:, :, :, component_index), ny, nz, nlines, nlines_z)
    call roctxPop("yslab compact transpose_from_full")

  end subroutine solve_compact_component_transposed_y

  subroutine apply_complex_derivative_transposed_y(src, dst)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: src(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(out) :: dst(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer(C_INT) :: nlines_z, nlines, my_first_line, my_line_count

    nlines_z = 2*nz + 1
    nlines = (nxN - nx0 + 1)*nlines_z
    call yslab_line_range(ipy, nlines, my_first_line, my_line_count)
    if (my_line_count <= 0) return

    call prepare_yslab_scratch(ny + 3, my_line_count)

    call roctxPush("yslab derivative transpose_to_full")
    call yslab_transpose_to_full(src, slab, ny, nz, nlines, nlines_z, .true.)
    call roctxPop("yslab derivative transpose_to_full")

    call solve_full_y_derivative_packed(my_line_count, "yslab derivative gpsv")

    call roctxPush("yslab derivative transpose_from_full")
    call yslab_transpose_from_full(slab, dst, ny, nz, nlines, nlines_z)
    call roctxPop("yslab derivative transpose_from_full")

  end subroutine apply_complex_derivative_transposed_y

  subroutine solve_full_y_derivative_packed(line_count, solve_label)
    use y_line_solvers, only: ys_prepare_gpsv_workspace, ys_solve_packed_gpsv, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, &
                              ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x
    implicit none
    integer(C_INT), intent(in) :: line_count
    character(len=*), intent(in) :: solve_label
    complex(C_DOUBLE_COMPLEX) :: rhs_value, lower_boundary_value, lower_ghost_value, upper_boundary_value, upper_ghost_value
    real(C_DOUBLE) :: c_m2, c_m1, c_0, c_p1, c_p2
    integer(C_INT) :: ilocal, iy, p

    if (line_count <= 0) return

    call ys_prepare_gpsv_workspace(ny - 1, line_count)

    call roctxPush("full-y derivative pack")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(slab, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, der, d140, d14m1, d14n, d14np1, &
    !$omp& line_count, ny) &
    !$omp private(ilocal, iy, p, rhs_value, c_m2, c_m1, c_0, c_p1, c_p2, lower_boundary_value, lower_ghost_value, upper_boundary_value, upper_ghost_value)
    do ilocal = 1, line_count
      do iy = 1, ny - 1
        lower_boundary_value = sum(d140(-2:2)*slab(1:5, ilocal))
        lower_ghost_value = sum(d14m1(-2:2)*slab(1:5, ilocal))
        upper_boundary_value = sum(d14n(-2:2)*slab(ny - 1:ny + 3, ilocal))
        upper_ghost_value = sum(d14np1(-2:2)*slab(ny - 1:ny + 3, ilocal))

        rhs_value = sum(der(iy, 1, -2:2)*slab(iy:iy + 4, ilocal))
        c_m2 = der(iy, 0, -2)
        c_m1 = der(iy, 0, -1)
        c_0 = der(iy, 0, 0)
        c_p1 = der(iy, 0, 1)
        c_p2 = der(iy, 0, 2)

        if (iy == 1) then
          rhs_value = rhs_value - c_m2*lower_ghost_value - c_m1*lower_boundary_value
          c_m2 = 0.0d0
          c_m1 = 0.0d0
        else if (iy == 2) then
          rhs_value = rhs_value - c_m2*lower_boundary_value
          c_m2 = 0.0d0
        else if (iy == ny - 2) then
          rhs_value = rhs_value - c_p2*upper_boundary_value
          c_p2 = 0.0d0
        else if (iy == ny - 1) then
          rhs_value = rhs_value - c_p1*upper_boundary_value - c_p2*upper_ghost_value
          c_p1 = 0.0d0
          c_p2 = 0.0d0
        end if

        p = (iy - 1)*line_count + ilocal
        ys_gpsv_x(p) = rhs_value
        ys_gpsv_ds(p) = cmplx(c_m2, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_dl(p) = cmplx(c_m1, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_d(p) = cmplx(c_0, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_du(p) = cmplx(c_p1, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_dw(p) = cmplx(c_p2, 0.0d0, kind=C_DOUBLE)
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("full-y derivative pack")

    call ys_solve_packed_gpsv(ny - 1, line_count, solve_label)

    call roctxPush("full-y derivative unpack")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(slab, ys_gpsv_x, d140, d14m1, d14n, d14np1, line_count, ny) &
    !$omp private(ilocal, iy, p, lower_boundary_value, lower_ghost_value, upper_boundary_value, upper_ghost_value)
    do ilocal = 1, line_count
      lower_boundary_value = sum(d140(-2:2)*slab(1:5, ilocal))
      lower_ghost_value = sum(d14m1(-2:2)*slab(1:5, ilocal))
      upper_boundary_value = sum(d14n(-2:2)*slab(ny - 1:ny + 3, ilocal))
      upper_ghost_value = sum(d14np1(-2:2)*slab(ny - 1:ny + 3, ilocal))
      do iy = 1, ny - 1
        p = (iy - 1)*line_count + ilocal
        slab(iy + 2, ilocal) = ys_gpsv_x(p)
      end do
      slab(1, ilocal) = lower_ghost_value
      slab(2, ilocal) = lower_boundary_value
      slab(ny + 2, ilocal) = upper_boundary_value
      slab(ny + 3, ilocal) = upper_ghost_value
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("full-y derivative unpack")
  end subroutine solve_full_y_derivative_packed
#endif

  SUBROUTINE scatter_full_y_line(full_line, local_line)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(in) :: full_line(-1:ny + 1)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: local_line(ny0 - 2:)
    integer(C_INT) :: local_y_first, local_y_last

    if (lbound(local_line, 1) <= -1 .and. ubound(local_line, 1) >= ny + 1) then
      local_line(-1:ny + 1) = full_line(-1:ny + 1)
      return
    end if

    local_y_first = lbound(local_line, 1)
    local_y_last = ubound(local_line, 1)
    local_line(local_y_first:local_y_last) = full_line(local_y_first:local_y_last)
  END SUBROUTINE scatter_full_y_line

  SUBROUTINE solve_mean_correction_line(x, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, lambda_coeff, diffusion_coeff)
    !!! THIS OPERATES ON A FULL Y LINE (ASSEMBLED FROM ALL PROCESSES) !!!
    use y_line_solvers, only: ys_solve_compact_system
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(-1:ny + 1)
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    real(C_DOUBLE), intent(in) :: lambda_coeff, diffusion_coeff
    real(C_DOUBLE) :: mat(1:ny + 1, -2:2)
    integer(C_INT) :: iy

    x(-1:0) = 0.0d0
    x(1:ny - 1) = 1.0d0
    x(ny:ny + 1) = 0.0d0
    mat = 0.0d0
    do iy = 1, ny - 1
      mat(iy, -2:2) = lambda_coeff*der(iy, 0, -2:2) - &
                      diffusion_coeff*ni*(der(iy, 2, -2:2) - k2(0, 0)*der(iy, 0, -2:2))
    end do
    call ys_solve_compact_system(x, mat, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                 (0.0d0, 0.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0), ny, 1_C_INT, ny - 1)
  END SUBROUTINE solve_mean_correction_line

! Orr-Sommerfeld and Squire opearators
#define OS(iy,j) (ni*(der(iy,3,j)-2.0d0*k2(iz,ix)*der(iy,2,j)+k2(iz,ix)*k2(iz,ix)*der(iy,0,j)))
#define SQ(iy,j) (ni*(der(iy,2,j)-k2(iz,ix)*der(iy,0,j)))

  SUBROUTINE linsolve(lambda)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: lambda
    integer(C_INT) :: ix, iz, i, j, iPhi, y_first, y_last
    complex(C_DOUBLE_COMPLEX) :: temp
    complex(C_DOUBLE_COMPLEX) :: zero_mode_u(-1:ny + 1), zero_mode_w(-1:ny + 1), zero_mode_ucor(-1:ny + 1)

    call roctxPush("linsolve solve_v")
    call solve_compact_component_current_layout(2_C_INT, v0bc, v0m1bc, vnbc, vnp1bc, &
                                                2_C_INT, 4_C_INT, -2_C_INT, -4_C_INT, lambda, 1.0d0)
    call roctxPop("linsolve solve_v")
    call roctxPush("linsolve solve_eta")
    call solve_compact_component_current_layout(1_C_INT, eta0bc, eta0m1bc, etanbc, etanp1bc, &
                                                5_C_INT, 0_C_INT, -5_C_INT, 0_C_INT, lambda, 1.0d0)
    call roctxPop("linsolve solve_eta")
    call roctxPush("linsolve d_v_dy")
    call apply_complex_derivative_current_layout(V(:, :, :, 2), V(:, :, :, 3))
    call roctxPop("linsolve d_v_dy")

    if (nx0 == 0) then
      call roctxPush("linsolve mean_correction")
      !$omp target update from(V(:, 0, 0, 1))
      call gather_full_y_line(ny, V(:, 0, 0, 1), zero_mode_u)
      zero_mode_w = dcmplx(dimag(zero_mode_u), 0.d0)
      zero_mode_u = dcmplx(dreal(zero_mode_u), 0.d0)
      call solve_mean_correction_line(zero_mode_ucor, eta0bc, eta0m1bc, etanbc, etanp1bc, lambda, 1.0d0)

      fr(1) = yintegr(zero_mode_u, y)
      fr(2) = yintegr(zero_mode_w, y)
      fr(3) = yintegr(zero_mode_ucor, y)
      IF (abs(meanflowx) > 1.0d-7) THEN
        corrpx = (meanflowx - fr(1))/fr(3)
        zero_mode_u = dcmplx(dreal(zero_mode_u) + corrpx*dreal(zero_mode_ucor), dimag(zero_mode_u))
      END IF
      IF (abs(meanflowz) > 1.0d-7) THEN
        corrpz = (meanflowz - fr(2))/fr(3)
        zero_mode_w = dcmplx(dreal(zero_mode_w) + corrpz*dreal(zero_mode_ucor), dimag(zero_mode_w))
      END IF
      call scatter_full_y_line(zero_mode_u, V(:, 0, 0, 1))
      call scatter_full_y_line(zero_mode_w, V(:, 0, 0, 3))
      !$omp target update to(V(:, 0, 0, 1), V(:, 0, 0, 3))
      call roctxPop("linsolve mean_correction")
    end if

    y_first = ny0
    y_last = nyN

    call roctxPush("linsolve recover_u_w")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ny, y_first, y_last, nz, nx0, nxN) &
    !$omp shared(V, k2, ialfa, ibeta) &
    !$omp private(ix, iz, j, temp)
    DO ix = nx0, nxN
      DO iz = -nz, nz
        if (ix == 0 .and. iz == 0) cycle
        do j = y_first - 2, y_last + 2
          temp = (ialfa(ix)*V(j, iz, ix, 3) - ibeta(iz)*V(j, iz, ix, 1))/k2(iz, ix)
          V(j, iz, ix, 3) = (ibeta(iz)*V(j, iz, ix, 3) + ialfa(ix)*V(j, iz, ix, 1))/k2(iz, ix)
          V(j, iz, ix, 1) = temp
        end do
      END DO
    END DO
    call roctxPop("linsolve recover_u_w")
  END SUBROUTINE linsolve

  SUBROUTINE linsolve_scalar(lambda, iPhi)
    IMPLICIT NONE
    integer(C_INT), intent(in) :: iPhi
    real(C_DOUBLE), intent(in) :: lambda
    integer(C_INT) :: ix, iz, i, j
    complex(C_DOUBLE_COMPLEX) :: temp
    complex(C_DOUBLE_COMPLEX) :: zero_mode_scalar(-1:ny + 1), zero_mode_tcor(-1:ny + 1)
    call solve_compact_component_current_layout(3_C_INT + iPhi, phi0bc, phi0m1bc, phinbc, phinp1bc, 5_C_INT + iPhi, 0_C_INT, &
                                                -(5_C_INT + iPhi), 0_C_INT, lambda, pra(iPhi))

    if (nx0 == 0) then
      !$omp target update from(V(:, 0, 0, 3 + iPhi))
      call gather_full_y_line(ny, V(:, 0, 0, 3 + iPhi), zero_mode_scalar)
      call solve_mean_correction_line(zero_mode_tcor, phi0bc, phi0m1bc, phinbc, phinp1bc, lambda, pra(iPhi))
      fr(3 + iPhi) = yintegr(zero_mode_scalar, y)
      fr(3 + nPhi + iPhi) = yintegr(zero_mode_tcor, y)
      IF (abs(meantb) > 1.0d-7) THEN
        corrtx(iPhi) = (meantb - fr(3 + iPhi))/fr(3 + nPhi + iPhi)
        zero_mode_scalar = dcmplx(dreal(zero_mode_scalar) + corrtx(iPhi)*dreal(zero_mode_tcor), dimag(zero_mode_scalar))
      END IF
      call scatter_full_y_line(zero_mode_scalar, V(:, 0, 0, 3 + iPhi))
      !$omp target update to(V(:, 0, 0, 3 + iPhi))
    end if

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ny, ny0, nyN, nz, nx0, nxN, nPhi, V) &
    !$omp private(ix, iz, temp)
    DO ix = nx0, nxN
      DO iz = -nz, nz
        if (ix == 0 .and. iz == 0) cycle
      END DO
    END DO
  END SUBROUTINE linsolve_scalar

  !--------------------------------------------------------------!
  !------------------------ convolutions ------------------------!

  SUBROUTINE assemble_vvdz(m, to)
    IMPLICIT NONE
    integer(C_INT) :: iy
    integer(C_INT) :: i, j, k, y_first, y_last
    integer(C_INT), intent(in) :: m, to
    y_first = ny0
    y_last = nyN

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, VVdz) shared(nxB, nzd, nx0, nxN, y_first, y_last, nz, m, to) private(i,j,k)
    DO i = y_first - 2, y_last + 2
      DO j = 1, nxB
        DO k = 1, nzd
          VVdz(k, j, i, to) = 0.0d0
        END DO
      END DO
    END DO
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, VVdz) shared(nx0, nxN, nz, nzd, y_first, y_last, m, to) private(i,j,k)
    DO i = y_first - 2, y_last + 2
      DO j = nx0, nxN
        DO k = 1, nzd
          IF (k <= nz + 1) THEN
            VVdz(k, j - nx0 + 1, i, to) = V(i, k - 1, j, m)
          ELSEIF (k >= nz + 2 .AND. k <= nzd - nz) THEN
            VVdz(k, j - nx0 + 1, i, to) = 0.0d0
          ELSE
            VVdz(k, j - nx0 + 1, i, to) = V(i, k - nzd - 1, j, m)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE assemble_vvdz

  SUBROUTINE zero_vvdx_hft(to)
    IMPLICIT NONE
    integer(C_INT) :: iy
    integer(C_INT) :: i, j, k, y_first, y_last
    integer(C_INT), intent(in) :: to
    y_first = ny0
    y_last = nyN
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdx) shared(nx, nxd, nzB, y_first, y_last, to) private(i,j,k)
    DO i = y_first - 2, y_last + 2
      DO j = 1, nzB
        DO k = nx + 2, nxd + 1
          VVdx(k, j, i, to) = 0.0
        END DO
      END DO
    END DO
  END SUBROUTINE zero_vvdx_hft

  subroutine compute_cfl()
    implicit none
    integer(C_INT) :: i, j, k, y_first, y_last
    real(C_DOUBLE) :: tmp
    y_first = ny0
    y_last = nyN
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp private(i,j,k,tmp) reduction(max:cfl) &
    !$omp shared(rVVdx, dx, dy, dz, ny, y_first, y_last, nxd, nzB)
    do j = 1, 2*nxd
      do k = 1, nzB
        do i = max(y_first - 2, 1_C_INT), min(y_last + 2, ny - 1)
          tmp = abs(rVVdx(j, k, i, 1))/dx + abs(rVVdx(j, k, i, 2))/dy(i) + abs(rVVdx(j, k, i, 3))/dz
          cfl = max(cfl, tmp)
        end do
      end do
    end do
  END SUBROUTINE compute_cfl

  SUBROUTINE build_products(m, to)
    IMPLICIT NONE
    integer(C_INT), intent(in) :: m, to
    integer(C_INT) :: iPhi, component, first, second, y_first, y_last
    integer(C_INT) :: i, j, k
    real(C_DOUBLE) :: a, b
    y_first = ny0
    y_last = nyN

    if (m > 6) then ! scalar cases 7, ...
      iPhi = (m - 4)/3
      component = mod(m - 4, 3) + 1
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(rVVdx, products) shared(nxd, nzB, y_first, y_last, factor, iPhi, component, to)
      DO i = y_first - 2, y_last + 2
        DO j = 1, nzB
          DO k = 1, 2*nxd
            products(k, j, i, to) = rVVdx(k, j, i, component)*rVVdx(k, j, i, 3 + iPhi)*factor
          END DO
        END DO
      END DO
    else if (m > 3) then !cases 4, 5, 6
      first = mod(m - 1, 3) + 1
      second = mod(m, 3) + 1
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(rVVdx, products) shared(nxd, nzB, y_first, y_last, factor, first, second, m, to) private(a, b)
      DO i = y_first - 2, y_last + 2
        DO j = 1, nzB
          DO k = 1, 2*nxd
            a = rVVdx(k, j, i, first)
            b = rVVdx(k, j, i, second)
            products(k, j, i, to) = a*b*factor
          END DO
        END DO
      END DO
    else ! cases 1, 2, 3
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(rVVdx, products) shared(nxd, nzB, y_first, y_last, factor, m, to) private(a)
      DO i = y_first - 2, y_last + 2
        DO j = 1, nzB
          DO k = 1, 2*nxd
            a = rVVdx(k, j, i, m)
            products(k, j, i, to) = a*a*factor
          END DO
        END DO
      END DO
    end if
  END SUBROUTINE build_products

  SUBROUTINE transform_to_physical()
    IMPLICIT NONE
    integer(C_INT) ::  m, to, from, mm1
    type(MPI_Request), dimension(:) :: requests(3 + nPhi)
    type(MPI_Status)  :: status
    DO m = 1, MERGE(3 + nPhi + 1, 3 + nPhi, overlapping)

      ! Compute indices depending on overlap mode
      mm1 = MERGE(m - 1, m, overlapping)
      to = MERGE(mod(m - 1, 2) + 1, 1, overlapping)
      from = MERGE(mod(mm1 - 1, 2) + 1, 1, overlapping)

      ! Step 1: assemble, pack, post alltoall (only if in range)
      if (m <= 3 + nPhi) then
        call roctxPush("transform_to_physical assemble_vvdz")
        CALL assemble_vvdz(m, to)
        call roctxPop("transform_to_physical assemble_vvdz")
        call roctxPush("transform_to_physical IFT")
        CALL IFT(VVdz(:, :, :, to))
        call roctxPop("transform_to_physical IFT")
        if (fft_transpose_is_local) then
          call roctxPush("transform_to_physical repack_zTOx_local")
          call repack_zTOx_local(VVdz(:, :, :, to), VVdx(:, :, :, to), ny)
          call roctxPop("transform_to_physical repack_zTOx_local")
        else
          call roctxPush("transform_to_physical pack_zTOx")
          CALL pack_zTOx(VVdz(:, :, :, to), sendbuf(:, to), ny)
          call roctxPop("transform_to_physical pack_zTOx")
          CALL alltoall(sendbuf(:, to), recvbuf(:, to), requests(m), "zTOx transform_to_physical")
        end if
      end if

      ! Step 2: wait, unpack, FFT (depending on overlap)
      if (MERGE(m > 1, .true., overlapping)) then
        if (.not. fft_transpose_is_local) then
          call roctxPush("MPI_Wait zTOx transform_to_physical")
          CALL MPI_WAIT(requests(mm1), status, ierr)
          call roctxPop("MPI_Wait zTOx transform_to_physical")
          call roctxPush("transform_to_physical unpack_zTOx")
          CALL unpack_zTOx(recvbuf(:, from), VVdx(:, :, :, from), ny)
          call roctxPop("transform_to_physical unpack_zTOx")
        end if
        call roctxPush("transform_to_physical zero_vvdx_hft")
        CALL zero_vvdx_hft(from)
        call roctxPop("transform_to_physical zero_vvdx_hft")
        call roctxPush("transform_to_physical RFT")
        CALL RFT(VVdx(:, :, :, from), rVVdx(:, :, :, mm1))
        call roctxPop("transform_to_physical RFT")
      end if
    END DO
  END SUBROUTINE transform_to_physical

  SUBROUTINE transform_back_and_build_rhs(ODE)
    IMPLICIT NONE
    integer(C_INT) ::  m, mm1, to, from
    real(C_DOUBLE), intent(in) :: ODE(1:3)
    type(MPI_Request), dimension(:) :: requests(6 + 3*nPhi)
    type(MPI_Status)  :: status

    ! Reverse pass to build rVVdx
    DO m = 1, MERGE(6 + 3*nPhi + 1, 6 + 3*nPhi, overlapping)

      ! Buffer indices depend on overlap mode

      mm1 = MERGE(m - 1, m, overlapping)
      to = MERGE(mod(m - 1, 2) + 1, 1, overlapping)
      from = MERGE(mod(m, 2) + 1, 1, overlapping)

      ! Step 1: Build, HFT, pack, and post alltoall
      if (m <= 6 + 3*nPhi) then
        call roctxPush("transform_back build_products")
        call build_products(m, to)
        call roctxPop("transform_back build_products")
        call roctxPush("transform_back HFT")
        call HFT(products(:, :, :, to), VVdx(:, :, :, to))
        call roctxPop("transform_back HFT")
        if (fft_transpose_is_local) then
          call roctxPush("transform_back repack_xTOz_local")
          call repack_xTOz_local(VVdx(:, :, :, to), VVdz(:, :, :, to), ny)
          call roctxPop("transform_back repack_xTOz_local")
        else
          call roctxPush("transform_back pack_xTOz")
          call pack_xTOz(VVdx(:, :, :, to), sendbuf(:, to), ny)
          call roctxPop("transform_back pack_xTOz")
          call alltoall(sendbuf(:, to), recvbuf(:, to), requests(m), "xTOz transform_back_and_build_rhs")
        end if
      end if

      ! Step 2: Wait, unpack, FFT, and build RHS
      if (MERGE(m > 1, .true., overlapping)) then
        if (.not. fft_transpose_is_local) then
          call roctxPush("MPI_Wait xTOz transform_back_and_build_rhs")
          call MPI_WAIT(requests(mm1), status, ierr)
          call roctxPop("MPI_Wait xTOz transform_back_and_build_rhs")
          call roctxPush("transform_back unpack_xTOz")
          call unpack_xTOz(recvbuf(:, from), VVdz(:, :, :, from), ny)
          call roctxPop("transform_back unpack_xTOz")
        end if
        call roctxPush("transform_back FFT")
        call FFT(VVdz(:, :, :, from))
        call roctxPop("transform_back FFT")
        call roctxPush("transform_back buildrhs")
        call buildrhs(ODE, mm1, from)
        call roctxPop("transform_back buildrhs")
      end if
    END DO
  END SUBROUTINE transform_back_and_build_rhs

  !--------------------------------------------------------------!
  !-------------------------- buildRHS --------------------------!
  ! (u,v,w) = (1,2,3)
  ! (uu,vv,ww,uv,vw,uw) = (1,2,3,4,5,6)
#define timescheme(rhs,old,unkn,impl,expl) rhs=ODE(1)*(unkn)/deltat+(impl)+ODE(2)*(expl)-ODE(3)*(old); old=expl
#define timescheme_accum(rhs,old,expl) rhs = rhs + ODE(2)*(expl); old = old + expl
#define DD(f,k) ( der(iy,f,-2)*VVdz(izd(iz)+1,ix+1-nx0,iy-2,k)+ \
  der(iy, f, -1)*VVdz(izd(iz) + 1, ix + 1 - nx0, iy - 1, k) + \
  der(iy, f, 0)*VVdz(izd(iz) + 1, ix + 1 - nx0, iy, k) + \
  der(iy, f, 1)*VVdz(izd(iz) + 1, ix + 1 - nx0, iy + 1, k) + \
  der(iy, f, 2)*VVdz(izd(iz) + 1, ix + 1 - nx0, iy + 2, k))
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
    integer(C_INT) :: iy, iz, ix, i, im2, im1, i0, i1, i2, k, iPhi, y_first, y_last
    complex(C_DOUBLE_COMPLEX) :: rhsu, rhsw, rhst, expl, tmp, unkn
    y_first = ny0
    y_last = nyN
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
    !$omp shared(nz, nx0, nxN, ny, y_first, y_last) shared(ialfa, ibeta) shared(k2, der) shared(memrhs, oldrhs) &
    !$omp shared(meanpx, meanpz, ni, deltat, ode) shared(vvdz, v, pra)
    DO iz = -nz, nz
      DO ix = nx0, nxN
        DO iy = y_first, y_last
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
    !$omp shared(nz, nx0, nxN, ny, y_first, y_last) shared(ialfa, ibeta) shared(k2, der) shared(memrhs, oldrhs) &
    !$omp shared(meanpx, meanpz, ni, deltat, ode) shared(vvdz, v, pra)
    DO iz = -nz, nz
      DO ix = nx0, nxN
        DO iy = y_first, y_last
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
    !$omp shared(nz, nx0, nxN, ny, y_first, y_last) shared(ialfa, ibeta) shared(k2, der) shared(memrhs, oldrhs) &
    !$omp shared(meanpx, meanpz, ni, deltat, ode) shared(vvdz, v, pra)
    DO iz = -nz, nz
      DO ix = nx0, nxN
        DO iy = y_first, y_last

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
      !$omp shared(nz, nx0, nxN, ny, y_first, y_last, pra, V, memrhs, ode, deltat, oldrhs, iphi, der, k2, ni)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = y_first, y_last
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
    !$omp shared(memrhs, V) shared(nz, nx0, nxN, ny, y_first, y_last, nPhi) private(iy, ix, iz, iPhi)
    DO iz = -nz, nz
    DO ix = nx0, nxN
    DO iy = y_first, y_last
      V(iy, iz, ix, 1) = newrhs(iy, iz, ix, 1); 
      V(iy, iz, ix, 2) = newrhs(iy, iz, ix, 2); 
      DO iPhi = 1, nPhi
        V(iy, iz, ix, 3 + iPhi) = newrhs(iy, iz, ix, 2 + iPhi); 
      END DO
    END DO
    END DO
    END DO

  END SUBROUTINE buildrhs_prepare

  SUBROUTINE buildrhs(ODE, component, from)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: ODE(1:3)
    integer(C_INT), intent(in) :: component, from
    integer(C_INT) :: iy, iz, ix, i, iPhi, y_first, y_last
    complex(C_DOUBLE_COMPLEX) :: rhsu, rhsw, rhst, expl, tmp, unkn
    y_first = ny0
    y_last = nyN

    SELECT CASE (component)
    CASE (1)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny, y_first, y_last, from) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = y_first, y_last

            rhsu = -ialfa(ix)*DD(0, from); rhsw = 0.0
            expl = ialfa(ix)*ialfa(ix)*DD(1, from)

            ACCUM_D2V(V, expl)
            ACCUM_ETA_UW(V, rhsu, rhsw)
          END DO
        END DO
      END DO
    CASE (2)
      !contribution from VVdz(:,:,:,2)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny, y_first, y_last, from) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = y_first, y_last

            expl = DD(1, from)*k2(iz, ix)
            ACCUM_D2V(V, expl)
          END DO
        END DO
      END DO
    CASE (3)
      !contribution from VVdz(:,:,:,3)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny, y_first, y_last, from) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = y_first, y_last
            rhsw = -ibeta(iz)*DD(0, from); rhsu = 0.0
            expl = ibeta(iz)*ibeta(iz)*DD(1, from)

            ACCUM_D2V(V, expl)
            ACCUM_ETA_UW(V, rhsu, rhsw)
          END DO
        END DO
      END DO
    CASE (4)
      !contribution from VVdz(:,:,:,4)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny, y_first, y_last, from) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = y_first, y_last

            rhsu = -DD(1, from); rhsw = 0.0
            expl = ialfa(ix)*DD(2, from) + ialfa(ix)*DD(0, from)*k2(iz, ix)

            ACCUM_D2V(V, expl)
            ACCUM_ETA_UW(V, rhsu, rhsw)
          END DO
        END DO
      END DO
    CASE (5)
      !contribution from VVdz(:,:,:,5)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny, y_first, y_last, from) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = y_first, y_last

            ! contribution from VVdz(:,:,:,5)
            rhsw = -DD(1, from); rhsu = 0.0
            expl = ibeta(iz)*DD(2, from) + ibeta(iz)*DD(0, from)*k2(iz, ix)

            ACCUM_D2V(V, expl)
            ACCUM_ETA_UW(V, rhsu, rhsw)
          END DO
        END DO
      END DO
    CASE (6)
      !contribution from VVdz(:,:,:,6)
      !$omp target teams distribute parallel do collapse(3) default(none)  &
      !$omp private(rhsu, rhsw, expl) private(iz, ix, iy) &
      !$omp shared(nz, nx0, nxN, ny, y_first, y_last, from) shared(ialfa, ibeta) shared(der) shared(V, oldrhs, izd, ode, k2) shared(vvdz)
      DO iz = -nz, nz
        DO ix = nx0, nxN
          DO iy = y_first, y_last

            rhsu = -ibeta(iz)*DD(0, from)
            rhsw = -ialfa(ix)*DD(0, from)
            expl = 2*ialfa(ix)*ibeta(iz)*DD(1, from)

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
        !$omp shared(nz, nx0, nxN, ny, y_first, y_last, ialfa, V, ode, der, oldrhs, iphi, vvdz, izd, from)
        DO iz = -nz, nz
          DO ix = nx0, nxN
            DO iy = y_first, y_last
              rhst = -ialfa(ix)*DD(0, from)
              timescheme_accum(V(iy, iz, ix, 3 + iPhi), oldrhs(iy, iz, ix, 2 + iPhi), rhst)
            END DO
          END DO
        END DO
      CASE (2)
        !$omp target teams distribute parallel do collapse(3) default(none)  &
        !$omp private(rhst, iz, ix, iy) &
        !$omp shared(nz, nx0, nxN, ny, y_first, y_last, V, ode, der, oldrhs, iphi, vvdz, izd, from)
        DO iz = -nz, nz
          DO ix = nx0, nxN
            DO iy = y_first, y_last

              rhst = -DD(1, from)
              timescheme_accum(V(iy, iz, ix, 3 + iPhi), oldrhs(iy, iz, ix, 2 + iPhi), rhst)
            END DO
          END DO
        END DO
      CASE (3)
        !$omp target teams distribute parallel do collapse(3) default(none)  &
        !$omp private(rhst, iz, ix, iy) &
        !$omp shared(nz, nx0, nxN, ny, y_first, y_last, ibeta,V, ode, der, oldrhs, iphi, vvdz, izd, from)
        DO iz = -nz, nz
          DO ix = nx0, nxN
            DO iy = y_first, y_last
              rhst = -ibeta(iz)*DD(0, from)
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
      call roctxPush("MPI_File_read_all restart")
      call MPI_file_read_all(fh, R, 1, vel_field_type, MPI_STATUS_IGNORE)
      call roctxPop("MPI_File_read_all restart")
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
          R(iy, iz, ix, 1) = 0.0000554*EXP(dcmplx(0, rn(1) - 0.5)); R(iy, iz, ix, 2) = 0.0000554*EXP(dcmplx(0, rn(2) - 0.5)); R(iy, iz, ix, 3) = 0.0000554*EXP(dcmplx(0, rn(3) - 0.5)); 
          !!R(iy,iz,ix,1) = 0.0001*EXP(dcmplx(0,rn(1)-0.5));  R(iy,iz,ix,2) = 0.0001*EXP(dcmplx(0,rn(2)-0.5));  R(iy,iz,ix,3) = 0.0001*EXP(dcmplx(0,rn(3)-0.5));
        END DO; END DO; END DO
      IF (has_average) THEN
        DO iy = ny0 - 2, nyN + 2
          R(iy, 0, 0, 1) = 3*0.5*y(iy)*(2 - y(iy))
          !R(iy, 0, 0, 1) = 3*0.5*y(iy)*(2 - y(iy)) + 0.01*SIN(8*y(iy)*2*PI)/ni
          !R(iy, 0, 0, 1) = y(iy)*(2 - y(iy))*3.d0/2.d0 + 0.001*SIN(8*y(iy)*2*PI);
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
    complex(C_DOUBLE_COMPLEX), intent(in) :: R(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3 + nPhi)
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
    call roctxPush("MPI_File_write_all restart")
    CALL MPI_File_write_all(fh, R, 1, owned2write_type, status)
    call roctxPop("MPI_File_write_all restart")

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
    complex(C_DOUBLE_COMPLEX) :: mean_line_u(-1:ny + 1), mean_line_w(-1:ny + 1), mean_line_scalar(-1:ny + 1)

    if (has_average) then
      !$omp target update from(V(ny0 - 2:nyN + 2, 0, 0, 1))
      !$omp target update from(V(ny0 - 2:nyN + 2, 0, 0, 3))
      call gather_full_y_line(ny, V(:, 0, 0, 1), mean_line_u)
      call gather_full_y_line(ny, V(:, 0, 0, 3), mean_line_w)
    end if
#ifdef HAVE_MPI
    call roctxPush("MPI_Allreduce outstats_cfl")
    CALL MPI_Allreduce(cfl, runtime_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD); cfl = 0; 
    call roctxPop("MPI_Allreduce outstats_cfl")
#else
    runtime_global = cfl
#endif
    IF (cflmax > 0) deltat = cflmax/runtime_global
    IF (has_average) THEN
      dudy(1, 2) = -sum(d14n(-2:2)*dreal(mean_line_u(ny - 3:ny + 1)))
      dudy(2, 2) = -sum(d14n(-2:2)*dreal(mean_line_w(ny - 3:ny + 1)))
      DO iPhi = 1, nPhi
        !$omp target update from(V(ny0 - 2:nyN + 2, 0, 0, 3 + iPhi))
        call gather_full_y_line(ny, V(:, 0, 0, 3 + iPhi), mean_line_scalar)
        dudy(2 + iPhi, 2) = sum(d14n(-2:2)*dreal(mean_line_scalar(ny - 3:ny + 1)))
        dudy(2 + iPhi, 1) = sum(d140(-2:2)*dreal(mean_line_scalar(-1:3)))
      END DO
      dudy(1, 1) = sum(d140(-2:2)*dreal(mean_line_u(-1:3)))
      dudy(2, 1) = sum(d140(-2:2)*dreal(mean_line_w(-1:3)))
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
    IF (.not. disable_restart_write .and. &
        (((FLOOR((time + 0.5*deltat)/dt_save) > FLOOR((time - 0.5*deltat)/dt_save)) .AND. (istep > 1)) .OR. istep == nstep)) THEN
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
