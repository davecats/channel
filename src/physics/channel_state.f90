#include "header.h"

! The state of a running simulation: the fields, the boundary planes, the
! physical parameters and the clock.
!
! Declarations only, and deliberately no `use` of any other module.  That is a
! hard constraint, not a style choice: nvfortran fails to resolve a used
! module's declare-target variables inside downstream target regions when the
! module it comes from both `use`s something and declares its own
! declare-target variables.  A use-free declaration module is the shape that
! is known to work.  Anything needing `use` belongs elsewhere.
!
! Everything here is filled in by case_setup (read_dnsin, init_memory) and
! consumed by the solver, the statistics and the I/O.
module channel_state

  use, intrinsic :: iso_c_binding

  implicit none

  !-------------------------------------------------------------------------
  ! Solution fields.
  !
  ! V(iy,iz,ix,i)      i = 1:u, 2:v, 3:w, 4..3+nPhi: scalars
  ! oldrhs(iy,iz,ix,i) i = 1:eta, 2:d2v, 3..2+nPhi: scalars
  ! memrhs             storage the new right-hand side is built into
  !-------------------------------------------------------------------------
  complex(C_DOUBLE_COMPLEX), allocatable, target :: V(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: oldrhs(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), pointer :: memrhs(:, :, :, :)

  !-------------------------------------------------------------------------
  ! Boundary planes.  bc0/bcn(iz,ix,i) with i = 1:u, 2:v, 3:w, 4:vy, 5:eta,
  ! 6..5+nPhi: scalars.  zero_bc is a homogeneous plane passed where a
  ! component has no inhomogeneous condition.
  !-------------------------------------------------------------------------
  complex(C_DOUBLE_COMPLEX), allocatable, target :: bc0(:, :, :), bcn(:, :, :), zero_bc(:, :)

  ! Wall coefficient rows for each component: the row enforcing the condition
  ! at the wall, and the row for the ghost node just outside it.  Computed by
  ! setup_boundary_conditions from the compact stencils.
  real(C_DOUBLE), dimension(-2:2) :: v0bc, v0m1bc, vnbc, vnp1bc
  real(C_DOUBLE), dimension(-2:2) :: eta0bc, eta0m1bc, etanbc, etanp1bc
  real(C_DOUBLE), dimension(-2:2) :: phi0bc, phi0m1bc, phinbc, phinp1bc

  !-------------------------------------------------------------------------
  ! Physical parameters.  ni is the inverse Reynolds number; pra holds the
  ! inverse Prandtl number of each scalar.  u0/uN and t0/tN are the wall
  ! values driving Couette-like and scalar boundary conditions.
  !-------------------------------------------------------------------------
  real(C_DOUBLE) :: ni, gamma
  !$omp declare target(ni)
  real(C_DOUBLE) :: u0, uN, t0, tN
  !$omp declare target(u0, uN, t0, tN)
  real(C_DOUBLE), allocatable :: pra(:)
  real(C_DOUBLE) :: perturbation_amplitude = 5.54d-5

  !-------------------------------------------------------------------------
  ! Mean flow forcing and the corrections that hold the flow rate fixed.
  !-------------------------------------------------------------------------
  real(C_DOUBLE) :: meanpx, meanpz, meanflowx, meanflowz, meantx, meantb
  real(C_DOUBLE) :: corrpx = 0.d0, corrpz = 0.d0
  real(C_DOUBLE), allocatable :: corrtx(:)
  complex(C_DOUBLE_COMPLEX), allocatable :: ucor(:), tcor(:, :)
  ! Flow rates: 1:x, 2:z, 3:correction, then one per scalar and its correction.
  real(C_DOUBLE), allocatable :: fr(:)

  !-------------------------------------------------------------------------
  ! Clock and step control.
  !-------------------------------------------------------------------------
  real(C_DOUBLE) :: time, time0 = 0, deltat, cflmax, cfl = 0.0d0
  real(C_DOUBLE) :: dt_field, dt_save, t_max
  integer(C_SIZE_T) :: istep, nstep, ifield
  logical :: time_from_restart
  logical :: disable_restart_write = .false.
  logical :: overlapping

  !-------------------------------------------------------------------------
  ! Three-stage Runge-Kutta coefficients, one column per substep.
  !-------------------------------------------------------------------------
  real(C_DOUBLE), dimension(3, 3) :: RK_rai = reshape( &
                                     (/120.0d0/32.0d0, 2.0d0, 0.0d0, &
                                       120.0d0/8.0d0, 50.0d0/8.0d0, 34.0d0/8.0d0, &
                                       120.0d0/20.0d0, 90.0d0/20.0d0, 50.0d0/20.0d0/), &
                                     shape=(/3, 3/))

end module channel_state
