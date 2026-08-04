#include "build_options.h"

! Wall boundary conditions -- both halves, and the only file you need to edit
! to change them.
!
! A wall condition is a *row* saying which combination of nodes is constrained,
! and a *value* saying what it is constrained to.  The two used to live apart:
! the rows here, the values written straight into bc0/bcn from inside the
! driver's time loop.  They are both here now.
!
!   setup_boundary_conditions  runs once, and builds the rows.
!   apply_wall_values          runs once per time step, and sets the values.
!
! Each component needs two rows at each wall: one enforcing the physical
! condition and one for the ghost node just outside it.  This module only
! builds them from the compact stencils; folding the ghost row into the wall
! row -- what leaves a banded system the wall-normal solver can take -- is
! done once at solve time, by the solver, for every component alike.  Doing it
! here as well would apply it twice.
!
!   v    : Dirichlet on v and on dv/dy (no-slip, impermeable)
!   eta  : Dirichlet on the wall-normal vorticity
!   phi  : Dirichlet by default, Neumann under -DphiNeumann
!
! Half-channel and body-force variants are selected by the cpp switches, as
! before.
!
! A time-dependent or wall-parallel-varying condition -- an oscillating wall,
! blowing and suction, a modulated scalar flux -- is a change to
! apply_wall_values alone.  It is called once per step, before the substeps,
! with the clock in `time`.
module channel_bcs

  use, intrinsic :: iso_c_binding
  use channel_grid
  use channel_state
  use stencil_coefficients
  use roctx, only: roctxPush, roctxPop

  implicit none

  public :: setup_boundary_conditions, apply_wall_values

contains

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

        IF (ix == 0 .AND. iz == 0) THEN
          bcn(iz, ix, 2) = 0
          bcn(iz, ix, 4) = 0
          bcn(iz, ix, 5) = dcmplx(dreal(bcn(iz, ix, 1)) - dimag(bcn(iz, ix, 3)), dimag(bcn(iz, ix, 1)) + dreal(bcn(iz, ix, 3)))
        ELSE
          bcn(iz, ix, 4) = -ialfa(ix)*bcn(iz, ix, 1) - ibeta(iz)*bcn(iz, ix, 3); bcn(iz, ix, 5) = ibeta(iz)*bcn(iz, ix, 1) - ialfa(ix)*bcn(iz, ix, 3)
        END IF
      END DO
    END DO
    !$omp target update from(bc0, bcn)
  END SUBROUTINE setup_boundary_conditions

  ! The inhomogeneous side of the wall conditions, refreshed once per time step.
  !
  ! Streamwise wall velocities u0/uN drive the Couette-like cases; the scalars
  ! are held at t0/tN on the mean mode and at zero on every other mode.  Both
  ! come from the input deck and are constant in time, so this recomputes the
  ! same planes every step -- which is what makes it the place to put a
  ! condition that is *not* constant.
  SUBROUTINE apply_wall_values()
    IMPLICIT NONE
    integer :: iPhi, ix, iz

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
  END SUBROUTINE apply_wall_values

end module channel_bcs
