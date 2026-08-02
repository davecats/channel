#include "header.h"

! Wall boundary conditions.
!
! Each component needs two rows at each wall: one enforcing the physical
! condition and one for the ghost node just outside it.  They are built from
! the compact stencils and then combined so the ghost unknown is eliminated,
! which is what leaves a banded system the wall-normal solver can take.
!
!   v    : Dirichlet on v and on dv/dy (no-slip, impermeable)
!   eta  : Dirichlet on the wall-normal vorticity
!   phi  : Dirichlet by default, Neumann under -Dphi Neumann
!
! Half-channel and body-force variants are selected by the cpp switches, as
! before.
module channel_bcs

  use, intrinsic :: iso_c_binding
  use channel_grid
  use channel_state
  use compact_stencils

  implicit none

  public :: setup_boundary_conditions

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

end module channel_bcs
