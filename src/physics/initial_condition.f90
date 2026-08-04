#include "build_options.h"

! The field a run starts from when there is no restart file.
!
! This is a physics choice, and the only one a user is likely to want to
! change without touching the equations: which base profile the transient
! grows out of, and how it is disturbed.  It used to live in the ELSE branch of
! the restart *reader*, in src/io, behind three generations of commented-out
! alternatives -- so the one place a user would look for it was the one place
! it was not.
!
! generate_initial_field is called by restart_read when the file is absent.  It
! writes the spectral field directly: the mean mode (iz = 0, ix = 0) carries
! the base profile, and every mode gets a small random phase perturbation to
! seed the transition.
!
! CAUTION -- the base profile below is written for a channel of height 2.
! `3*0.5*y*(2 - y)` is the parabolic profile with unit bulk velocity over
! ymin = 0, ymax = 2, which is the convention every bundled deck uses.  It does
! NOT rescale itself to other ymin/ymax: with a different channel height it is
! still a parabola through y = 0 and y = 2, which is not the profile you asked
! for and may not even vanish at the walls.  Runs from a restart file are
! unaffected.  Rescaling this is a deliberate behaviour change -- every
! restartless run would start from a different field -- so it is left as it
! was found rather than fixed silently.
!
! Variants that were carried as commented-out lines at the old site, kept here
! because they are the examples a user starting from this file would want:
!
!   base profile with a wall-normal wave superposed
!     R(iy, 0, 0, 1) = 3*0.5*y(iy)*(2 - y(iy)) + 0.01*SIN(8*y(iy)*2*PI)/ni
!     R(iy, 0, 0, 1) = y(iy)*(2 - y(iy))*3.d0/2.d0 + 0.001*SIN(8*y(iy)*2*PI)
!   plane Couette base profile instead of the parabola
!     R(iy, 0, 0, 1) = y(iy) - 1
!   fixed perturbation amplitude, ignoring the deck value
!     R(iy, iz, ix, 1) = 0.0001*EXP(dcmplx(0, rn(1) - 0.5))
module initial_condition

  use, intrinsic :: iso_c_binding
  use channel_grid

  implicit none
  private

  public :: generate_initial_field

contains

  SUBROUTINE generate_initial_field(R, perturbation_amplitude)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(INOUT) :: R(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3 + nPhi)
    real(C_DOUBLE), intent(in) :: perturbation_amplitude
    integer(C_SIZE_T) :: ix, iy, iz, iPhi
    real(C_DOUBLE) :: rn(1:3)

    R = 0
    IF (has_terminal) WRITE (*, *) "Generating initial field..."
    DO iy = ny0 - 2, nyN + 2; DO ix = nx0, nxN; DO iz = -nz, nz
        CALL RANDOM_NUMBER(rn)
        R(iy, iz, ix, 1) = perturbation_amplitude*EXP(dcmplx(0, rn(1) - 0.5))
        R(iy, iz, ix, 2) = perturbation_amplitude*EXP(dcmplx(0, rn(2) - 0.5))
        R(iy, iz, ix, 3) = perturbation_amplitude*EXP(dcmplx(0, rn(3) - 0.5))
      END DO; END DO; END DO
    IF (has_average) THEN
      DO iy = ny0 - 2, nyN + 2
        R(iy, 0, 0, 1) = 3*0.5*y(iy)*(2 - y(iy))
        DO iPhi = 1, nPhi
          R(iy, 0, 0, 3 + iPhi) = 3*0.5*y(iy)*(2 - y(iy))
        END DO
      END DO
    END IF
  END SUBROUTINE generate_initial_field

end module initial_condition
