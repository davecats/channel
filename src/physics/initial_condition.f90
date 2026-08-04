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
  use channel_state, only: INITIAL_SEED_UNSET

  implicit none
  private

  public :: generate_initial_field, uniform_from_key

contains

  SUBROUTINE generate_initial_field(R, perturbation_amplitude, seed)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(INOUT) :: R(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3 + nPhi)
    real(C_DOUBLE), intent(in) :: perturbation_amplitude
    integer(C_INT), intent(in) :: seed
    integer(C_SIZE_T) :: ix, iy, iz, iPhi
    real(C_DOUBLE) :: rn(1:3)

    R = 0
    IF (has_terminal) WRITE (*, *) "Generating initial field..."
    if (seed == INITIAL_SEED_UNSET) then
      IF (has_terminal) WRITE (*, *) "  perturbation: unseeded, this field will not reproduce"
      DO iy = ny0 - 2, nyN + 2; DO ix = nx0, nxN; DO iz = -nz, nz
          CALL RANDOM_NUMBER(rn)
          R(iy, iz, ix, 1) = perturbation_amplitude*EXP(dcmplx(0, rn(1) - 0.5))
          R(iy, iz, ix, 2) = perturbation_amplitude*EXP(dcmplx(0, rn(2) - 0.5))
          R(iy, iz, ix, 3) = perturbation_amplitude*EXP(dcmplx(0, rn(3) - 0.5))
        END DO; END DO; END DO
    else
      IF (has_terminal) WRITE (*, '(A,I0)') "   perturbation: seed = ", seed
      DO iy = ny0 - 2, nyN + 2; DO ix = nx0, nxN; DO iz = -nz, nz
          rn(1) = uniform_from_key(seed, 1_C_INT, int(iy, C_INT), int(iz, C_INT), int(ix, C_INT))
          rn(2) = uniform_from_key(seed, 2_C_INT, int(iy, C_INT), int(iz, C_INT), int(ix, C_INT))
          rn(3) = uniform_from_key(seed, 3_C_INT, int(iy, C_INT), int(iz, C_INT), int(ix, C_INT))
          R(iy, iz, ix, 1) = perturbation_amplitude*EXP(dcmplx(0, rn(1) - 0.5))
          R(iy, iz, ix, 2) = perturbation_amplitude*EXP(dcmplx(0, rn(2) - 0.5))
          R(iy, iz, ix, 3) = perturbation_amplitude*EXP(dcmplx(0, rn(3) - 0.5))
        END DO; END DO; END DO
    end if
    IF (has_average) THEN
      DO iy = ny0 - 2, nyN + 2
        R(iy, 0, 0, 1) = 3*0.5*y(iy)*(2 - y(iy))
        DO iPhi = 1, nPhi
          R(iy, 0, 0, 3 + iPhi) = 3*0.5*y(iy)*(2 - y(iy))
        END DO
      END DO
    END IF
  END SUBROUTINE generate_initial_field

  ! A uniform deviate in [0, 1) determined by the seed and the *global* mode
  ! indices alone -- not by any running stream.
  !
  ! That is the whole point, and it is why this does not just call
  ! RANDOM_SEED.  A seeded intrinsic generator gives each rank a sequence, and
  ! which value lands on which mode then depends on how many ranks there are
  ! and which slab each one owns, so a seeded run at np=2 would not match the
  ! same seed at np=4.  Keying on (component, iy, iz, ix) instead makes the
  ! field a pure function of the deck, bit-for-bit at any decomposition.  It
  ! also makes each rank's ghost rows agree with the neighbour's interior rows,
  ! which the old loop did not.
  !
  ! Across *compilers* the agreement is a few ULP rather than bit-for-bit, and
  ! that is not this function's doing -- it returns identical bits under
  ! gfortran and nvfortran, measured.  The residue comes from the caller's
  ! EXP(dcmplx(0, phase)) and from the tanh mesh stretching, whose last bits are
  ! libm's.  Measured at 4 ULP over a small case.
  !
  ! The mixer is an fmix32-shaped avalanche, and every intermediate value is
  ! masked back to 32 bits so that NOTHING here can overflow.  That is not
  ! fastidiousness.  The obvious implementation is the splitmix64 finalizer,
  ! whose multiplications rely on 64-bit two's-complement wraparound -- which
  ! is signed integer overflow, and gfortran is entitled to assume it does not
  ! happen.  It duly miscompiles it: built -O0 that version agreed with a
  ! reference implementation to the bit, and built -O2 it returned the same
  ! constant for every key, i.e. a start field with no perturbation in it at
  ! all.  The CPU build happens to use -O0 today, so it would have worked until
  ! somebody raised the optimisation level.
  !
  ! Here the state is held in an integer(C_INT64_T) but is never allowed above
  ! 2^32, and the multipliers are odd and below 2^31, so every product is at
  ! most (2^32 - 1)*(2^31 - 1) < 2^63 and stays in range by construction.  A
  ! 32-bit deviate is ample: the value is used as a phase in [-0.5, 0.5].
  pure function uniform_from_key(seed, component, iy, iz, ix) result(u)
    IMPLICIT NONE
    integer(C_INT), intent(in) :: seed, component, iy, iz, ix
    real(C_DOUBLE) :: u
    integer(C_INT64_T) :: h
    integer(C_INT64_T), parameter :: MASK32 = 4294967295_C_INT64_T   ! 2^32 - 1
    integer(C_INT64_T), parameter :: MIX_A = 2146121005_C_INT64_T    ! odd, < 2^31
    integer(C_INT64_T), parameter :: MIX_B = 2032289749_C_INT64_T    ! odd, < 2^31

    ! Distinct odd strides per index, so no two modes share a key.  One operand
    ! of each product is an index of at most a few thousand, so these cannot
    ! overflow either.
    h = int(seed, C_INT64_T)
    h = h + 2654435761_C_INT64_T*int(component, C_INT64_T)
    h = h + 2246822519_C_INT64_T*int(iy, C_INT64_T)
    h = h + 3266489917_C_INT64_T*int(iz, C_INT64_T)
    h = h + 668265263_C_INT64_T*int(ix, C_INT64_T)
    h = iand(h, MASK32)

    h = ieor(h, ishft(h, -16))
    h = iand(h*MIX_A, MASK32)
    h = ieor(h, ishft(h, -13))
    h = iand(h*MIX_B, MASK32)
    h = ieor(h, ishft(h, -16))

    ! h is now a non-negative 32-bit value, which converts to double exactly.
    u = real(h, C_DOUBLE)*(1.0d0/4294967296.0d0)
  end function uniform_from_key

end module initial_condition
