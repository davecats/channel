!============================================!
!                                            !
!     Linear Algebra and LU Decomposition    !
!         (OpenMP GPU-accelerated)           !
!                  for the                   !
!     Direct Numerical Simulation (DNS)      !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: Dr. Davide Gatti
! GPU Adaptation: Gemini
!

MODULE rbmat

  USE, intrinsic :: iso_c_binding

  !- Overload division operator -!
  !------------------------------!
  INTERFACE operator (.bs.)
    PROCEDURE LLUdiv
  END INTERFACE OPERATOR (.bs.)

  INTERFACE operator (.bsr.)
    PROCEDURE LLU5div
  END INTERFACE OPERATOR (.bsr.)

  CONTAINS

  !- in-place LU Decomposition of a real square matrix -!
  !-----------------------Doolittle---------------------!
  !$OMP DECLARE TARGET
  SUBROUTINE LUdecompD(A)
    real(C_DOUBLE), intent(inout) :: A(:,:)
    integer(C_INT) :: HI, i, j
    real(C_DOUBLE) :: piv
    HI=SIZE(A,1)
    DO i=1,HI
     FORALL (j=i:HI)   A(i,j)=A(i,j)-sum(A(i,1:i-1)*A(1:i-1,j))
     piv=1/A(i,i)
     FORALL (j=i+1:HI) A(j,i)=(A(j,i)-sum(A(j,1:i-1)*A(1:i-1,i)))*piv
    END DO
  END SUBROUTINE LUdecompD

  !---- in-place LU Decomposition of a banded matrix ---!
  !-----------------------------------------------------!
  !$OMP DECLARE TARGET
  SUBROUTINE LU5decomp(A)
	  real(C_DOUBLE), intent(inout) :: A(:,-2:)
	  integer(C_INT) :: HI1,HI2, i,k,j
	  real(C_DOUBLE) :: piv
	  HI1=SIZE(A,1); HI2=2;
	  ! This algorithm has data dependencies, parallelizing inner loops is incorrect.
	  ! The parallelism comes from calling this routine for many independent
	  ! systems in linsolve.
	  DO i=HI1-HI2,1,-1
	    DO k=HI2,1,-1
	      piv=A(i,k)
	      DO j=-1,-2,-1
		      A(i,j+k)=A(i,j+k)-piv*A(i+k,j)
	      END DO
	    END DO
	    piv=1.0d0/A(i,0); A(i,0)=piv; A(i,-2:-1)=A(i,-2:-1)*piv
	  END DO
  END SUBROUTINE LU5decomp

  !- Left LU division of a banded matrix -!
  !------------- for complex vectors -----!
  !$OMP DECLARE TARGET
  SUBROUTINE LeftLU5div_c(x,A,b)
  complex(C_DOUBLE_COMPLEX), intent(out) :: x(:)
  complex(C_DOUBLE_COMPLEX), intent(in) :: b(:)
  real(C_DOUBLE), intent(in)  :: A(:,-2:)
  integer(C_INT) :: HI, i, j
  HI=SIZE(A,1)
  x=b
  ! Backward substitution (dependency)
  DO i=HI-2,1,-1
    x(i)=(x(i)-sum(A(i,1:2)*x(i+1:i+2)))*A(i,0)
  END DO
  ! Forward substitution (dependency)
  DO i=1,HI
    j = max(-2, 1-i)
    x(i)=x(i)-sum(A(i,j:-1)*x(i+j:i-1))
  END DO
  END SUBROUTINE LeftLU5div_c

  !- Left LU division of a banded matrix -!
  !--------------- for real vectors ------!
  !$OMP DECLARE TARGET
  SUBROUTINE LeftLU5div_s(x,A,b)
  real(C_DOUBLE), intent(out) :: x(:)
  real(C_DOUBLE), intent(in) :: b(:)
  real(C_DOUBLE), intent(in)  :: A(:,-2:)
  integer(C_INT) :: HI, i, j
  HI=SIZE(A,1)
  x=b
  DO i=HI-2,1,-1
    x(i)=(x(i)-sum(A(i,1:2)*x(i+1:i+2)))*A(i,0)
  END DO
  DO i=1,HI
    j = max(-2, 1-i)
    x(i)=x(i)-sum(A(i,j:-1)*x(i+j:i-1))
  END DO
  END SUBROUTINE LeftLU5div_s

  !---- OPERATOR Functions are not used in performance-critical paths, left on host ---!
  FUNCTION LLUdiv(A,b) RESULT(x)
    real(C_DOUBLE), intent(in) :: b(:)
    real(C_DOUBLE), intent(in) :: A(:,:)
    real(C_DOUBLE), allocatable :: x(:)
    integer(C_INT) :: HI, i
    HI=SIZE(A,1)
    ALLOCATE(x(1:HI))
    x(HI)=b(HI)*A(HI,HI)
    DO i=HI-1,1,-1
      x(i)=(b(i)-sum(A(i,i+1:HI)*x(i+1:HI)))*A(i,i)
    END DO
    DO i=2,HI
      x(i)=x(i)-sum(A(i,1:i-1)*x(1:i-1))
    END DO
  END FUNCTION LLUdiv

  FUNCTION LLU5div(A,b) RESULT(x)
    real(C_DOUBLE), intent(in) :: b(:)
    real(C_DOUBLE), intent(in) :: A(:,-2:)
    real(C_DOUBLE), allocatable :: x(:)
    integer(C_INT) :: HI,i,j
    HI=SIZE(A,1)
    ALLOCATE(x(1:HI))
    DO i=HI,1,-1
      j=MIN(2,HI-i)
      x(i)=(b(i)-sum(A(i,1:j)*x(i+1:i+j)))*A(i,0)
    END DO
    DO i=1,HI
      j=MAX(-2,1-i)
      x(i)=x(i)-sum(A(i,j:-1)*x(i+j:i-1))
    END DO
  END FUNCTION LLU5div


END MODULE rbmat


