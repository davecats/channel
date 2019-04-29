!============================================!
!                                            !
!     Linear Algebra and LU Decomposition    !
!                  for the                   !
!     Direct Numerical Simulation (DNS)      !
!        of a turbulent channel flow         !
!                                            !
!============================================!
! 
! Author: Dr. Davide Gatti
! Date  : 28/Jul/2015
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
  SUBROUTINE LUdecompD(A)
    real(C_DOUBLE), intent(inout) :: A(:,:)
    integer(C_INT) :: HI
    real(C_DOUBLE) :: piv
    HI=SIZE(A,1)
    DO i=1,HI
     FORALL (j=i:HI)   A(i,j)=A(i,j)-sum(A(i,1:i-1)*A(1:i-1,j))
     piv=1/A(i,i)
     FORALL (j=i+1:HI) A(j,i)=(A(j,i)-sum(A(j,1:i-1)*A(1:i-1,i)))*piv
    END DO
  END SUBROUTINE LUdecompD

  !- in-place LU Decomposition of a real square matrix -!
  !-------------------------Crout-----------------------!
  SUBROUTINE LUdecompC(A)
    real(C_DOUBLE), intent(inout) :: A(:,:)
    integer(C_INT) :: HI
    real(C_DOUBLE) :: piv
    HI=SIZE(A,1)
    DO k=1,HI  
      FORALL (i=k:HI)   A(i,k)=A(i,k)-sum(A(i,1:k-1)*A(1:k-1,k))
      piv=1/A(k,k)
      FORALL (j=k+1:HI) A(k,j)=(A(k,j)-sum(A(k,1:k-1)*A(1:k-1,j)))*piv      
    END DO
  END SUBROUTINE LUdecompC

  !- in-place LU Decomposition of a real square matrix -!
  !-----------------------------------------------------!
  SUBROUTINE LUdecomp(A)
    real(C_DOUBLE), intent(inout) :: A(:,:)
    integer(C_INT) :: HI
    real(C_DOUBLE) :: piv
    HI=SIZE(A,1)
    DO i=HI,2,-1 
      piv=1.0d0/A(i,i); A(i,i)=piv; A(i,1:i-1)=A(i,1:i-1)*piv
      DO k=1,i-1
        piv=A(k,i)
        A(k,1:i-1)=A(k,1:i-1)-piv*A(i,1:i-1)
        !DO j=i-1,1,-1
        !  A(k,j)=A(k,j)-piv*A(i,j)
        !END DO
      END DO  
    END DO
    A(1,1)=1.0d0/A(1,1)
  END SUBROUTINE LUdecomp

  !---- in-place LU Decomposition of a banded matrix ---!
  !-----------------------------------------------------!
  SUBROUTINE LU5decomp(A)
    real(C_DOUBLE), intent(inout) :: A(0:,-2:)
    integer(C_INT) :: HI
    real(C_DOUBLE) :: piv
    HI=SIZE(A,1)-1;
    DO i=HI,0,-1 
      piv=1.0d0/A(i,0); A(i,0)=piv; A(i,-2:-1)=A(i,-2:-1)*piv
      DO k=max(-2,-i),-1 
        piv=A(i+k,-k)
        DO j=-1,-2,-1
        !A(i+k,M:2*M-1)=A(i+k,M:2*M-1)-piv*A(i,1:M)
          A(i+k,j-k)=A(i+k,j-k)-piv*A(i,j)
        END DO
      END DO  
    END DO
  END SUBROUTINE LU5decomp

  !- Left LU division of a square matrix -!
  !---------------Doolittle---------------!
  SUBROUTINE LeftLUdivD(x,A)
    real(C_DOUBLE), intent(inout) :: x(:)
    real(C_DOUBLE), intent(in) :: A(:,:)
    integer(C_INT) :: HI
    HI=SIZE(A,1)
    DO k=2,HI
      x(k)=x(k)-sum(A(k,1:k-1)*x(1:k-1))
    END DO
    x(HI)=x(HI)/A(HI,HI)
    DO k=HI-1,1,-1
      x(k)=(x(k)-sum(A(k,k+1:HI)*x(k+1:HI)))/A(k,k)
    END DO
  END SUBROUTINE LeftLUdivD

  !- Left LU division of a square matrix -!
  !-----------------Crout-----------------!
  SUBROUTINE LeftLUdivC(x,A)
    real(C_DOUBLE), intent(inout) :: x(:)
    real(C_DOUBLE), intent(in) :: A(:,:)
    integer(C_INT) :: HI
    HI=SIZE(A,1)
    x(1)=x(1)/A(1,1)
    DO i=2,HI
       x(i)=(x(i)-sum(A(i,1:i-1)*x(1:i-1)))/A(i,i)
    END DO 
    DO i=HI-1,1,-1
      x(i)=x(i)-sum(A(i,i+1:HI)*x(i+1:HI))
    END DO
  END SUBROUTINE LeftLUdivC

  !- Left LU division of a square matrix -!
  !---------------------------------------!
  SUBROUTINE LeftLUdiv(x,A)
    real(C_DOUBLE), intent(inout) :: x(:)
    real(C_DOUBLE), intent(in) :: A(:,:)
    integer(C_INT) :: HI
    HI=SIZE(A,1)
    x(HI)=x(HI)*A(HI,HI)
    DO i=HI-1,1,-1
      x(i)=(x(i)-sum(A(i,i+1:HI)*x(i+1:HI)))*A(i,i)
    END DO
    DO i=2,HI
      x(i)=x(i)-sum(A(i,1:i-1)*x(1:i-1))
    END DO
  END SUBROUTINE LeftLUdiv

  !- Left LU division of a banded matrix -!
  !---------------------------------------!
  SUBROUTINE LeftLU5div(A,b)
    complex(C_DOUBLE_COMPLEX), intent(inout)  :: b(:)
    real(C_DOUBLE), intent(in)  :: A(:,-2:)
    integer(C_INT) :: HI
    HI=SIZE(A,1)
    DO i=HI,1,-1
      j=MIN(2,HI-i)
      b(i)=(b(i)-sum(A(i,1:j)*b(i+1:i+j)))*A(i,0)
    END DO
    DO i=1,HI
      j=MAX(-2,1-i)
      b(i)=b(i)-sum(A(i,j:-1)*b(i+j:i-1))
    END DO
  END SUBROUTINE LeftLU5div

  !- Left LU division of a square matrix -!
  !----OPERATOR----Doolittle--------------!
  FUNCTION LLUdivD(A,b) RESULT(x)
    real(C_DOUBLE), intent(in) :: b(:)
    real(C_DOUBLE), intent(in) :: A(:,:)
    real(C_DOUBLE), allocatable :: x(:)
    integer(C_INT) :: HI
    HI=SIZE(A,1)
    ALLOCATE(x(1:HI))
    x(1)=b(1)
    DO k=2,HI
      x(k)=b(k)-sum(A(k,1:k-1)*x(1:k-1))
    END DO
    x(HI)=x(HI)/A(HI,HI)
    DO k=HI-1,1,-1
      x(k)=(x(k)-sum(A(k,k+1:HI)*x(k+1:HI)))/A(k,k)
    END DO
  END FUNCTION LLUdivD

  !- Left LU division of a square matrix -!
  !----OPERATOR-----Crout-----------------!
  FUNCTION LLUdivC(A,b) RESULT(x)
    real(C_DOUBLE), intent(in) :: b(:)
    real(C_DOUBLE), intent(in) :: A(:,:)
    real(C_DOUBLE), allocatable :: x(:)
    integer(C_INT) :: HI
    HI=SIZE(A,1)
    ALLOCATE(x(1:HI))
    x(1)=b(1)/A(1,1)
    DO i=2,HI
       x(i)=(b(i)-sum(A(i,1:i-1)*x(1:i-1)))/A(i,i)
    END DO 
    DO i=HI-1,1,-1
      x(i)=x(i)-sum(A(i,i+1:HI)*x(i+1:HI))
    END DO
  END FUNCTION LLUdivC

  !- Left LU division of a square matrix -!
  !----OPERATOR---------------------------!
  FUNCTION LLUdiv(A,b) RESULT(x)
    real(C_DOUBLE), intent(in) :: b(:)
    real(C_DOUBLE), intent(in) :: A(:,:)
    real(C_DOUBLE), allocatable :: x(:)
    integer(C_INT) :: HI
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

  !- Left LU division of a banded matrix -!
  !----OPERATOR---------------------------!
  FUNCTION LLU5div(A,b) RESULT(x)
    real(C_DOUBLE), intent(in) :: b(:)
    real(C_DOUBLE), intent(in) :: A(:,-2:)
    real(C_DOUBLE), allocatable :: x(:)
    integer(C_INT) :: HI
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
