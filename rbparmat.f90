!============================================!
!                                            !
!                  Parallel                  !  
!     Linear Algebra and LU Decomposition    !
!                  for the                   !
!     Direct Numerical Simulation (DNS)      !
!        of a turbulent channel flow         !
!                                            !
!============================================!
! 
! Author: Dr. Davide Gatti
! Date  : 16/Apr/2019
! 


!---- in-place LU Decomposition of a banded matrix ---!
!-----------------------------------------------------!
SUBROUTINE LU5decompStep(A)
  real(C_DOUBLE), intent(inout) :: A(0:,-2:)
  integer(C_INT) :: HI1,HI2
  real(C_DOUBLE) :: piv
  real(C_DOUBLE) :: buf(0:1,-2:2)
  TYPE(MPI_REQUEST) :: R
  TYPE(MPI_STATUS) :: S
  HI1=SIZE(A,1)-1; HI2=SIZE(A,2)-3;
  IF (last) THEN 
    A(HI1-2,1:2)=0; A(HI1-3,2)=0
  ELSE 
    CALL MPI_IRecv(buf,10,MPI_DOUBLE_PRECISION,ipy+1,TAG_LUDECOMP,MPI_COMM_Y,R)
    CALL MPI_Wait(R,S)
    A(HI1-1:HI1,:)=buf
  END IF
  DO i=HI1-HI2,0,-1 
    DO k=HI2,1,-1
      piv=A(i,k)
      DO j=-1,-2,-1
        A(i,j+k)=A(i,j+k)-piv*A(i+k,j)
      END DO
    END DO  
    piv=1.0d0/A(i,0); A(i,0)=piv; A(i,-2:-1)=A(i,-2:-1)*piv
  END DO
  IF (first) THEN 
    A(0,-2:-1)=0; A(1,-2)=0
  ELSE
    buf=A(0:1,:)
    CALL MPI_ISend(buf,10,MPI_DOUBLE_PRECISION,ipy-1,TAG_LUDECOMP,MPI_COMM_Y,R)
    CALL MPI_Wait(R,S)
  END IF
END SUBROUTINE LU5decompStep
!-----------------------------------------------------!


!------- Left LU division of a banded matrix ---------!
!------------------------step 1-----------------------!
SUBROUTINE LeftLU5divStep1(x,A,b)
  complex(C_DOUBLE_COMPLEX), intent(in) :: b(-2:)
  complex(C_DOUBLE_COMPLEX), intent(out) :: x(-2:)
  real(C_DOUBLE), intent(in)  :: A(0:,-2:)
  integer(C_INT) :: HI1,HI2
  TYPE(MPI_REQUEST) :: R
  TYPE(MPI_STATUS) :: S
  HI1=SIZE(A,1)-1; HI2=SIZE(A,2)-3;
  x=b
  IF (.NOT. last) THEN 
    CALL MPI_IRecv(x(HI1-1:HI1),2,cmpl,ipy+1,TAG_LUDIVSTEP1,MPI_COMM_Y,R)
    CALL MPI_Wait(R,S)
  END IF
  DO i=HI1-HI2,0,-1 
    x(i)=(x(i)-sum(A(i,1:2)*x(i+1:i+2)))*A(i,0)
  END DO
  IF (.NOT. first) THEN 
    CALL MPI_ISend(x(0:1),2,cmpl,ipy-1,TAG_LUDIVSTEP1,MPI_COMM_Y,R)
    CALL MPI_Wait(R,S)
  END IF
END SUBROUTINE LeftLU5divStep1
!-----------------------------------------------------!

!------- Left LU division of a banded matrix ---------!
!------------------------step 2-----------------------!
SUBROUTINE LeftLU5divStep2(A,b)
  complex(C_DOUBLE_COMPLEX), intent(inout) :: b(-2:)
  real(C_DOUBLE), intent(in)  :: A(0:,-2:)
  integer(C_INT) :: HI1,HI2
  TYPE(MPI_REQUEST) :: R
  TYPE(MPI_STATUS) :: S
  HI1=SIZE(A,1)-1; HI2=SIZE(A,2)-3;
  IF (.NOT. first) THEN 
    CALL MPI_IRecv(b(-2:-1),2,cmpl,ipy-1,TAG_LUDIVSTEP2,MPI_COMM_Y,R)
    CALL MPI_Wait(R,S)
  END IF
  DO i=0,HI1
    b(i)=b(i)-sum(A(i,-2:-1)*b(i-2:i-1))
  END DO
  IF (.NOT. last) THEN 
    CALL MPI_ISend(b(HI1-3:HI1-2),2,cmpl,ipy+1,TAG_LUDIVSTEP2,MPI_COMM_Y,R)
    CALL MPI_Wait(R,S)
  END IF
END SUBROUTINE LeftLU5divStep2
!-----------------------------------------------------!

