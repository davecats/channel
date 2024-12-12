!============================================!
!                                            !
!              Type definitions              !
!                  for the                   !
!     Direct Numerical Simulation (DNS)      !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: M.Sc. Davide Gatti
! Date  : 11/Dec/2024
! 

#include "header.h"

MODULE typedef

  USE, intrinsic :: iso_c_binding
  IMPLICIT NONE

  TYPE :: RHSTYPE
    complex(C_DOUBLE_COMPLEX) :: eta,d2v
    complex(C_DOUBLE_COMPLEX) :: phi(1:nPhi)
  END TYPE RHSTYPE

  TYPE :: Di
    real(C_DOUBLE), dimension(-2:2) :: d0,d1,d2,d4
  END TYPE Di

  TYPE :: BCOND
    complex(C_DOUBLE_COMPLEX) :: u,v,w,vy,eta
    complex(C_DOUBLE_COMPLEX) :: phi(1:nPhi)
  END TYPE BCOND

END MODULE typedef
