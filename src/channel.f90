!============================================!
!                                            !
!     Direct Numerical Simulation (DNS)      !
!       of a turbulent channel flow          !
!                                            !
!============================================!
!
! This program has been written following the
! KISS (Keep it Simple and Stupid) philosophy
!
! Author: Dr.-Ing. Davide Gatti
!

#include "header.h"

PROGRAM channel
  USE driver
  IMPLICIT NONE

  CALL initialize("dns.in", "Dati.cart.out")
  CALL timeloop()
  CALL finalize()
END PROGRAM channel
