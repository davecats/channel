!============================================!
!                                            !
!           Fast Fourier Transforms          !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: Dr. Davide Gatti
! Date  : 28/Jul/2015
!

#include "header.h"

MODULE ffts

  USE, intrinsic :: iso_c_binding
  USE typedef
  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

#ifdef no_optimising_code
  integer, save        :: plan_type=FFTW_ESTIMATE
#else
  integer, save        :: plan_type=FFTW_PATIENT
#endif
  TYPE(C_PTR), save    :: pFFT,pIFT,pRFT,pHFT,ptrVVdx,ptrVVdz,ptrFdx,ptrFdz

CONTAINS

#ifdef ibm
   SUBROUTINE init_fft(VVdz,VVdx,rVVdx,Fdz,Fdx,rFdx,nxd,nxB,nzd,nzB,odd_n_real,s)
#else
   SUBROUTINE init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB,odd_n_real,s)
#endif
     integer(C_INT), intent(in) :: nxd,nxB,nzd,nzB
     complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:), intent(out) :: VVdx, VVdz
     real(C_DOUBLE), pointer, dimension(:,:,:,:), intent(out) :: rVVdx
     logical, optional, intent(in) :: odd_n_real
     integer, dimension(2), optional :: s
     integer, dimension(2) :: sn = 6
#ifdef ibm
     complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:), intent(out) :: Fdx, Fdz
     real(C_DOUBLE), pointer, dimension(:,:,:,:), intent(out) :: rFdx
#endif
     integer(C_INT), dimension(1) :: n_z, n_x, rn_x
     n_z=[nzd]; n_x=[nxd]; rn_x=[2*nxd];
     if (present(odd_n_real)) then
       ! notice: odd_n_real is basically .FALSE. by default
       ! meaning that by default the logical size of the real transform is even
       if (odd_n_real .eqv. .TRUE.) rn_x=rn_x - 1
     endif
     if (present(s)) sn = s
     !Allocate aligned memory
     ptrVVdz=fftw_alloc_complex(int(nxB*nzd*sn(1)*sn(2), C_SIZE_T))
     ptrVVdx=fftw_alloc_complex(int((nxd+1)*nzB*sn(1)*sn(2), C_SIZE_T))
#ifdef ibm
     ptrFdz=fftw_alloc_complex(int(nxB*nzd*3, C_SIZE_T))
     ptrFdx=fftw_alloc_complex(int((nxd+1)*nzB*3, C_SIZE_T))
#endif
     !Convert C to F pointer
     CALL c_f_pointer(ptrVVdz, VVdz, [nzd,nxB,sn(1),sn(2)]);  CALL c_f_pointer(ptrVVdx,   VVdx, [nxd+1,nzB,sn(1),sn(2)])
                                                      CALL c_f_pointer(ptrVVdx,  rVVdx, [2*(nxd+1),nzB,sn(1),sn(2)])
#ifdef ibm
     CALL c_f_pointer(ptrFdz,  Fdz,  [nzd,nxB,3,1]);  CALL c_f_pointer(ptrFdx,     Fdx, [nxd+1,nzB,3,1])
                                                      CALL c_f_pointer(ptrFdx,    rFdx, [2*(nxd+1),nzB,3,1])
#endif
     !FFTs plans
     pFFT=fftw_plan_many_dft(1, n_z, nxB, VVdz(:,:,1,1), n_z, 1, nzd, VVdz(:,:,1,1), n_z, 1, nzd, FFTW_FORWARD,  plan_type)
     pIFT=fftw_plan_many_dft(1, n_z, nxB, VVdz(:,:,1,1), n_z, 1, nzd, VVdz(:,:,1,1), n_z, 1, nzd, FFTW_BACKWARD, plan_type)
     pRFT=fftw_plan_many_dft_c2r(1, rn_x, nzB,  VVdx(:,:,1,1),  n_x+1, 1,   (nxd+1), &
                                               rVVdx(:,:,1,1), 2*(n_x+1), 1, 2*(nxd+1), plan_type)
     pHFT=fftw_plan_many_dft_r2c(1, rn_x, nzB, rVVdx(:,:,1,1), 2*(n_x+1), 1, 2*(nxd+1), &
                                                VVdx(:,:,1,1),  n_x+1, 1,   (nxd+1), plan_type)
   END SUBROUTINE init_fft

   LOGICAL FUNCTION fftFIT(i) result(isFIT)
     integer(C_INT), intent(in) :: i
     integer(C_INT) :: j
     j=i
     DO WHILE ( MOD(j,2)==0 )
       j=SHIFTA(j,1)
     END DO
     isFIT=( (j==1) .OR. (j==3) )
   END FUNCTION fftFIT

   SUBROUTINE FFT(x)
     complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:,:)
     CALL fftw_execute_dft(pFFT,x,x)
   END SUBROUTINE FFT

   SUBROUTINE IFT(x)
     complex(C_DOUBLE_COMPLEX), intent(inout) :: x(:,:)
     CALL fftw_execute_dft(pIFT,x,x)
   END SUBROUTINE IFT

   SUBROUTINE RFT(x,rx)
     complex(C_DOUBLE_COMPLEX) :: x(:,:)
     real(C_DOUBLE) :: rx(:,:)
     CALL fftw_execute_dft_c2r(pRFT,x,rx)
   END SUBROUTINE RFT

   SUBROUTINE HFT(rx,x)
     complex(C_DOUBLE_COMPLEX) :: x(:,:)
     real(C_DOUBLE) :: rx(:,:)
     CALL fftw_execute_dft_r2c(pHFT,rx,x)
   END SUBROUTINE HFT

   SUBROUTINE free_fft()
     CALL fftw_free(ptrVVdx); CALL fftw_free(ptrVVdz); 
#ifdef ibm
     CALL fftw_free(ptrFdx); CALL fftw_free(ptrFdz);
#endif
   END SUBROUTINE free_fft


END MODULE ffts
