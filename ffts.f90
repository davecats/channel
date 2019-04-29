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

MODULE ffts

  USE, intrinsic :: iso_c_binding
  USE typedef
  IMPLICIT NONE
  INCLUDE 'fftw3.f03'

  integer, save        :: plan_type=FFTW_PATIENT
  TYPE(C_PTR), save    :: pFFT,pIFT,pRFT,pHFT,ptrVVdx,ptrVVdz


CONTAINS

  SUBROUTINE init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB)
     integer(C_INT), intent(in) :: nxd,nxB,nzd,nzB
     complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:), intent(out) :: VVdx, VVdz
     real(C_DOUBLE), pointer, dimension(:,:,:,:), intent(out) :: rVVdx
     integer(C_INT), dimension(1) :: n_z, n_x, rn_x
     n_z=[nzd]; n_x=[nxd]; rn_x=[2*nxd];
     !Allocate aligned memory
     ptrVVdz=fftw_alloc_complex(int(nxB*nzd*6*6, C_SIZE_T))
     ptrVVdx=fftw_alloc_complex(int((nxd+1)*nzB*6*6, C_SIZE_T))
     !Convert C to F pointer
     CALL c_f_pointer(ptrVVdz, VVdz, [nzd,nxB,6,6]);  CALL c_f_pointer(ptrVVdx,   VVdx, [nxd+1,nzB,6,6])
                                                      CALL c_f_pointer(ptrVVdx,  rVVdx, [2*(nxd+1),nzB,6,6])
     !FFTs plans
     pFFT=fftw_plan_many_dft(1, n_z, nxB, VVdz(:,:,1,1), n_z, 1, nzd, VVdz(:,:,1,1), n_z, 1, nzd, FFTW_FORWARD,  plan_type)
     pIFT=fftw_plan_many_dft(1, n_z, nxB, VVdz(:,:,1,1), n_z, 1, nzd, VVdz(:,:,1,1), n_z, 1, nzd, FFTW_BACKWARD, plan_type)
     pRFT=fftw_plan_many_dft_c2r(1, rn_x, nzB,  VVdx(:,:,1,1),  n_x+1, 1,   (nxd+1), &
                                               rVVdx(:,:,1,1), rn_x+2, 1, 2*(nxd+1), plan_type)
     pHFT=fftw_plan_many_dft_r2c(1, rn_x, nzB, rVVdx(:,:,1,1), rn_x+2, 1, 2*(nxd+1), &
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
   END SUBROUTINE free_fft


END MODULE ffts
