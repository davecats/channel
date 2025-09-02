!============================================!
!                                            !
!        Distributed Traspose (GPU)          !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: Dr. Davide Gatti
! GPU Adaptation: Gemini
!

#include "header.h"

MODULE mpi_transpose

  USE, intrinsic :: iso_c_binding
  USE, intrinsic :: iso_fortran_env
  USE mpi_f08

  IMPLICIT NONE

  TYPE(MPI_Comm) :: MPI_CART_COMM, MPI_COMM_X, MPI_COMM_Y
  integer(C_INT), save :: nproc,iproc,ierr,nzd,nx
  integer(C_INT), save :: nx0,nxN,nxB,nz0,nzN,nzB,ny0,nyN,miny,maxy,block
  complex(C_DOUBLE_COMPLEX), allocatable :: Ain(:),Aout(:)
  logical, save :: has_terminal,has_average
  TYPE(MPI_Datatype), save :: writeview_type,owned2write_type,vel_read_type,vel_field_type

CONTAINS

  !-------------- Transpose: Z to X --------------!
  !-----------------------------------------------!
  !$OMP DECLARE TARGET
  SUBROUTINE zTOx(Vz, Vx)
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vz(:,:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(:,:)
    integer(C_SIZE_T) :: i,j, ix, one, num, zero, npxt

    one = 1
    num = nzd/nproc
    zero = 0
    npxt = nproc

    ! Ain and Aout are host arrays for MPI communication.
    ! The packing/unpacking is done on the GPU for performance.

    ! Pack data on the GPU into the host-mapped Ain buffer
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    DO j=zero,npxt-one
      DO ix=nx0,nxN
        Ain( (j*nxB + (ix-nx0))*num + 1 : (j*nxB + (ix-nx0)+1)*num ) = Vz(one+j*num:(j+one)*num,ix-nx0+one)
      END DO
    END DO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

    ! Data is now packed in Ain on the host. Perform MPI transpose.
!$OMP END TARGET
    CALL MPI_Alltoall(Ain, nxB*nzB, MPI_DOUBLE_COMPLEX, Aout, nxB*nzB, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD)
!$OMP TARGET

    ! Unpack data on the GPU from the host-mapped Aout buffer
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    DO ix=0,nx
        Vx(ix+one,1:num) = Aout( ix*num + 1 : (ix+1)*num )
    END DO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

  END SUBROUTINE zTOx

  !-------------- Transpose: X to Z --------------!
  !-----------------------------------------------!
  !$OMP DECLARE TARGET
  SUBROUTINE xTOz(Vx,Vz)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(:,:)
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vx(:,:)
    integer(C_SIZE_T) :: i,j, iz, one, num, zero, npxt

    one = 1
    num = (nx+1)/nproc
    zero = 0
    npxt = nproc

    ! Pack data on the GPU into the host-mapped Ain buffer
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    DO j=zero,npxt-one
      DO iz=nz0,nzN
         Ain( (j*nzB + (iz-nz0))*num + 1 : (j*nzB + (iz-nz0)+1)*num ) = Vx(one+j*num:(j+one)*num,iz-nz0+one)
      END DO
    END DO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

!$OMP END TARGET
    ! Perform MPI transpose
    CALL MPI_Alltoall(Ain, nxB*nzB, MPI_DOUBLE_COMPLEX, Aout, nxB*nzB, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD)
!$OMP TARGET

    ! Unpack data on the GPU from the host-mapped Aout buffer
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
    DO iz=0,nzd-1
      Vz(iz+one,1:num) = Aout( iz*num+1 : (iz+1)*num )
    END DO
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

  END SUBROUTINE xTOz

  !------- Divide the problem in 1D slices -------!
  !-----------------------------------------------!
  SUBROUTINE init_MPI(nxpp,nz,ny,nxd_in,nzd_in)
    integer(C_INT), intent(in)  :: nxpp,nz,ny,nxd_in,nzd_in
    integer, parameter :: ndims = 4
    integer :: i
    integer :: array_of_sizes(ndims), array_of_subsizes(ndims), array_of_starts(ndims), ierror
    ! Define which process write on screen
    has_terminal=(iproc==0)
    ! Calculate domain division in wall-normal direction
    ny0=1; nyN=ny-1; miny=ny0-2; maxy=nyN+2
    ! Calculate domain division
    nx0=iproc*(nxpp)/nproc;  nxN=(iproc+1)*(nxpp)/nproc-1;  nxB=nxN-nx0+1;
    nz0=iproc*nzd_in/nproc;     nzN=(iproc+1)*nzd_in/nproc-1;     nzB=nzN-nz0+1;
    nx = nxpp - 1
    nzd = nzd_in
    block=max(nxB*nzd_in,nxpp*nzB)
    has_average=(nx0==0)
#ifdef mpiverbose
    DO i=0,nproc-1
       IF (iproc==i) WRITE(*,*) "iproc=",iproc,"nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, "ny0=", ny0, "nyN=", nyN
       CALL MPI_Barrier(MPI_COMM_WORLD)
    END DO
    FLUSH(output_unit)
#endif
    ! Allocate buffers for transposes. Use managed memory for easy GPU access.
    ALLOCATE(Ain(0:block-1)); Ain=0
    ALLOCATE(Aout(0:block-1)); Aout=0
!$OMP TARGET ENTER DATA MAP(TO: Ain, Aout)

    ! For READING VELOCITY, SETTING VIEW: datatype that maps velocity on disk to memory (it differs from writing: halo cells are read twice!)
    CALL MPI_Type_create_subarray(ndims, [ny+3, 2*nz+1, nxpp, 3], [nyN-ny0+5, 2*nz+1, nxB, 3], [ny0-1,0,nx0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_read_type, ierror)
    CALL MPI_Type_commit(vel_read_type, ierror)
    ! For READING VELOCITY, datatype in memory (avoids overflow) - NOTICE THAT THIS ARRAY IS FULL (NOT REALLY SUBARRAY)
    CALL MPI_Type_create_subarray(ndims, [nyN-ny0+5, 2*nz+1, nxB, 3], [nyN-ny0+5, 2*nz+1, nxB, 3], [0,0,0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_field_type, ierror)
    CALL MPI_Type_commit(vel_field_type, ierror)
    ! For WRITING VELOCITY, SETTING VIEW: datatype to map distributed velocity array to file
    array_of_sizes = [ny+3, 2*nz+1, nxpp, 3] ! size along each dimension of the WHOLE array ON DISK
    array_of_subsizes = [maxy-miny+1, 2*nz+1, nxB, 3] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny+1, 0, nx0, 0] ! starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, writeview_type, ierror)
    CALL MPI_Type_commit(writeview_type, ierror)
    ! For WRITING VELOCITY, SKIPPING HALO CELLS: datatype with holes to skip halo cells and select only data to be written
    array_of_sizes = [(nyN+2)-(ny0-2)+1, 2*nz+1, nxB, 3] ! size along each dimension of the array IN MEMORY owned by each process
    array_of_subsizes = [maxy-miny+1, 2*nz+1, nxB, 3] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny-(ny0-2), 0, 0, 0] ! starting position of each component; !!! IT'S ZERO BASED AND WRT TO ARRAY IN MEMORY !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, owned2write_type, ierror)
    CALL MPI_Type_commit(owned2write_type, ierror)
  END SUBROUTINE init_MPI


END MODULE mpi_transpose


