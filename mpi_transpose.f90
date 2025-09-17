!============================================!
!                                            !
!            Distributed Traspose            !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: Dr. Davide Gatti
! Date  : 15/Apr/2019
!

#include "header.h"

MODULE mpi_transpose

  USE, intrinsic :: iso_c_binding
  USE, intrinsic :: iso_fortran_env
  USE mpi_f08

  IMPLICIT NONE

  TYPE(MPI_Comm) :: MPI_CART_COMM, MPI_COMM_X, MPI_COMM_Y
  integer(C_INT), save :: nproc, iproc, ierr, nzd, nx
  integer(C_INT), save :: nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN, miny, maxy
  !$omp declare target(ny0, nyN)
  complex(C_DOUBLE_COMPLEX), allocatable :: sendbuf(:), recvbuf(:)
  logical, save :: has_terminal, has_average
  TYPE(MPI_Datatype), save :: writeview_type, owned2write_type, vel_read_type, vel_field_type

CONTAINS

  !-------------- Transpose: Z to X --------------!
  !-----------------------------------------------!
  SUBROUTINE zTOx(Vz, Vx, ny, nPhi)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(inout)  :: Vz(1:, 1:, :, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(1:, 1:, :, :)
    integer(C_INT), intent(in) :: ny, nPhi
    integer(C_SIZE_T) :: iy, iV, ix, iz, dest, src
    integer :: nField
    integer :: sendcount, p

    IF (nproc == 1) THEN
      !$omp target teams distribute parallel do collapse(4) defaultmap(none) default(none) &
      !$omp shared(Vx, Vz) shared(ny, nzd, nx, nPhi) private(iy, iV, ix, iz)
      DO iy = 1, ny + 3
        DO iV = 1, 3 + nPhi
          DO ix = 1, nx + 1
            DO iz = 1, nzd
              Vx(ix, iz, iy, iV) = Vz(iz, ix, iy, iV)
            END DO
          END DO
        END DO
      END DO

    ELSE

      nField = 3 + nPhi        ! same as your iV loop upper bound
      sendcount = nxB*nzB*(ny + 3)*nField

      !$omp target teams distribute parallel do collapse(5) defaultmap(none) default(none) &
      !$omp shared(Vz, sendbuf) shared(ny, nxB, nzB, nproc, nField, sendcount) private(iy, iV, ix, iz, dest, p)
      do dest = 0, nproc - 1
        do iV = 1, nField
          do iy = 1, ny + 3
            do ix = 1, nxB
              do iz = 1, nzB
                p = dest*sendcount + iz + (nzB*(ix - 1)) + (nzB*nxB*(iy - 1)) + (nzB*nxB*(ny + 3)*(iV - 1))
                sendbuf(p) = Vz(dest*nzB + iz, ix, iy, iV)
              end do
            end do
          end do
        end do
      end do

      !$omp target data use_device_ptr(sendbuf, recvbuf)
      call MPI_Alltoall(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, &
                        recvbuf, sendcount, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

      !$omp end target data
      !$omp target teams distribute parallel do collapse(5) defaultmap(none) default(none) &
      !$omp shared(Vx, recvbuf) shared(ny, nxB, nzB, nproc, nField, sendcount) private(iy, iV, ix, iz, dest, p)
      do src = 0, nproc - 1
        do iV = 1, nField
          do iy = 1, ny + 3
            do ix = 1, nxB
              do iz = 1, nzB
                p = src*sendcount + iz + (nzB*(ix - 1)) + (nzB*nxB*(iy - 1)) + (nzB*nxB*(ny + 3)*(iV - 1))
                Vx(ix + src*nxB, iz, iy, iV) = recvbuf(p)
              end do
            end do
          end do
        end do
      end do

    END IF

  END SUBROUTINE zTOx

  !-------------- Transpose: X to Z --------------!
  !-----------------------------------------------!
  SUBROUTINE xTOz(Vx, Vz, ny, nPhi)
    use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
    implicit none

    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(1:, 1:, :, :)
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vx(1:, 1:, :, :)
    integer(C_INT), intent(in) :: ny, nPhi
    integer(C_SIZE_T) :: iy, iV, ix, iz, dest, src
    integer :: nField
    integer :: sendcount, p

    ! --- single-rank short-circuit ---
    IF (nproc == 1) THEN
      !$omp target teams distribute parallel do collapse(4) defaultmap(none) default(none) &
      !$omp shared(Vx, Vz) shared(ny, nzd, nx, nPhi) private(iy, iV, ix, iz)
      DO iy = 1, ny + 3
        DO iV = 1, 6 + 3*nPhi
          DO iz = 1, nzd
            DO ix = 1, nx + 1
              Vz(iz, ix, iy, iV) = Vx(ix, iz, iy, iV)
            END DO
          END DO
        END DO
      END DO

    ELSE
      nField = 6 + 3*nPhi
      sendcount = nxB*nzB*(ny + 3)*nField

      !$omp target teams distribute parallel do collapse(5) defaultmap(none) default(none) &
      !$omp shared(Vx, sendbuf) shared(ny, nxB, nzB, nproc, nField, sendcount) private(iy, iV, ix, iz, dest, p)
      do dest = 0, nproc - 1
        do iV = 1, nField
          do iy = 1, ny + 3
            do iz = 1, nzB
              do ix = 1, nxB
                p = dest*sendcount + ix + (nxB*(iz - 1)) + (nxB*nzB*(iy - 1)) + (nxB*nzB*(ny + 3)*(iV - 1))
                sendbuf(p) = Vx(dest*nxB + ix, iz, iy, iV)
              end do
            end do
          end do
        end do
      end do

      !$omp target data use_device_ptr(sendbuf, recvbuf)
      call MPI_Alltoall(sendbuf, sendcount, MPI_DOUBLE_COMPLEX, &
                        recvbuf, sendcount, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)
      !$omp end target data

      !$omp target teams distribute parallel do collapse(5) defaultmap(none) default(none) &
      !$omp shared(Vz, recvbuf) shared(ny, nxB, nzB, nproc, nField, sendcount) private(iy, iV, ix, iz, dest, p)
      do src = 0, nproc - 1
        do iV = 1, nField
          do iy = 1, ny + 3
            do iz = 1, nzB
              do ix = 1, nxB
                p = src*sendcount + ix + (nxB*(iz - 1)) + (nxB*nzB*(iy - 1)) + (nxB*nzB*(ny + 3)*(iV - 1))
                Vz(iz + src*nzB, ix, iy, iV) = recvbuf(p)
              end do
            end do
          end do
        end do
      end do
    END IF

  END SUBROUTINE xTOz

  !------- Divide the problem in 1D slices -------!
  !-----------------------------------------------!
  SUBROUTINE init_MPI(nxpp, nz, ny, nxd, nzd, nPhi)
    integer(C_INT), intent(in)  :: nxpp, nz, ny, nxd, nzd, nPhi
    integer, parameter :: ndims = 4
    integer :: i
    integer :: array_of_sizes(ndims), array_of_subsizes(ndims), array_of_starts(ndims), ierror
    ! Define which process write on screen
    has_terminal = (iproc == 0)
    ! Calculate domain division in wall-normal direction
    ny0 = 1; nyN = ny - 1; miny = ny0 - 2; maxy = nyN + 2
    !$omp target update to(ny0, nyN)

    ! Calculate domain division
    nx0 = iproc*(nxpp)/nproc; nxN = (iproc + 1)*(nxpp)/nproc - 1; nxB = nxN - nx0 + 1; 
    nz0 = iproc*nzd/nproc; nzN = (iproc + 1)*nzd/nproc - 1; nzB = nzN - nz0 + 1; 
    has_average = (nx0 == 0)
#ifdef mpiverbose
    DO i = 0, nproc - 1
       IF (iproc==i) WRITE(*,*) "iproc=",iproc,"nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, "ny0=", ny0, "nyN=", nyN 
      CALL MPI_Barrier(MPI_COMM_WORLD)
    END DO
    FLUSH (output_unit)
#endif
    if (int(nproc, 8)*int(nxB, 8)*int(nzB, 8)*int(ny + 3, 8)*int(6 + 3*nPhi, 8) > huge(0_C_INT)) then
      if (has_terminal) then
        print *, "Error: problem too large for MPI transpose (integer overflow). Try to increase the number of processes."
      end if
      CALL MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    end if

    ! Allocate buffers for transposes
    ALLOCATE (sendbuf(nproc*nxB*nzB*(ny + 3)*(6 + 3*nPhi))); sendbuf = 0
    ALLOCATE (recvbuf(nproc*nxB*nzB*(ny + 3)*(6 + 3*nPhi))); recvbuf = 0
    !$omp target enter data map(alloc: sendbuf, recvbuf)
    ! For READING VELOCITY, SETTING VIEW: datatype that maps velocity on disk to memory (it differs from writing: halo cells are read twice!)
    CALL MPI_Type_create_subarray(ndims, [ny+3, 2*nz+1, nxpp, 3+nPhi], [nyN-ny0+5, 2*nz+1, nxB, 3+nPhi], [ny0-1,0,nx0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_read_type, ierror)
    CALL MPI_Type_commit(vel_read_type, ierror)
    ! For READING VELOCITY, datatype in memory (avoids overflow) - NOTICE THAT THIS ARRAY IS FULL (NOT REALLY SUBARRAY)
    CALL MPI_Type_create_subarray(ndims, [nyN-ny0+5, 2*nz+1, nxB, 3+nPhi], [nyN-ny0+5, 2*nz+1, nxB, 3+nPhi], [0,0,0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_field_type, ierror)
    CALL MPI_Type_commit(vel_field_type, ierror)
    ! For WRITING VELOCITY, SETTING VIEW: datatype to map distributed velocity array to file
    array_of_sizes = [ny + 3, 2*nz + 1, nxpp, 3 + nPhi] ! size along each dimension of the WHOLE array ON DISK
    array_of_subsizes = [maxy - miny + 1, 2*nz + 1, nxB, 3 + nPhi] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny + 1, 0, nx0, 0] ! starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, writeview_type, ierror)
    CALL MPI_Type_commit(writeview_type, ierror)
    ! For WRITING VELOCITY, SKIPPING HALO CELLS: datatype with holes to skip halo cells and select only data to be written
    array_of_sizes = [(nyN + 2) - (ny0 - 2) + 1, 2*nz + 1, nxB, 3 + nPhi] ! size along each dimension of the array IN MEMORY owned by each process
    array_of_subsizes = [maxy - miny + 1, 2*nz + 1, nxB, 3] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny - (ny0 - 2), 0, 0, 0] ! starting position of each component; !!! IT'S ZERO BASED AND WRT TO ARRAY IN MEMORY !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, owned2write_type, ierror)
    CALL MPI_Type_commit(owned2write_type, ierror)
  END SUBROUTINE init_MPI

END MODULE mpi_transpose
