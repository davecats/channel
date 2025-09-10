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
  integer(C_INT), save :: nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN, miny, maxy, block
  !$omp declare target(ny0, nyN)
  complex(C_DOUBLE_COMPLEX), allocatable :: Ain(:), Aout(:)
  logical, save :: has_terminal, has_average
  TYPE(MPI_Datatype), save :: writeview_type, owned2write_type, vel_read_type, vel_field_type

CONTAINS

  !-------------- Transpose: Z to X --------------!
  !-----------------------------------------------!
  SUBROUTINE zTOx(Vz, Vx)
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vz(1:, 1:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(1:, 1:)
    integer(C_SIZE_T) :: i, j
    integer(C_SIZE_T) :: ix, iz, one, num, zero, npxt

    one = 1
    num = nzd/nproc
    zero = 0
    npxt = nproc

    Ain = 0; Aout = 0

    !Pack
    i = zero
    DO j = zero, npxt - one
      DO ix = nx0, nxN
        Ain(i:i + num - one) = Vz(one + j*num:(j + one)*num, ix - nx0 + one)
        i = i + num
      END DO
    END DO

    !Send
    CALL MPI_Alltoall(Ain, nxB*nzB, MPI_DOUBLE_COMPLEX, Aout, nxB*nzB, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD)

    !Unpack
    i = zero
    DO ix = 0, nx
      Vx(ix + one, one:num) = Aout(i:i + num - one)
      i = i + num
    END DO
  END SUBROUTINE zTOx

  !-------------- Transpose: X to Z --------------!
  !-----------------------------------------------!
  SUBROUTINE xTOz(Vx, Vz)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(1:, 1:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vx(1:, 1:)
    integer(C_SIZE_T) :: i, j
    integer(C_SIZE_T) :: ix, iz, one, num, zero, npxt

    one = 1
    num = (nx + 1)/nproc
    zero = 0
    npxt = nproc

    !Pack
    i = zero
    DO j = zero, npxt - one
      DO iz = nz0, nzN
        Ain(i:i + num - one) = Vx(one + j*num:(j + one)*num, iz - nz0 + one)
        i = i + num
      END DO
    END DO

    !Send
    CALL MPI_Alltoall(Ain, nxB*nzB, MPI_DOUBLE_COMPLEX, Aout, nxB*nzB, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD)

    !Unpack
    i = zero
    DO iz = 0, nzd - 1
      Vz(iz + one, one:num) = Aout(i:i + num - one)
      i = i + num
    END DO
  END SUBROUTINE xTOz

  !------- Divide the problem in 1D slices -------!
  !-----------------------------------------------!
  SUBROUTINE init_MPI(nxpp, nz, ny, nxd, nzd)
    integer(C_INT), intent(in)  :: nxpp, nz, ny, nxd, nzd
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
    block = max(nxB*nzd, nxpp*nzB)
    has_average = (nx0 == 0)
#ifdef mpiverbose
    DO i = 0, nproc - 1
       IF (iproc==i) WRITE(*,*) "iproc=",iproc,"nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, "ny0=", ny0, "nyN=", nyN 
      CALL MPI_Barrier(MPI_COMM_WORLD)
    END DO
    FLUSH (output_unit)
#endif
    ! Allocate buffers for transposes
    ALLOCATE (Ain(0:block - 1)); Ain = 0
    ALLOCATE (Aout(0:block - 1)); Aout = 0
    ! For READING VELOCITY, SETTING VIEW: datatype that maps velocity on disk to memory (it differs from writing: halo cells are read twice!)
    CALL MPI_Type_create_subarray(ndims, [ny+3, 2*nz+1, nxpp, 3], [nyN-ny0+5, 2*nz+1, nxB, 3], [ny0-1,0,nx0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_read_type, ierror)
    CALL MPI_Type_commit(vel_read_type, ierror)
    ! For READING VELOCITY, datatype in memory (avoids overflow) - NOTICE THAT THIS ARRAY IS FULL (NOT REALLY SUBARRAY)
    CALL MPI_Type_create_subarray(ndims, [nyN-ny0+5, 2*nz+1, nxB, 3], [nyN-ny0+5, 2*nz+1, nxB, 3], [0,0,0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_field_type, ierror)
    CALL MPI_Type_commit(vel_field_type, ierror)
    ! For WRITING VELOCITY, SETTING VIEW: datatype to map distributed velocity array to file
    array_of_sizes = [ny + 3, 2*nz + 1, nxpp, 3] ! size along each dimension of the WHOLE array ON DISK
    array_of_subsizes = [maxy - miny + 1, 2*nz + 1, nxB, 3] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny + 1, 0, nx0, 0] ! starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, writeview_type, ierror)
    CALL MPI_Type_commit(writeview_type, ierror)
    ! For WRITING VELOCITY, SKIPPING HALO CELLS: datatype with holes to skip halo cells and select only data to be written
    array_of_sizes = [(nyN + 2) - (ny0 - 2) + 1, 2*nz + 1, nxB, 3] ! size along each dimension of the array IN MEMORY owned by each process
    array_of_subsizes = [maxy - miny + 1, 2*nz + 1, nxB, 3] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny - (ny0 - 2), 0, 0, 0] ! starting position of each component; !!! IT'S ZERO BASED AND WRT TO ARRAY IN MEMORY !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, owned2write_type, ierror)
    CALL MPI_Type_commit(owned2write_type, ierror)
  END SUBROUTINE init_MPI

END MODULE mpi_transpose
