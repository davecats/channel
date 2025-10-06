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

#ifdef HAVE_HIP
module roctx

  implicit none

  private
  integer :: n = 0

  public :: roctxpush, roctxpop

  interface
    subroutine roctxrangepush(message) bind(c, name="roctxRangePushA")
      use iso_c_binding, only: c_char
      implicit none
      character(c_char) :: message(*)
    end subroutine roctxrangepush

    subroutine roctxrangepop() bind(c, name="roctxRangePop")
      implicit none
    end subroutine roctxrangepop

  end interface

contains

  subroutine roctxPush(name)
    character(len=*), intent(in) :: name
    n = n + 1
    call roctxRangePush(name)
  end subroutine roctxPush

  subroutine roctxPop(name)
    character(len=*), intent(in) :: name
    n = n - 1
    ! Print the marker name if there are more pop calls than push calls
    if (n < 0) then
      print *, "invalid pop for: ", name
      return
    end if
    call roctxRangePop()
  end subroutine roctxPop

end module
#endif

MODULE mpi_transpose

  USE, intrinsic :: iso_c_binding
  USE, intrinsic :: iso_fortran_env
#ifdef HAVE_MPI
  USE mpi_f08
#endif
#if defined(HAVE_HIP)
  use omp_lib
  use roctx
#endif

  IMPLICIT NONE

#ifdef HAVE_MPI
  TYPE(MPI_Comm) :: MPI_CART_COMM, MPI_COMM_X, MPI_COMM_Y
#endif
#if defined(HAVE_HIP)
  complex(C_DOUBLE_COMPLEX), pointer:: sendbuf(:, :), recvbuf(:, :)
#else
  complex(C_DOUBLE_COMPLEX), allocatable :: sendbuf(:, :), recvbuf(:, :)
#endif
  integer(C_INT), save :: nproc, iproc, ierr, nzd, nx
  integer(C_INT), save :: nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN, miny, maxy, sendcount
  !$omp declare target(ny0, nyN)

  logical, save :: has_terminal, has_average
#ifdef HAVE_MPI
  TYPE(MPI_Datatype), save :: writeview_type, owned2write_type, vel_read_type, vel_field_type
#endif

CONTAINS

  SUBROUTINE pack_zTOx(Vz, send, ny)
    use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vz(1:, 1:, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: send(:)
    integer(C_INT), intent(in)  :: ny
    integer(C_SIZE_T) :: iy, ix, iz, dest, p

    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vz, send) shared(ny, nxB, nzB, nproc, sendcount) private(iy, ix, iz, dest, p)
    do dest = 0, nproc - 1
      do iy = 1, ny + 3
        do ix = 1, nxB
          do iz = 1, nzB
            p = dest*sendcount + iz + (nzB*(ix - 1)) + (nzB*nxB*(iy - 1))
            send(p) = Vz(dest*nzB + iz, ix, iy)
          end do
        end do
      end do
    end do

  END SUBROUTINE pack_zTOx

  SUBROUTINE unpack_zTOx(recv, Vx, ny)
    use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in)  :: recv(:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(1:, 1:, :)
    integer(C_INT), intent(in)  :: ny
    integer(C_SIZE_T) :: iy, ix, iz, src, p

    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vx, recv) shared(ny, nxB, nzB, nproc, sendcount) private(iy, ix, iz, src, p)
    do src = 0, nproc - 1
      do iy = 1, ny + 3
        do ix = 1, nxB
          do iz = 1, nzB
            p = src*sendcount + iz + (nzB*(ix - 1)) + (nzB*nxB*(iy - 1))
            Vx(ix + src*nxB, iz, iy) = recv(p)
          end do
        end do
      end do
    end do
  END SUBROUTINE unpack_zTOx

  SUBROUTINE pack_xTOz(Vx, send, ny)
    use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vx(1:, 1:, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: send(:)
    integer(C_INT), intent(in)  :: ny
    integer(C_SIZE_T) :: iy, ix, iz, dest, p

    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vx, send) shared(ny, nxB, nzB, nproc, sendcount) private(iy, ix, iz, dest, p)
    do dest = 0, nproc - 1
      do iy = 1, ny + 3
        do iz = 1, nzB
          do ix = 1, nxB
            p = dest*sendcount + ix + (nxB*(iz - 1)) + (nxB*nzB*(iy - 1))
            send(p) = Vx(dest*nxB + ix, iz, iy)
          end do
        end do
      end do
    end do
  END SUBROUTINE pack_xTOz

  SUBROUTINE unpack_xTOz(recv, Vz, ny)
    use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in)  :: recv(:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(1:, 1:, :)
    integer(C_INT), intent(in)  :: ny
    integer(C_SIZE_T) :: iy, ix, iz, src, p

    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vz, recv) shared(ny, nxB, nzB, nproc, sendcount) private(iy, ix, iz, src, p)
    do src = 0, nproc - 1
      do iy = 1, ny + 3
        do iz = 1, nzB
          do ix = 1, nxB
            p = src*sendcount + ix + (nxB*(iz - 1)) + (nxB*nzB*(iy - 1))
            Vz(iz + src*nzB, ix, iy) = recv(p)
          end do
        end do
      end do
    end do
  END SUBROUTINE unpack_xTOz

  SUBROUTINE alltoall(send, recv, request)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(out) :: recv(:)
    complex(C_DOUBLE_COMPLEX), intent(in)  :: send(:)
    type(MPI_Request), intent(inout) :: request
#ifndef HAVE_HIP
    !$omp target data use_device_ptr(send, recv)
#endif
    call MPI_IALLTOALL(send, sendcount, MPI_DOUBLE_COMPLEX, &
                       recv, sendcount, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, request, ierr)
#ifndef HAVE_HIP
    !$omp end target data
#endif

  END SUBROUTINE alltoall

  !------- Divide the problem in 1D slices -------!
  !-----------------------------------------------!
  SUBROUTINE init_MPI(nxpp, nz, ny, nxd, nzd, nPhi)
    integer(C_INT), intent(in)  :: nxpp, nz, ny, nxd, nzd, nPhi
    integer, parameter :: ndims = 4
    integer :: i
    integer :: array_of_sizes(ndims), array_of_subsizes(ndims), array_of_starts(ndims), ierror
    type(c_ptr) :: sendptr, recvptr
    integer(c_size_t) :: sendsize, recvsize
    ! Define which process write on screen
    has_terminal = (iproc == 0)
    ! Calculate domain division in wall-normal direction
    ny0 = 1; nyN = ny - 1; miny = ny0 - 2; maxy = nyN + 2
    !$omp target update to(ny0, nyN)

    ! Calculate domain division
    nx0 = iproc*(nxpp)/nproc; nxN = (iproc + 1)*(nxpp)/nproc - 1; nxB = nxN - nx0 + 1; 
    nz0 = iproc*nzd/nproc; nzN = (iproc + 1)*nzd/nproc - 1; nzB = nzN - nz0 + 1; 
    has_average = (nx0 == 0)
#ifdef HAVE_MPI
#ifdef mpiverbose
    DO i = 0, nproc - 1
       IF (iproc==i) WRITE(*,*) "iproc=",iproc,"nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, "ny0=", ny0, "nyN=", nyN 
      CALL MPI_Barrier(MPI_COMM_WORLD)
    END DO
    FLUSH (output_unit)
#endif
    if (int(nproc, 8)*int(nxB, 8)*int(nzB, 8)*int(ny + 3, 8) > huge(0_C_INT)) then
      if (has_terminal) then
        print *, "Error: problem too large for MPI transpose (integer overflow). Try to increase the number of processes."
      end if
      CALL MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
    sendsize = nproc*nxB*nzB*(ny + 3)
    recvsize = nproc*nxB*nzB*(ny + 3)

    sendcount = nxB*nzB*(ny + 3)

    ! Allocate buffers for transposes*int(16, c_size_t)
#if defined(HAVE_HIP)
    ! On HIP with HSA_XNACK=1, the use_device_ptr statements around the MPI calls are ignored.
    ! Hence, MPI does a CPU mpi copy! So we need to allocate it explicity on the device.
    sendptr = omp_target_alloc(sendsize*int(16*(2), c_size_t), omp_get_default_device())
    recvptr = omp_target_alloc(recvsize*int(16*(2), c_size_t), omp_get_default_device())
    call c_f_pointer(sendptr, sendbuf, [sendsize, int(2, c_size_t)])
    call c_f_pointer(recvptr, recvbuf, [recvsize, int(2, c_size_t)])
#else
    ALLOCATE (sendbuf(sendsize, 2)); sendbuf = 0
    ALLOCATE (recvbuf(recvsize, 2)); recvbuf = 0
    !$omp target enter data map(alloc: sendbuf, recvbuf)
#endif

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
#endif
  END SUBROUTINE init_MPI

END MODULE mpi_transpose
