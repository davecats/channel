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
#ifdef HAVE_MPI
  USE mpi_f08
#endif
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  use omp_lib
#endif
  use roctx

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
  integer(C_INT), save :: npy_grid = 1, npxz = 1, ipy = 0, ipxz = 0
  integer(C_INT), save :: nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN, sendcount
  !$omp declare target(npy_grid, npxz, ipy, ipxz, nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN, sendcount)

  logical, save :: has_terminal, has_average, fft_transpose_is_local
#ifdef HAVE_MPI
  TYPE(MPI_Datatype), save :: writeview_type, owned2write_type, vel_read_type, vel_field_type
#endif

CONTAINS

  SUBROUTINE split_block(total, nparts, part, start, count)
    integer(C_INT), intent(in) :: total, nparts, part
    integer(C_INT), intent(out) :: start, count
    integer(C_INT) :: base, rem

    base = total/nparts
    rem = mod(total, nparts)
    count = base
    if (part < rem) count = count + 1
    start = part*base + min(part, rem) + 1
  END SUBROUTINE split_block

  SUBROUTINE repack_zTOx_local(Vz, Vx, ny)
    use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vz(1:, 1:, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(1:, 1:, :)
    integer(C_INT), intent(in) :: ny
    integer(C_SIZE_T) :: iy, ix, iz
    integer(C_INT) :: ny_batch

    ny_batch = size(Vz, 3)
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(Vz, Vx) shared(ny_batch, nxB, nzd) private(iy, ix, iz)
    do iy = 1, ny_batch
      do ix = 1, nxB
        do iz = 1, nzd
          Vx(ix, iz, iy) = Vz(iz, ix, iy)
        end do
      end do
    end do
  END SUBROUTINE repack_zTOx_local

  SUBROUTINE repack_xTOz_local(Vx, Vz, ny)
    use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vx(1:, 1:, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(1:, 1:, :)
    integer(C_INT), intent(in) :: ny
    integer(C_SIZE_T) :: iy, ix, iz
    integer(C_INT) :: ny_batch

    ny_batch = size(Vx, 3)
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(Vx, Vz) shared(ny_batch, nxB, nzd) private(iy, ix, iz)
    do iy = 1, ny_batch
      do iz = 1, nzd
        do ix = 1, nxB
          Vz(iz, ix, iy) = Vx(ix, iz, iy)
        end do
      end do
    end do
  END SUBROUTINE repack_xTOz_local

  SUBROUTINE pack_zTOx(Vz, send, ny)
    use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vz(1:, 1:, :)
    complex(C_DOUBLE_COMPLEX), intent(out) :: send(:)
    integer(C_INT), intent(in)  :: ny
    integer(C_SIZE_T) :: iy, ix, iz, dest, p
    integer(C_INT) :: ny_batch

    ny_batch = size(Vz, 3)
    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vz, send) shared(ny_batch, nxB, nzB, npxz, sendcount) private(iy, ix, iz, dest, p)
    do dest = 0, npxz - 1
      do iy = 1, ny_batch
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
    integer(C_INT) :: ny_batch

    ny_batch = size(Vx, 3)
    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vx, recv) shared(ny_batch, nxB, nzB, npxz, sendcount) private(iy, ix, iz, src, p)
    do src = 0, npxz - 1
      do iy = 1, ny_batch
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
    integer(C_INT) :: ny_batch

    ny_batch = size(Vx, 3)
    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vx, send) shared(ny_batch, nxB, nzB, npxz, sendcount) private(iy, ix, iz, dest, p)
    do dest = 0, npxz - 1
      do iy = 1, ny_batch
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
    integer(C_INT) :: ny_batch

    ny_batch = size(Vz, 3)
    !$omp target teams distribute parallel do collapse(4) default(none) &
    !$omp shared(Vz, recv) shared(ny_batch, nxB, nzB, npxz, sendcount) private(iy, ix, iz, src, p)
    do src = 0, npxz - 1
      do iy = 1, ny_batch
        do iz = 1, nzB
          do ix = 1, nxB
            p = src*sendcount + ix + (nxB*(iz - 1)) + (nxB*nzB*(iy - 1))
            Vz(iz + src*nzB, ix, iy) = recv(p)
          end do
        end do
      end do
    end do
  END SUBROUTINE unpack_xTOz

  SUBROUTINE alltoall(send, recv, request, label)
    IMPLICIT NONE
    complex(C_DOUBLE_COMPLEX), intent(out) :: recv(:)
    complex(C_DOUBLE_COMPLEX), intent(in)  :: send(:)
    type(MPI_Request), intent(inout) :: request
    character(len=*), intent(in), optional :: label
    character(len=96) :: range_name

    range_name = "MPI_Ialltoall fft_transpose"
    if (present(label)) range_name = "MPI_Ialltoall "//trim(label)
    call roctxPush(range_name)
#ifndef HAVE_HIP
    !$omp target data use_device_ptr(send, recv)
#endif
    call MPI_IALLTOALL(send, sendcount, MPI_DOUBLE_COMPLEX, &
                       recv, sendcount, MPI_DOUBLE_COMPLEX, MPI_COMM_X, request, ierr)
#ifndef HAVE_HIP
    !$omp end target data
#endif
    call roctxPop(range_name)

  END SUBROUTINE alltoall

  SUBROUTINE gather_full_y_line(ny, local_line, full_line)
    IMPLICIT NONE
    integer(C_INT), intent(in) :: ny
    complex(C_DOUBLE_COMPLEX), intent(in) :: local_line(ny0 - 2:nyN + 2)
    complex(C_DOUBLE_COMPLEX), intent(out) :: full_line(-1:ny + 1)

#ifdef HAVE_MPI
    complex(C_DOUBLE_COMPLEX) :: send_line(-1:ny + 1)
    integer(C_INT) :: send_start, send_end
#else
    integer(C_INT) :: send_start, send_end
#endif

    send_start = ny0
    send_end = nyN

    if (ipy == 0) send_start = -1
    if (ipy == npy_grid - 1) send_end = ny + 1

#ifdef HAVE_MPI
    send_line = (0.0d0, 0.0d0)

    send_line(send_start:send_end) = local_line(send_start:send_end)

    call roctxPush("MPI_Allreduce gather_full_y_line")
    call MPI_Allreduce(send_line, full_line, ny + 3, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_Y, ierr)
    call roctxPop("MPI_Allreduce gather_full_y_line")
#else
    full_line(-1:ny + 1) = local_line(-1:ny + 1)
#endif
  END SUBROUTINE gather_full_y_line

  !------- Divide the problem in pencils ---------!
  !-----------------------------------------------!
  SUBROUTINE init_MPI(nxpp, nz, ny, nxd, nzd, nPhi, overlapping, npy_requested)
    integer(C_INT), intent(in)  :: nxpp, nz, ny, nxd, nzd, nPhi, npy_requested
    logical, intent(in) :: overlapping
    integer, parameter :: ndims = 4
    integer :: i, color, key
    integer :: array_of_sizes(ndims), array_of_subsizes(ndims), array_of_starts(ndims), ierror
    integer(C_INT) :: miny_local, maxy_local
    type(c_ptr) :: sendptr, recvptr
    integer(c_size_t) :: sendsize, recvsize
    ! Define which process write on screen
    has_terminal = (iproc == 0)
    npy_grid = npy_requested
#ifdef HAVE_MPI
    if (npy_grid < 1) then
      if (has_terminal) print *, "Error: npy must be >= 1."
      call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
    if (mod(nproc, npy_grid) /= 0) then
      if (has_terminal) then
        print *, "Error: nproc must be divisible by npy."
        print *, "       Received nproc=", nproc, " npy=", npy_grid
      end if
      call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
#else
    if (npy_grid /= 1) error stop "init_MPI: npy > 1 requires MPI"
#endif
    npxz = nproc/npy_grid
    ipy = iproc/npxz
    ipxz = mod(iproc, npxz)

#ifdef HAVE_MPI
    color = ipy
    key = ipxz
    call MPI_Comm_split(MPI_COMM_WORLD, color, key, MPI_COMM_X, ierr)
    color = ipxz
    key = ipy
    call MPI_Comm_split(MPI_COMM_WORLD, color, key, MPI_COMM_Y, ierr)
#endif
    ! Calculate domain division in wall-normal direction.
    ! ny0:nyN are the only persistent y-partition variables.
    ny0 = 1 + ipy*(ny - 1)/npy_grid
    nyN = (ipy + 1)*(ny - 1)/npy_grid

    ! Local write extent, including physical ghost rows only on the end ranks.
    miny_local = ny0
    maxy_local = nyN
    if (ipy == 0) miny_local = ny0 - 2
    if (ipy == npy_grid - 1) maxy_local = nyN + 2
    if (npy_grid == 1) then
      miny_local = ny0 - 2
      maxy_local = nyN + 2
    end if

    ! Calculate domain division
    nx0 = ipxz*(nxpp)/npxz; nxN = (ipxz + 1)*(nxpp)/npxz - 1; nxB = nxN - nx0 + 1; 
    nz0 = ipxz*nzd/npxz; nzN = (ipxz + 1)*nzd/npxz - 1; nzB = nzN - nz0 + 1; 
    has_average = (nx0 == 0)
    !$omp target update to(npy_grid, npxz, ipy, ipxz, nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN)
    fft_transpose_is_local = (nzB == nzd)
#ifdef HAVE_MPI
#ifdef mpiverbose
    DO i = 0, nproc - 1
       IF (iproc==i) WRITE(*,*) "iproc=",iproc," ipxz=",ipxz," ipy=",ipy," nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, "ny0=", ny0, "nyN=", nyN
      CALL MPI_Barrier(MPI_COMM_WORLD)
    END DO
    FLUSH (output_unit)
#endif
    ! The pairwise all-to-all transpose uses one shared sendcount for every rank,
    ! so both decomposed dimensions must divide evenly across MPI ranks.
    if (mod(nxpp, npxz) /= 0 .or. mod(nzd, npxz) /= 0) then
      if (has_terminal) then
        print *, "Error: FFT transpose requires npxz to divide both nx+1 and nzd."
        print *, "       Received nx+1=", nxpp, " nzd=", nzd, " npxz=", npxz, " npy=", npy_grid
        print *, "       This run would create uneven transpose counts and can fail in MPI_Ialltoall."
      end if
      CALL MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
    if (int(npxz, 8)*int(nxB, 8)*int(nzB, 8)*int(nyN - ny0 + 5, 8) > huge(0_C_INT)) then
      if (has_terminal) then
        print *, "Error: problem too large for MPI transpose (integer overflow). Try to increase the number of processes."
      end if
      CALL MPI_Abort(MPI_COMM_WORLD, 1, ierror)
    end if
    sendsize = npxz*nxB*nzB*(nyN - ny0 + 5)
    recvsize = npxz*nxB*nzB*(nyN - ny0 + 5)

    sendcount = nxB*nzB*(nyN - ny0 + 5)
    !$omp target update to(sendcount)

    ! Allocate buffers for transposes*int(16, c_size_t)
#if defined(HAVE_HIP)
    ! On HIP with HSA_XNACK=1, the use_device_ptr statements around the MPI calls are ignored.
    ! Hence, MPI does a CPU mpi copy! So we need to allocate it explicity on the device.
    sendptr = omp_target_alloc(sendsize*int(16*merge(2, 1, overlapping), c_size_t), omp_get_default_device())
    recvptr = omp_target_alloc(recvsize*int(16*merge(2, 1, overlapping), c_size_t), omp_get_default_device())
    call c_f_pointer(sendptr, sendbuf, [sendsize, int(merge(2, 1, overlapping), c_size_t)])
    call c_f_pointer(recvptr, recvbuf, [recvsize, int(merge(2, 1, overlapping), c_size_t)])
#else
    ALLOCATE (sendbuf(sendsize, merge(2, 1, overlapping))); sendbuf = 0
    ALLOCATE (recvbuf(recvsize, merge(2, 1, overlapping))); recvbuf = 0
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
    array_of_subsizes = [maxy_local - miny_local + 1, 2*nz + 1, nxB, 3 + nPhi] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny_local + 1, 0, nx0, 0] ! starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, writeview_type, ierror)
    CALL MPI_Type_commit(writeview_type, ierror)
    ! For WRITING VELOCITY, SKIPPING HALO CELLS: datatype with holes to skip halo cells and select only data to be written
    array_of_sizes = [(nyN + 2) - (ny0 - 2) + 1, 2*nz + 1, nxB, 3 + nPhi] ! size along each dimension of the array IN MEMORY owned by each process
    array_of_subsizes = [maxy_local - miny_local + 1, 2*nz + 1, nxB, 3 + nPhi] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny_local - (ny0 - 2), 0, 0, 0] ! starting position of each component; !!! IT'S ZERO BASED AND WRT TO ARRAY IN MEMORY !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, owned2write_type, ierror)
    CALL MPI_Type_commit(owned2write_type, ierror)
#endif
  END SUBROUTINE init_MPI

END MODULE mpi_transpose
