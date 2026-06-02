#include "header.h"

program test_xz_transpose
  use, intrinsic :: iso_c_binding
  use mpi_transpose
#ifdef HAVE_MPI
  use mpi_f08
#endif
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  use omp_lib
#endif
  implicit none

  integer(C_INT), parameter :: nx_test = 5
  integer(C_INT), parameter :: ny_test = 16
  integer(C_INT), parameter :: nz_test = 8
  integer(C_INT), parameter :: nxd_test = 9
  integer(C_INT), parameter :: nzd_test = 24
  complex(C_DOUBLE_COMPLEX), allocatable :: vz(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: vx(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: vz_back(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: vx_expected(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: vz_expected(:, :, :)
  real(C_DOUBLE) :: local_err_ztox, local_err_xtoz
  real(C_DOUBLE) :: global_err_ztox, global_err_xtoz, tol
  integer(C_INT) :: ix, iz, iy, src
#ifdef HAVE_MPI
  type(MPI_Request) :: request
  type(MPI_Status) :: status
#endif
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  integer :: num_dev, dev
#endif

#ifndef HAVE_MPI
  error stop "test_xz_transpose requires MPI"
#else
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nproc /= 2) then
    if (iproc == 0) write (*, *) "test_xz_transpose requires 2 MPI ranks"
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  num_dev = omp_get_num_devices()
  if (num_dev > 0) then
    dev = mod(iproc, num_dev)
    call omp_set_default_device(dev)
  end if
#endif

  call init_MPI(nx_test + 1, nz_test, ny_test, nxd_test + 1, nzd_test, 0_C_INT, .false., 1_C_INT)

  allocate (vz(nzd_test, nxB, ny_test + 3))
  allocate (vx(nx_test + 1, nzB, ny_test + 3))
  allocate (vz_back(nzd_test, nxB, ny_test + 3))
  allocate (vx_expected(nx_test + 1, nzB, ny_test + 3))
  allocate (vz_expected(nzd_test, nxB, ny_test + 3))

  do iy = 1, ny_test + 3
    do ix = 1, nxB
      do iz = 1, nzd_test
        vz(iz, ix, iy) = make_value(iy, iz, ix + nx0)
        vz_expected(iz, ix, iy) = vz(iz, ix, iy)
      end do
    end do
  end do

  do iy = 1, ny_test + 3
    do iz = 1, nzB
      do ix = 1, nx_test + 1
        src = (ix - 1)/nxB
        vx_expected(ix, iz, iy) = make_value(iy, iz + ipxz*nzB, modulo(ix - 1, nxB) + 1 + src*nxB)
      end do
    end do
  end do

  vx = (0.0d0, 0.0d0)
  vz_back = (0.0d0, 0.0d0)

  !$omp target enter data map(to: vz, vx_expected, vz_expected) map(alloc: vx, vz_back)

  call pack_zTOx(vz, sendbuf(:, 1), ny_test)
  call alltoall(sendbuf(:, 1), recvbuf(:, 1), request)
  call MPI_Wait(request, status, ierr)
  call unpack_zTOx(recvbuf(:, 1), vx, ny_test)

  call pack_xTOz(vx, sendbuf(:, 1), ny_test)
  call alltoall(sendbuf(:, 1), recvbuf(:, 1), request)
  call MPI_Wait(request, status, ierr)
  call unpack_xTOz(recvbuf(:, 1), vz_back, ny_test)

  !$omp target update from(vx, vz_back)
  !$omp target exit data map(delete: vz, vx, vz_back, vx_expected, vz_expected)

  local_err_ztox = 0.0d0
  local_err_xtoz = 0.0d0
  do iy = 1, ny_test + 3
    do iz = 1, nzB
      do ix = 1, nx_test + 1
        local_err_ztox = max(local_err_ztox, abs(vx(ix, iz, iy) - vx_expected(ix, iz, iy)))
      end do
    end do
    do ix = 1, nxB
      do iz = 1, nzd_test
        local_err_xtoz = max(local_err_xtoz, abs(vz_back(iz, ix, iy) - vz_expected(iz, ix, iy)))
      end do
    end do
  end do

  call MPI_Allreduce(local_err_ztox, global_err_ztox, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(local_err_xtoz, global_err_xtoz, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  tol = 1.0d-12
  if (iproc == 0) then
    write (*, *) "xz transpose max error z->x = ", global_err_ztox
    write (*, *) "xz transpose max error x->z = ", global_err_xtoz
  end if

  if (global_err_ztox > tol .or. global_err_xtoz > tol) then
    if (iproc == 0) write (*, *) "xz transpose FAILED"
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  if (iproc == 0) write (*, *) "xz transpose PASSED"

  deallocate (vz, vx, vz_back, vx_expected, vz_expected)
  call MPI_Finalize()
#endif

contains

  pure complex(C_DOUBLE_COMPLEX) function make_value(iy, iz, ix_global)
    integer(C_INT), intent(in) :: iy, iz, ix_global

    make_value = cmplx(dble(1000*ix_global + 100*iz + iy), &
                       dble(-2000*ix_global + 10*iz - iy), kind=C_DOUBLE)
  end function make_value

end program test_xz_transpose
