#include "header.h"

program test_y_pencil_transpose
  use, intrinsic :: iso_c_binding
  use mpi_transpose
#ifdef HAVE_MPI
  use mpi_f08
#endif
  implicit none

  integer(C_INT), parameter :: nx_test = 5
  integer(C_INT), parameter :: ny_test = 16
  integer(C_INT), parameter :: nz_test = 8
  integer(C_INT), parameter :: nxd_test = 9
  integer(C_INT), parameter :: nzd_test = 24
  integer(C_INT), parameter :: npy_test = 2
  complex(C_DOUBLE_COMPLEX), allocatable :: xz_field(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: y_pencil(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: xz_roundtrip(:, :, :)
  real(C_DOUBLE) :: local_err, global_err, tol
  integer(C_INT) :: iy, iz, ix
  integer(C_INT) :: global_y, global_z, global_x

#ifndef HAVE_MPI
  error stop "test_y_pencil_transpose requires MPI"
#else
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nproc /= 2) then
    if (iproc == 0) write (*, *) "test_y_pencil_transpose requires 2 MPI ranks"
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  call init_MPI(nx_test + 1, nz_test, ny_test, nxd_test, nzd_test, 0_C_INT, .false., npy_test)

  if (npy_grid /= npy_test .or. npxz /= nproc/npy_test) then
    write (*, *) "Unexpected process-grid metadata on rank ", iproc, ": npy=", npy_grid, " npxz=", npxz
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  allocate (xz_field(ylB, 2*nz_test + 1, nxB))
  allocate (y_pencil(ny_test + 3, zpyB, nxB))
  allocate (xz_roundtrip(ylB, 2*nz_test + 1, nxB))

  do ix = 1, nxB
    global_x = nx0 + ix - 1
    do iz = 1, 2*nz_test + 1
      global_z = iz
      do iy = 1, ylB
        global_y = yl0 + iy - 1
        xz_field(iy, iz, ix) = expected_value(global_y, global_z, global_x)
      end do
    end do
  end do

  y_pencil = (0.0d0, 0.0d0)
  xz_roundtrip = (0.0d0, 0.0d0)
  call transpose_xz_to_y_pencil(xz_field, y_pencil)

  local_err = 0.0d0
  do ix = 1, nxB
    global_x = nx0 + ix - 1
    do iz = 1, zpyB
      global_z = zpy0 + iz - 1
      do iy = 1, ny_test + 3
        local_err = max(local_err, abs(y_pencil(iy, iz, ix) - expected_value(iy, global_z, global_x)))
      end do
    end do
  end do

  call transpose_y_pencil_to_xz(y_pencil, xz_roundtrip)

  do ix = 1, nxB
    do iz = 1, 2*nz_test + 1
      do iy = 1, ylB
        local_err = max(local_err, abs(xz_roundtrip(iy, iz, ix) - xz_field(iy, iz, ix)))
      end do
    end do
  end do

  call MPI_Allreduce(local_err, global_err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  tol = 1.0d-12
  if (global_err < tol) then
    if (iproc == 0) write (*, *) "Y-pencil transpose roundtrip PASSED, max error = ", global_err
  else
    if (iproc == 0) write (*, *) "Y-pencil transpose roundtrip FAILED, max error = ", global_err
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  deallocate (xz_field, y_pencil, xz_roundtrip)
  call MPI_Finalize()
#endif

contains

  complex(C_DOUBLE_COMPLEX) function expected_value(global_y, global_z, global_x)
    integer(C_INT), intent(in) :: global_y, global_z, global_x

    expected_value = cmplx(dble(global_y + 100*global_z + 10000*global_x), &
                           dble(1000000 + global_y + 10*global_z + global_x), kind=C_DOUBLE)
  end function expected_value

end program test_y_pencil_transpose
