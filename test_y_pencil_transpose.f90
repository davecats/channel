#include "header.h"

program test_y_pencil_transpose
  use, intrinsic :: iso_c_binding
  use mpi_transpose
  use y_line_solvers, only: ys_solve_ghost_field_reduced
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
  complex(C_DOUBLE_COMPLEX), allocatable :: exact_field(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: dummy_field(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: reduced_field(:, :, :)
  real(C_DOUBLE) :: local_err, global_err, tol
  integer(C_INT) :: iy, iz, ix
  integer(C_INT) :: global_z, global_x
  complex(C_DOUBLE_COMPLEX) :: exact(-1:ny_test + 1)

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

  allocate (exact_field(ny_test + 3, 2*nz_test + 1, nxB))
  allocate (dummy_field(ny_test + 3, 2*nz_test + 1, nxB))
  allocate (reduced_field(ny_test + 3, 2*nz_test + 1, nxB))

  local_err = 0.0d0

  dummy_field = (0.0d0, 0.0d0)
  do ix = 1, nxB
    global_x = nx0 + ix - 1
    do iz = 1, 2*nz_test + 1
      global_z = iz
      call fill_exact_line(exact, 1_C_INT, global_z, global_x)
      do iy = 1, ny_test + 3
        exact_field(iy, iz, ix) = exact(iy - 2)
      end do
    end do
  end do

  call ys_solve_ghost_field_reduced(exact_field, dummy_field, reduced_field, ny_test, nz_test, build_field_line)

  do ix = 1, nxB
    global_x = nx0 + ix - 1
    do iz = 1, 2*nz_test + 1
      global_z = iz
      call fill_exact_line(exact, 1_C_INT, global_z, global_x)
      do iy = 1, ny_test + 3
        local_err = max(local_err, abs(reduced_field(iy, iz, ix) - exact(iy - 2)))
      end do
    end do
  end do

  call MPI_Allreduce(local_err, global_err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  tol = 1.0d-12
  if (global_err < tol) then
    if (iproc == 0) write (*, *) "Reduced ghost backend PASSED, max error = ", global_err
  else
    if (iproc == 0) write (*, *) "Reduced ghost backend FAILED, max error = ", global_err
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  deallocate (exact_field, dummy_field, reduced_field)
  call MPI_Finalize()
#endif

contains

  subroutine fill_exact_line(line, field_index, global_z, global_x)
    complex(C_DOUBLE_COMPLEX), intent(out) :: line(-1:ny_test + 1)
    integer(C_INT), intent(in) :: field_index, global_z, global_x
    integer(C_INT) :: jy

    do jy = -1, ny_test + 1
      line(jy) = cmplx(dble(jy + 3*field_index + 7*global_z + 11*global_x), &
                       dble(2*jy + 5*field_index + global_z + 13*global_x), kind=C_DOUBLE)
    end do
  end subroutine fill_exact_line

  subroutine apply_test_system(exact_line, rhs_line)
    complex(C_DOUBLE_COMPLEX), intent(in) :: exact_line(-1:ny_test + 1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: rhs_line(-1:ny_test + 1)
    integer(C_INT) :: jy

    rhs_line(-1) = exact_line(-1)
    rhs_line(0) = exact_line(0)
    rhs_line(ny_test) = exact_line(ny_test)
    rhs_line(ny_test + 1) = exact_line(ny_test + 1)
    do jy = 1, ny_test - 1
      rhs_line(jy) = (-0.05d0)*exact_line(jy - 2) + (-0.25d0)*exact_line(jy - 1) + 2.5d0*exact_line(jy) + &
                     (-0.25d0)*exact_line(jy + 1) + (-0.05d0)*exact_line(jy + 2)
    end do
  end subroutine apply_test_system

  subroutine build_field_line(global_x, global_z, src0_line, src1_line, line, a, eqm1, eq0, eqn, eqnp1, ny_line)
    integer(C_INT), intent(in) :: global_x, global_z, ny_line
    complex(C_DOUBLE_COMPLEX), intent(in) :: src0_line(:), src1_line(:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: line(-1:ny_line + 1)
    real(C_DOUBLE), intent(out) :: a(1:ny_line + 1, -2:2)
    real(C_DOUBLE), intent(out) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    complex(C_DOUBLE_COMPLEX) :: exact_line(-1:ny_line + 1)
    integer(C_INT) :: jy

    a = 0.0d0
    do jy = 1, ny_line - 1
      a(jy, -2) = -0.05d0
      a(jy, -1) = -0.25d0
      a(jy, 0) = 2.5d0
      a(jy, 1) = -0.25d0
      a(jy, 2) = -0.05d0
    end do
    a(ny_line:ny_line + 1, :) = 0.0d0

    eqm1 = 0.0d0
    eq0 = 0.0d0
    eqn = 0.0d0
    eqnp1 = 0.0d0
    eqm1(-2) = 1.0d0
    eq0(-1) = 1.0d0
    eqn(1) = 1.0d0
    eqnp1(2) = 1.0d0
    do jy = -1, ny_line + 1
      exact_line(jy) = src0_line(jy + 2)
    end do
    call apply_test_system(exact_line, line)
  end subroutine build_field_line

end program test_y_pencil_transpose
