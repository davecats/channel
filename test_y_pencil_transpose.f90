#include "header.h"

program test_y_pencil_transpose
  use, intrinsic :: iso_c_binding
  use mpi_transpose
  use y_line_solvers, only: ys_prepare_ghost_field_workspace, ys_solve_ghost_field_reduced, ys_rhs_store, ys_matrix_store, &
                            ys_eqm1_store, ys_eq0_store, ys_eqn_store, ys_eqnp1_store
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
  allocate (reduced_field(ny_test + 3, 2*nz_test + 1, nxB))

  local_err = 0.0d0

  call ys_prepare_ghost_field_workspace(ny_test, nz_test, nxB)
  do ix = 1, nxB
    global_x = nx0 + ix - 1
    do iz = 1, 2*nz_test + 1
      global_z = iz
      call fill_exact_line(exact, 1_C_INT, global_z, global_x)
      exact_field(:, iz, ix) = exact(-1:ny_test + 1)
      call assemble_test_line(exact, ix, iz)
    end do
  end do

  call ys_solve_ghost_field_reduced(reduced_field, ny_test, nz_test)

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

  deallocate (exact_field, reduced_field)
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

  subroutine assemble_test_line(exact_line, ix_local, iz_local)
    complex(C_DOUBLE_COMPLEX), intent(in) :: exact_line(-1:ny_test + 1)
    integer(C_INT), intent(in) :: ix_local, iz_local
    complex(C_DOUBLE_COMPLEX) :: rhs_line(-1:ny_test + 1)
    integer(C_INT) :: jy, iline

    iline = (ix_local - 1)*(2*nz_test + 1) + iz_local
    call apply_test_system(exact_line, rhs_line)
    ys_rhs_store(-1:ny_test + 1, iline) = rhs_line
    do jy = 1, ny_test - 1
      ys_matrix_store(jy, -2, iline) = -0.05d0
      ys_matrix_store(jy, -1, iline) = -0.25d0
      ys_matrix_store(jy, 0, iline) = 2.5d0
      ys_matrix_store(jy, 1, iline) = -0.25d0
      ys_matrix_store(jy, 2, iline) = -0.05d0
    end do
    ys_eqm1_store(-2, iline) = 1.0d0
    ys_eq0_store(-1, iline) = 1.0d0
    ys_eqn_store(1, iline) = 1.0d0
    ys_eqnp1_store(2, iline) = 1.0d0
  end subroutine assemble_test_line

end program test_y_pencil_transpose
