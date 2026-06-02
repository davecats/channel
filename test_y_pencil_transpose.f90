#include "header.h"

program test_y_pencil_transpose
  use, intrinsic :: iso_c_binding
  use mpi_transpose
  use y_line_solvers, only: ys_solve_ghost_system, ys_solve_ghost_field_with_y_pencil, ys_solve_ghost_field_reduced
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
  complex(C_DOUBLE_COMPLEX), allocatable :: rhs_xz(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: rhs_y(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: sol_y(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: sol_xz(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: exact_field(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: dummy_field(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: gathered_field(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: reduced_field(:, :, :)
  real(C_DOUBLE) :: local_err, global_err, tol
  integer(C_INT) :: iy, iz, ix
  integer(C_INT) :: global_y, global_z, global_x
  integer(C_INT), parameter :: nsolve_fields = 2
  integer(C_INT) :: ifield
  real(C_DOUBLE) :: amat(1:ny_test + 1, -2:2)
  real(C_DOUBLE) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
  complex(C_DOUBLE_COMPLEX) :: work(-1:ny_test + 1)
  complex(C_DOUBLE_COMPLEX) :: exact(-1:ny_test + 1)
  complex(C_DOUBLE_COMPLEX) :: rhs_line(-1:ny_test + 1)
  real(C_DOUBLE) :: amat_work(1:ny_test + 1, -2:2)
  real(C_DOUBLE) :: eqm1_work(-2:2), eq0_work(-2:2), eqn_work(-2:2), eqnp1_work(-2:2)

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
  allocate (rhs_xz(ylB, 2*nz_test + 1, nxB, nsolve_fields))
  allocate (rhs_y(ny_test + 3, zpyB, nxB, nsolve_fields))
  allocate (sol_y(ny_test + 3, zpyB, nxB, nsolve_fields))
  allocate (sol_xz(ylB, 2*nz_test + 1, nxB, nsolve_fields))
  allocate (exact_field(ny_test + 3, 2*nz_test + 1, nxB))
  allocate (dummy_field(ny_test + 3, 2*nz_test + 1, nxB))
  allocate (gathered_field(ny_test + 3, 2*nz_test + 1, nxB))
  allocate (reduced_field(ny_test + 3, 2*nz_test + 1, nxB))

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

  call build_test_system(amat, eqm1, eq0, eqn, eqnp1)
  rhs_xz = (0.0d0, 0.0d0)
  do ifield = 1, nsolve_fields
    do ix = 1, nxB
      global_x = nx0 + ix - 1
      do iz = 1, 2*nz_test + 1
        global_z = iz
        call fill_exact_line(exact, ifield, global_z, global_x)
        call apply_test_system(exact, rhs_line)
        do iy = 1, ylB
          global_y = yl0 + iy - 1
          rhs_xz(iy, iz, ix, ifield) = rhs_line(global_y - 2)
        end do
      end do
    end do
  end do

  call transpose_xz_to_y_pencil_fields(rhs_xz, rhs_y)
  sol_y = rhs_y
  do ifield = 1, nsolve_fields
    do ix = 1, nxB
      do iz = 1, zpyB
        do iy = -1, ny_test + 1
          work(iy) = sol_y(iy + 2, iz, ix, ifield)
        end do
        amat_work = amat
        eqm1_work = eqm1
        eq0_work = eq0
        eqn_work = eqn
        eqnp1_work = eqnp1
        call ys_solve_ghost_system(work, amat_work, eqm1_work, eq0_work, eqn_work, eqnp1_work, ny_test)
        do iy = -1, ny_test + 1
          sol_y(iy + 2, iz, ix, ifield) = work(iy)
        end do
      end do
    end do
  end do

  call transpose_y_pencil_to_xz_fields(sol_y, sol_xz)
  do ifield = 1, nsolve_fields
    do ix = 1, nxB
      global_x = nx0 + ix - 1
      do iz = 1, 2*nz_test + 1
        global_z = iz
        call fill_exact_line(exact, ifield, global_z, global_x)
        do iy = 1, ylB
          global_y = yl0 + iy - 1
          local_err = max(local_err, abs(sol_xz(iy, iz, ix, ifield) - exact(global_y - 2)))
        end do
      end do
    end do
  end do

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

  call ys_solve_ghost_field_with_y_pencil(exact_field, dummy_field, gathered_field, ny_test, nz_test, build_field_line)
  call ys_solve_ghost_field_reduced(exact_field, dummy_field, reduced_field, ny_test, nz_test, build_field_line)

  do ix = 1, nxB
    global_x = nx0 + ix - 1
    do iz = 1, 2*nz_test + 1
      global_z = iz
      call fill_exact_line(exact, 1_C_INT, global_z, global_x)
      do iy = 1, ny_test + 3
        local_err = max(local_err, abs(gathered_field(iy, iz, ix) - exact(iy - 2)))
        local_err = max(local_err, abs(reduced_field(iy, iz, ix) - exact(iy - 2)))
        local_err = max(local_err, abs(reduced_field(iy, iz, ix) - gathered_field(iy, iz, ix)))
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

  deallocate (xz_field, y_pencil, xz_roundtrip, rhs_xz, rhs_y, sol_y, sol_xz)
  deallocate (exact_field, dummy_field, gathered_field, reduced_field)
  call MPI_Finalize()
#endif

contains

  complex(C_DOUBLE_COMPLEX) function expected_value(global_y, global_z, global_x)
    integer(C_INT), intent(in) :: global_y, global_z, global_x

    expected_value = cmplx(dble(global_y + 100*global_z + 10000*global_x), &
                           dble(1000000 + global_y + 10*global_z + global_x), kind=C_DOUBLE)
  end function expected_value

  subroutine build_test_system(a, eqm1, eq0, eqn, eqnp1)
    real(C_DOUBLE), intent(out) :: a(1:ny_test + 1, -2:2)
    real(C_DOUBLE), intent(out) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)
    integer(C_INT) :: jy

    a = 0.0d0
    do jy = 1, ny_test - 1
      a(jy, -2) = -0.05d0
      a(jy, -1) = -0.25d0
      a(jy, 0) = 2.5d0
      a(jy, 1) = -0.25d0
      a(jy, 2) = -0.05d0
    end do
    a(ny_test:ny_test + 1, :) = 0.0d0

    eqm1 = 0.0d0
    eq0 = 0.0d0
    eqn = 0.0d0
    eqnp1 = 0.0d0
    eqm1(-2) = 1.0d0
    eq0(-1) = 1.0d0
    eqn(1) = 1.0d0
    eqnp1(2) = 1.0d0
  end subroutine build_test_system

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

    call build_test_system(a, eqm1, eq0, eqn, eqnp1)
    do jy = -1, ny_line + 1
      exact_line(jy) = src0_line(jy + 2)
    end do
    call apply_test_system(exact_line, line)
  end subroutine build_field_line

end program test_y_pencil_transpose
