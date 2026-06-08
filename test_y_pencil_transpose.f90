#include "header.h"

program test_y_reduced_known_good_local_dst
  use, intrinsic :: iso_c_binding
  use mpi_transpose
  use y_line_solvers, only: ys_prepare_ghost_field_workspace, ys_solve_ghost_field_reduced_const_operator, &
                            ys_local_rhs, ys_local_operator, &
                            ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
                            ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row
#ifdef HAVE_MPI
  use mpi_f08
#endif

  implicit none

  integer(C_INT), parameter :: ny_test = 16
  integer(C_INT), parameter :: nz_test = 0
  integer(C_INT), parameter :: nxpp_test = 2
  integer(C_INT), parameter :: nxd_test = 2
  integer(C_INT), parameter :: nzd_test = 1
  integer(C_INT), parameter :: npy_test = 2

  complex(C_DOUBLE_COMPLEX), allocatable :: reduced_local(:, :, :)
  complex(C_DOUBLE_COMPLEX) :: exact(-1:ny_test + 1)

  real(C_DOUBLE) :: local_err, global_err, tol
  integer(C_INT) :: ix, iz, iy, iline
  integer(C_INT) :: global_x, global_z

#ifndef HAVE_MPI
  error stop "test_y_reduced_known_good_local_dst requires MPI"
#else

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, iproc, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nproc /= 2) then
    if (iproc == 0) write (*, *) "test_y_reduced_known_good_local_dst requires exactly 2 MPI ranks"
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  call init_MPI(nxpp_test, nz_test, ny_test, nxd_test, nzd_test, 0_C_INT, .false., npy_test)

  if (npy_grid /= npy_test .or. npxz /= nproc/npy_test) then
    write (*, *) "Unexpected process-grid metadata on rank ", iproc, ": npy=", npy_grid, " npxz=", npxz
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  if (iproc == 0) then
    write (*, *) "Known-good reduced y-solver test with local-only dst"
    write (*, *) "ny_test = ", ny_test
    write (*, *) "npy     = ", npy_grid
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  write (*, '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
    "rank ", iproc, ": ipy=", ipy, " ny0=", ny0, " nyN=", nyN, " nx0=", nx0, " nxN=", nxN
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! Local-only output.
  !
  ! Important:
  !
  !   first dimension bounds are ny0:nyN
  !
  ! Therefore inside the reduced ghost solve:
  !
  !   dst_full_start = lbound(dst,1) + 2 = ny0 + 2
  !   dst_full_end   = ubound(dst,1) + 2 = nyN + 2
  !
  ! Since the gathered full-y buffer is indexed as y+2, this means:
  !
  !   reduced_local(iy, iz, ix) receives global y = iy
  !
  allocate (reduced_local(ny0:nyN, 1:2*nz_test + 1, 1:nxB))
  reduced_local = cmplx(-999.0d0, -999.0d0, kind=C_DOUBLE)

  call ys_prepare_ghost_field_workspace(ny_test, nz_test, nxB)

  do ix = 1, nxB
    global_x = nx0 + ix - 1

    do iz = 1, 2*nz_test + 1
      global_z = iz

      iline = (ix - 1)*(2*nz_test + 1) + iz

      call clear_line_storage(iline)
      call fill_exact_line(exact, global_z, global_x)
      call assemble_owned_rows_from_exact(exact, iline)
    end do
  end do

  call ys_solve_ghost_field_reduced_const_operator(reduced_local, ny_test, nz_test)

  local_err = 0.0d0

  do ix = 1, nxB
    global_x = nx0 + ix - 1

    do iz = 1, 2*nz_test + 1
      global_z = iz

      call fill_exact_line(exact, global_z, global_x)

      do iy = ny0, nyN
        local_err = max(local_err, abs(reduced_local(iy, iz, ix) - exact(iy)))
      end do
    end do
  end do

  call MPI_Allreduce(local_err, global_err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  tol = 1.0d-12

  if (global_err < tol) then
    if (iproc == 0) write (*, *) "Known-good local-dst reduced y-solver test PASSED, max error = ", global_err
  else
    if (iproc == 0) write (*, *) "Known-good local-dst reduced y-solver test FAILED, max error = ", global_err
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  deallocate (reduced_local)

  call MPI_Finalize()

#endif

contains

  pure subroutine fill_exact_line(line, global_z, global_x)
    complex(C_DOUBLE_COMPLEX), intent(out) :: line(-1:ny_test + 1)
    integer(C_INT), intent(in) :: global_z, global_x

    integer(C_INT) :: jy
    real(C_DOUBLE) :: yr, zr, xr

    zr = dble(global_z)
    xr = dble(global_x)

    do jy = -1, ny_test + 1
      yr = dble(jy)

      line(jy) = cmplx( &
                 1.0d0 + 0.37d0*yr - 0.013d0*yr*yr + 0.002d0*yr*yr*yr + 0.11d0*zr + 0.07d0*xr, &
                 -0.5d0 + 0.19d0*yr + 0.017d0*yr*yr - 0.001d0*yr*yr*yr + 0.03d0*zr - 0.09d0*xr, &
                 kind=C_DOUBLE)
    end do
  end subroutine fill_exact_line

  subroutine clear_line_storage(iline)
    integer(C_INT), intent(in) :: iline

    ys_local_rhs(:, iline) = (0.0d0, 0.0d0)
    ys_local_operator(:, :, iline) = 0.0d0

    ys_lower_ghost_rhs(iline) = (0.0d0, 0.0d0)
    ys_lower_boundary_rhs(iline) = (0.0d0, 0.0d0)
    ys_upper_boundary_rhs(iline) = (0.0d0, 0.0d0)
    ys_upper_ghost_rhs(iline) = (0.0d0, 0.0d0)

    ys_lower_ghost_row(:, iline) = 0.0d0
    ys_lower_boundary_row(:, iline) = 0.0d0
    ys_upper_boundary_row(:, iline) = 0.0d0
    ys_upper_ghost_row(:, iline) = 0.0d0
  end subroutine clear_line_storage

  subroutine assemble_owned_rows_from_exact(exact_line, iline)
    complex(C_DOUBLE_COMPLEX), intent(in) :: exact_line(-1:ny_test + 1)
    integer(C_INT), intent(in) :: iline

    integer(C_INT) :: jy

    do jy = ny0, nyN
      ys_local_operator(jy, -2, iline) = -0.05d0
      ys_local_operator(jy, -1, iline) = -0.25d0
      ys_local_operator(jy, 0, iline) = 2.50d0
      ys_local_operator(jy, 1, iline) = -0.25d0
      ys_local_operator(jy, 2, iline) = -0.05d0

      ys_local_rhs(jy, iline) = &
        cmplx(ys_local_operator(jy, -2, iline), 0.0d0, kind=C_DOUBLE)*exact_line(jy - 2) + &
        cmplx(ys_local_operator(jy, -1, iline), 0.0d0, kind=C_DOUBLE)*exact_line(jy - 1) + &
        cmplx(ys_local_operator(jy, 0, iline), 0.0d0, kind=C_DOUBLE)*exact_line(jy) + &
        cmplx(ys_local_operator(jy, 1, iline), 0.0d0, kind=C_DOUBLE)*exact_line(jy + 1) + &
        cmplx(ys_local_operator(jy, 2, iline), 0.0d0, kind=C_DOUBLE)*exact_line(jy + 2)
    end do

    if (ny0 == 1) then
      ys_lower_ghost_rhs(iline) = exact_line(-1)
      ys_lower_boundary_rhs(iline) = exact_line(0)

      ys_lower_ghost_row(:, iline) = 0.0d0
      ys_lower_ghost_row(-2, iline) = 1.0d0

      ys_lower_boundary_row(:, iline) = 0.0d0
      ys_lower_boundary_row(-1, iline) = 1.0d0
    end if

    if (nyN == ny_test - 1) then
      ys_upper_boundary_rhs(iline) = exact_line(ny_test)
      ys_upper_ghost_rhs(iline) = exact_line(ny_test + 1)

      ys_upper_boundary_row(:, iline) = 0.0d0
      ys_upper_boundary_row(1, iline) = 1.0d0

      ys_upper_ghost_row(:, iline) = 0.0d0
      ys_upper_ghost_row(2, iline) = 1.0d0
    end if
  end subroutine assemble_owned_rows_from_exact

end program test_y_reduced_known_good_local_dst
