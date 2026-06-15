#include "header.h"

program test_y_reduced_known_good_local_dst
  use, intrinsic :: iso_c_binding
  use dnsdata, only: eliminate_assembled_boundaries, nz
  use mpi_transpose, only: ierr, iproc, nproc, init_MPI, npy_grid, npxz, ipy, ny0, nyN, nx0, nxN, nxB, MPI_COMM_WORLD, &
                           MPI_Init, MPI_Comm_rank, MPI_Comm_size, MPI_Abort, MPI_Barrier, MPI_Allreduce, MPI_Finalize, &
                           MPI_DOUBLE_PRECISION, MPI_MAX
  use y_line_solvers, only: ys_prepare_assembled_workspace, ys_solve_endpoint_schur, &
                            ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, &
                            ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
                            ys_eqm1, ys_eq0, ys_eqn, ys_eqnp1

  implicit none

  integer(C_INT), parameter :: ny_test = 16
  integer(C_INT), parameter :: nz_test = 0
  integer(C_INT), parameter :: nxpp_test = 2
  integer(C_INT), parameter :: nzd_test = 1
  integer(C_INT), parameter :: npy_test = 2
  integer(C_INT), parameter :: schur_pass_counts(1) = (/2_C_INT/)

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

  call init_MPI(nxpp_test, nz_test, ny_test, nzd_test, 0_C_INT, .false., npy_test)

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

  call ys_prepare_assembled_workspace(ny_test, nz_test, ny0, nyN, 1_C_INT, nxB*(2*nz_test + 1), &
                                      schur_pass_counts)
  nz = nz_test

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

  call eliminate_assembled_boundaries(ny0, nyN, ny0 == 1_C_INT, nyN == ny_test - 1)
  call ys_solve_endpoint_schur(reduced_local, .false.)

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
    integer(C_INT) :: jy, p

    do jy = ny0, nyN
      p = (jy - ny0)*nxB*(2*nz_test + 1) + iline
      ys_gpsv_x(p) = (0.0d0, 0.0d0)
      ys_gpsv_ds(p) = (0.0d0, 0.0d0)
      ys_gpsv_dl(p) = (0.0d0, 0.0d0)
      ys_gpsv_d(p) = (0.0d0, 0.0d0)
      ys_gpsv_du(p) = (0.0d0, 0.0d0)
      ys_gpsv_dw(p) = (0.0d0, 0.0d0)
    end do

    ys_lower_ghost_rhs(iline) = (0.0d0, 0.0d0)
    ys_lower_boundary_rhs(iline) = (0.0d0, 0.0d0)
    ys_upper_boundary_rhs(iline) = (0.0d0, 0.0d0)
    ys_upper_ghost_rhs(iline) = (0.0d0, 0.0d0)

    ys_eqm1(:, iline) = 0.0d0
    ys_eq0(:, iline) = 0.0d0
    ys_eqn(:, iline) = 0.0d0
    ys_eqnp1(:, iline) = 0.0d0
  end subroutine clear_line_storage

  subroutine assemble_owned_rows_from_exact(exact_line, iline)
    complex(C_DOUBLE_COMPLEX), intent(in) :: exact_line(-1:ny_test + 1)
    integer(C_INT), intent(in) :: iline

    integer(C_INT) :: jy, p

    do jy = ny0, nyN
      p = (jy - ny0)*nxB*(2*nz_test + 1) + iline
      ys_gpsv_ds(p) = (-0.05d0, 0.0d0)
      ys_gpsv_dl(p) = (-0.25d0, 0.0d0)
      ys_gpsv_d(p) = (2.50d0, 0.0d0)
      ys_gpsv_du(p) = (-0.25d0, 0.0d0)
      ys_gpsv_dw(p) = (-0.05d0, 0.0d0)

      ys_gpsv_x(p) = ys_gpsv_ds(p)*exact_line(jy - 2) + ys_gpsv_dl(p)*exact_line(jy - 1) + ys_gpsv_d(p)*exact_line(jy) + &
                     ys_gpsv_du(p)*exact_line(jy + 1) + ys_gpsv_dw(p)*exact_line(jy + 2)
    end do

    if (ny0 == 1) then
      ys_lower_ghost_rhs(iline) = exact_line(-1)
      ys_lower_boundary_rhs(iline) = exact_line(0)

      ys_eqm1(:, iline) = 0.0d0
      ys_eqm1(-2, iline) = 1.0d0

      ys_eq0(:, iline) = 0.0d0
      ys_eq0(-1, iline) = 1.0d0
    end if

    if (nyN == ny_test - 1) then
      ys_upper_boundary_rhs(iline) = exact_line(ny_test)
      ys_upper_ghost_rhs(iline) = exact_line(ny_test + 1)

      ys_eqn(:, iline) = 0.0d0
      ys_eqn(1, iline) = 1.0d0

      ys_eqnp1(:, iline) = 0.0d0
      ys_eqnp1(2, iline) = 1.0d0
    end if
  end subroutine assemble_owned_rows_from_exact

end program test_y_reduced_known_good_local_dst
