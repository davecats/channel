program test_mean_correction
  use, intrinsic :: iso_c_binding
  use dnsdata
  use driver
  use mpi_transpose, only: ipy, npy_grid, MPI_COMM_Y, ierr
  use compact_line_solvers, only: solve_full_line_compact
  implicit none

  character(len=256) :: config_file, restart_in
  complex(C_DOUBLE_COMPLEX), allocatable :: x_ref(:), x_test(:)
  real(C_DOUBLE), allocatable :: mat_ref(:, :), mat_test(:, :)
  real(C_DOUBLE) :: local_err, global_err, tol
  integer(C_INT) :: iy

  config_file = "tests/data/dns_test_npy2.in"
  restart_in = "tests/data/start_field.out"

  call initialize(config_file, restart_in)

  allocate (x_ref(-1:ny + 1), x_test(-1:ny + 1))
  allocate (mat_ref(1:ny + 1, -2:2), mat_test(1:ny + 1, -2:2))

  x_ref(-1:0) = 0.0d0
  x_ref(1:ny - 1) = 1.0d0
  x_ref(ny:ny + 1) = 0.0d0
  x_test = x_ref

  mat_ref = 0.0d0
  do iy = 1, ny - 1
    mat_ref(iy, -2:2) = der(iy, 0, -2:2) - ni*(der(iy, 2, -2:2) - k2(0, 0)*der(iy, 0, -2:2))
  end do
  mat_test = mat_ref

  call solve_full_line_compact(x_ref, mat_ref, eta0bc, eta0m1bc, etanbc, etanp1bc, &
                               (0.0d0, 0.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0), ny)

  if (npy_grid == 1 .or. ipy == 0) then
    call solve_full_line_compact(x_test, mat_test, eta0bc, eta0m1bc, etanbc, etanp1bc, &
                                 (0.0d0, 0.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0), ny)
  end if
#ifdef HAVE_MPI
  if (npy_grid > 1) call MPI_Bcast(x_test, ny + 3, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_Y, ierr)
#endif

  local_err = maxval(abs(x_ref - x_test))
#ifdef HAVE_MPI
  call MPI_Allreduce(local_err, global_err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
  global_err = local_err
#endif

  tol = 1.0e-12_C_DOUBLE
  if (global_err < tol) then
    write (*, *) "Mean correction test PASSED, max error = ", global_err
  else
    write (*, *) "Mean correction test FAILED, max error = ", global_err
#ifdef HAVE_MPI
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  end if

  call free_memory(.true.)
#ifdef HAVE_MPI
  call MPI_Finalize()
#endif
end program test_mean_correction
