#include "build_options.h"

! The scalar wall closure, both conditions, at machine precision.
!
! This is the test that would have caught `phinp1bc = d04n`: that row carries a
! zero coefficient on the ghost node ny+1, and four places in the closure divide
! by exactly that entry, so the Neumann leg produced NaN rather than a number.
! `ctest` was green throughout, because nothing solved a scalar line with
! Neumann rows and nothing ran the binary with -DphiNeumann.
!
! Why a quartic and not the cosh/sinh solution of the Helmholtz problem.  The
! compact stencils are built by requiring exactness on polynomials of degree 4
! (case_setup.fypp, setup_derivatives), and that holds for every row the closure
! uses:
!
!   sum_j der(iy,2,j)*p(y_j)  =  sum_j der(iy,0,j)*p''(y_j)   exactly
!   sum_j der(iy,3,j)*p(y_j)  =  24*c4                        exactly
!   sum_j d140(j)*p(y_{j-1})  =  p'(0),   sum_j d14n(j)*p(y_{ny-1+j}) = p'(2)
!
! so a quartic satisfies the *discrete* system exactly, with a right-hand side
! computed from the continuous data alone.  The solve is unique (lambda > 0
! makes the operator coercive even with Neumann at both walls), so the solver
! must return the nodal values of that quartic to round-off.  A cosh/sinh
! reference would only be recovered to the truncation error of the scheme and
! of the fourth-derivative ghost closure, i.e. against an eyeballed tolerance
! instead of against zero.
!
! Eight cases: {Dirichlet, Neumann} x {k = 0, k /= 0} x {lambda = 1,
! lambda = 3.75/deltat}.  The stiff lambda is the one the time stepper actually
! passes, and the one that makes the wall layer unresolved on this mesh.
! Dirichlet is included so a passing Neumann leg means something: the same
! assertions, the same right-hand side, only the two physical rows swapped.
program test_scalar_neumann
  use, intrinsic :: iso_c_binding
  use case_setup, only: ny, der, ni, y, d040, d140, d04n, d14n, iproc, &
                        phi0bc, phi0m1bc, phinbc, phinp1bc, free_memory
  use driver, only: initialize
  use mpi_transpose, only: ierr, MPI_Abort, MPI_Finalize, MPI_COMM_WORLD
  use compact_line_solvers, only: solve_full_line_compact
  implicit none

  ! p(y) = sum_i c(i) y^i -- any quartic will do; this one has no zero
  ! derivative at either wall, so an unenforced row cannot pass by accident.
  real(C_DOUBLE), parameter :: c(0:4) = [0.3d0, 0.7d0, -0.4d0, 0.25d0, -0.1d0]
  real(C_DOUBLE), parameter :: tol = 1.0d-12

  character(len=256) :: config_file, restart_in
  complex(C_DOUBLE_COMPLEX), allocatable :: x(:)
  real(C_DOUBLE), allocatable :: mat(:, :)
  real(C_DOUBLE) :: lower_bc(-2:2), upper_bc(-2:2)
  ! der is (iy, n, j), so der(iy, 3, :) is a stride-2 slice; copy the two ghost
  ! rows into contiguous locals rather than pass a temporary at every call.
  real(C_DOUBLE) :: lower_ghost(-2:2), upper_ghost(-2:2)
  real(C_DOUBLE) :: lambda, alpha, k2v, worst, err, g0, gn
  real(C_DOUBLE) :: rhs_lower, rhs_upper, rhs_ghost_lower, rhs_ghost_upper
  integer(C_INT) :: iy, ikind, ik, ilam
  logical :: failed
  character(len=9) :: kind_name(2) = ["Dirichlet", "Neumann  "]

  config_file = "tests/data/dns_test_npy2.in"
  restart_in = "tests/data/start_field.out"

  call initialize(config_file, restart_in)

  allocate (x(-1:ny + 1), mat(1:ny + 1, -2:2))
  failed = .false.
  lower_ghost = der(1, 3, :); upper_ghost = der(ny - 1, 3, :)

  ! The rows the production code assembles must be the rows tested below.
  call check_row("phi0bc  ", phi0bc, merge_row(d140, d040))
  call check_row("phinbc  ", phinbc, merge_row(d14n, d04n))
  call check_row("phi0m1bc", phi0m1bc, lower_ghost)
  call check_row("phinp1bc", phinp1bc, upper_ghost)

  do ikind = 1, 2
    if (ikind == 1) then
      lower_bc = d040; upper_bc = d04n          ! value rows
    else
      lower_bc = d140; upper_bc = d14n          ! gradient rows
    end if
    do ik = 1, 2
      k2v = merge(0.0d0, 0.5d0**2 + 1.0d0**2, ik == 1)   ! alfa0 = 0.5, beta0 = 1
      do ilam = 1, 2
        lambda = merge(1.0d0, 3.75d0/1.0d-3, ilam == 1)
        alpha = ni/0.025d0                      ! the stiffest Prandtl of the campaign

        ! lambda*phi - alpha*(D^2 - k^2)phi = 0, so the interior right-hand side
        ! is the compact operator applied to the exact solution -- assembled from
        ! p and p'' at the nodes, not from the matrix.
        mat = 0.0d0
        x = (0.0d0, 0.0d0)
        do iy = 1, ny - 1
          mat(iy, -2:2) = lambda*der(iy, 0, -2:2) - alpha*(der(iy, 2, -2:2) - k2v*der(iy, 0, -2:2))
          x(iy) = dcmplx(sum(der(iy, 0, -2:2)*((lambda + alpha*k2v)*p(y(iy - 2:iy + 2)) &
                                               - alpha*p2(y(iy - 2:iy + 2)))), 0.0d0)
        end do

        if (ikind == 1) then
          rhs_lower = p(y(0)); rhs_upper = p(y(ny))
        else
          rhs_lower = p1(y(0)); rhs_upper = p1(y(ny))
        end if
        ! Both ghost rows are the fourth-derivative relation, which a quartic
        ! satisfies with 24*c4 rather than with zero.
        rhs_ghost_lower = 24.0d0*c(4); rhs_ghost_upper = 24.0d0*c(4)

        call solve_full_line_compact(x, mat, lower_bc, lower_ghost, upper_bc, upper_ghost, &
                                     dcmplx(rhs_lower, 0.0d0), dcmplx(rhs_ghost_lower, 0.0d0), &
                                     dcmplx(rhs_upper, 0.0d0), dcmplx(rhs_ghost_upper, 0.0d0), ny)

        ! Comparisons are written as `.not. (a <= b)`, never as max() or `>`.
        ! The failure this test exists for is a division by zero, so every
        ! quantity below can be NaN -- and gfortran's max(x, NaN) returns x, so
        ! a max-based reduction silently reports 0.0 and passes.  The negated
        ! form is true for NaN, which is what a test needs.
        worst = 0.0d0
        do iy = -1, ny + 1
          err = abs(x(iy) - dcmplx(p(y(iy)), 0.0d0))
          if (.not. (err <= worst)) worst = err
        end do
        ! The gradients the wall rows deliver, read back off the solution the
        ! way runtime_diagnostics reads them.
        g0 = sum(d140(-2:2)*dreal(x(-1:3)))
        gn = sum(d14n(-2:2)*dreal(x(ny - 3:ny + 1)))
        err = worst
        if (.not. (abs(g0 - p1(y(0))) <= err)) err = abs(g0 - p1(y(0)))
        if (.not. (abs(gn - p1(y(ny))) <= err)) err = abs(gn - p1(y(ny)))

        if (.not. (err <= tol)) failed = .true.
        if (iproc == 0) write (*, "(A,A,A,F6.3,A,ES9.2,A,ES10.3)") &
          "  ", kind_name(ikind), "  k2 = ", k2v, "  lambda = ", lambda, "  max error = ", err
      end do
    end do
  end do

  if (.not. failed) then
    if (iproc == 0) write (*, *) "Scalar wall closure test PASSED"
  else
    if (iproc == 0) write (*, *) "Scalar wall closure test FAILED, tolerance = ", tol
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  call free_memory(.true.)
  call MPI_Finalize()

contains

  elemental real(C_DOUBLE) function p(yy)
    real(C_DOUBLE), intent(in) :: yy
    p = c(0) + yy*(c(1) + yy*(c(2) + yy*(c(3) + yy*c(4))))
  end function p

  elemental real(C_DOUBLE) function p1(yy)
    real(C_DOUBLE), intent(in) :: yy
    p1 = c(1) + yy*(2*c(2) + yy*(3*c(3) + yy*4*c(4)))
  end function p1

  elemental real(C_DOUBLE) function p2(yy)
    real(C_DOUBLE), intent(in) :: yy
    p2 = 2*c(2) + yy*(6*c(3) + yy*12*c(4))
  end function p2

  ! The row the build in hand is expected to use: the first under phiNeumann,
  ! the second otherwise.
  function merge_row(neumann_row, dirichlet_row) result(expected)
    real(C_DOUBLE), intent(in) :: neumann_row(-2:2), dirichlet_row(-2:2)
    real(C_DOUBLE) :: expected(-2:2)
#ifdef phiNeumann
    expected = neumann_row
#else
    expected = dirichlet_row
#endif
  end function merge_row

  subroutine check_row(name, actual, expected)
    character(len=*), intent(in) :: name
    real(C_DOUBLE), intent(in) :: actual(-2:2), expected(-2:2)
    if (maxval(abs(actual - expected)) > 0.0d0) then
      if (iproc == 0) write (*, *) "Scalar wall row ", name, " is not the one this test pins"
      failed = .true.
    end if
  end subroutine check_row

end program test_scalar_neumann
