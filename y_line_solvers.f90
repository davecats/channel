#include "header.h"

module y_line_solvers

  use, intrinsic :: iso_c_binding
  use mpi_transpose, only: ny0, nyN, nx0, nxN, nxB, npy_grid, ierr, ipy, MPI_COMM_Y
  use roctx, only: roctxPush, roctxPop
#ifdef HAVE_MPI
  use mpi_f08
#endif
#ifdef HAVE_CUDA
  use cusparse
#endif

  implicit none
  private
  integer(C_INT), parameter :: YS_ENDPOINT_RESPONSE_FULL = 0_C_INT
  integer(C_INT), parameter :: YS_ENDPOINT_RESPONSE_CONST = 1_C_INT
  integer(C_INT), parameter :: YS_ENDPOINT_RESPONSE_SYMMETRIC = 2_C_INT

  public :: ys_lu5decomp, ys_leftlu5div
  public :: ys_solve_compact_derivative, ys_solve_compact_system, ys_solve_ghost_system
  public :: ys_prepare_ghost_field_workspace, ys_release_ghost_field_workspace, ys_solve_ghost_field_reduced
  public :: ys_local_rhs, ys_local_operator
  public :: ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs
  public :: ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row
  public :: ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x
  public :: ys_prepare_gpsv_workspace, ys_solve_packed_pentadiagonal
  public :: ys_solve_ghost_field_reduced_const_operator
  public :: ys_solve_ghost_field_reduced_symmetric_operator
  public :: ys_solve_full_y_lines_packed

  integer(C_INT), save :: ys_workspace_ny = -1
  integer(C_INT), save :: ys_workspace_nz = -1
  integer(C_INT), save :: ys_workspace_nx = -1
  integer(C_INT), save :: ys_workspace_nlines = 0
  integer(C_INT), save :: ys_workspace_active_n = 0
  integer(C_INT), save :: ys_workspace_npy = -1
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_local_rhs(:, :)
  real(C_DOUBLE), allocatable, save :: ys_local_operator(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_lower_ghost_rhs(:), ys_lower_boundary_rhs(:), ys_upper_boundary_rhs(:), ys_upper_ghost_rhs(:)
  real(C_DOUBLE), allocatable, save :: ys_lower_ghost_row(:, :), ys_lower_boundary_row(:, :), ys_upper_boundary_row(:, :), ys_upper_ghost_row(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_boundary_lower_rhs0(:), ys_boundary_upper_rhsn(:)
  real(C_DOUBLE), allocatable, save :: ys_boundary_lower_eq(:, :), ys_boundary_upper_eq(:, :)
  real(C_DOUBLE), allocatable, save :: ys_interior_lu(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_interior_response_columns(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rows_send(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_left_interface_values(:, :), ys_right_interface_values(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rows_recv(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_matrix_lu(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_reduced_rhs(:, :)
#ifdef HAVE_CUDA
  type(cusparseHandle), save :: ys_gpsv_handle
  logical, save :: ys_gpsv_handle_created = .false.
#endif
  integer(C_INT), save :: ys_gpsv_n = -1, ys_gpsv_batch = -1
#ifdef HAVE_CUDA
  integer(8), save :: ys_gpsv_buffer_size = 0_8
#endif
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ys_gpsv_ds(:), ys_gpsv_dl(:), ys_gpsv_d(:), ys_gpsv_du(:), ys_gpsv_dw(:), ys_gpsv_x(:)
#ifdef HAVE_CUDA
  character(c_char), allocatable, save :: ys_gpsv_buffer(:)
#endif

contains

  subroutine ys_prepare_ghost_field_workspace(ny, nz, nx_lines)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, nx_lines
    integer(C_INT) :: row_start, row_end, active_n, nlines

    ! The reduced interface solve keeps two boundary-adjacent unknowns per side,
    ! so each rank must own at least four physical/ghost rows.
    if (npy_grid > 1 .and. nyN - ny0 + 1 < 4) error stop "ys_solve_ghost_field requires at least four y rows per rank"
    row_start = ny0
    row_end = nyN
    active_n = row_end - row_start + 1
    if (active_n < 4) error stop "ys_solve_ghost_field requires at least four active y rows"

    nlines = nx_lines*(2*nz + 1)
    if (allocated(ys_local_rhs)) then
      if (ys_workspace_ny /= ny .or. ys_workspace_nz /= nz .or. ys_workspace_nx /= nx_lines .or. &
          ys_workspace_active_n /= active_n .or. ys_workspace_npy /= npy_grid) then
        !$omp target exit data map(delete: ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
        !$omp& ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ys_boundary_lower_rhs0, ys_boundary_upper_rhsn, &
        !$omp& ys_boundary_lower_eq, ys_boundary_upper_eq, ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, &
        !$omp& ys_right_interface_values, ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
  deallocate (ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs)
        deallocate (ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row)
        deallocate (ys_boundary_lower_rhs0, ys_boundary_upper_rhsn, ys_boundary_lower_eq, ys_boundary_upper_eq)
deallocate (ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values)
        deallocate (ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
      end if
    end if

    if (.not. allocated(ys_local_rhs)) then
      allocate (ys_local_rhs(ny0:nyN, nlines), ys_local_operator(ny0:nyN, -2:2, nlines))
     allocate (ys_lower_ghost_rhs(nlines), ys_lower_boundary_rhs(nlines), ys_upper_boundary_rhs(nlines), ys_upper_ghost_rhs(nlines))
      allocate (ys_lower_ghost_row(-2:2, nlines), ys_lower_boundary_row(-2:2, nlines), ys_upper_boundary_row(-2:2, nlines), ys_upper_ghost_row(-2:2, nlines))
      allocate (ys_boundary_lower_rhs0(nlines), ys_boundary_upper_rhsn(nlines), ys_boundary_lower_eq(-1:2, nlines), ys_boundary_upper_eq(-2:1, nlines))
      allocate (ys_interior_lu(0:active_n - 1, -2:2, nlines), ys_interior_response_columns(0:active_n - 1, 5, nlines))
      allocate (ys_reduced_rows_send(20, nlines), ys_left_interface_values(2, nlines), ys_right_interface_values(2, nlines))
      allocate (ys_reduced_rows_recv(20, nlines, npy_grid), ys_reduced_matrix_lu(4*npy_grid, 11, nlines), ys_reduced_rhs(4*npy_grid, nlines))
      !$omp target enter data map(alloc: ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
      !$omp& ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ys_boundary_lower_rhs0, ys_boundary_upper_rhsn, &
      !$omp& ys_boundary_lower_eq, ys_boundary_upper_eq, ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, &
      !$omp& ys_right_interface_values, ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
    end if

    ys_workspace_ny = ny
    ys_workspace_nz = nz
    ys_workspace_nx = nx_lines
    ys_workspace_nlines = nlines
    ys_workspace_active_n = active_n
    ys_workspace_npy = npy_grid
    if (npy_grid == 1) call ys_prepare_gpsv_workspace(active_n, nlines)

    ys_local_rhs = (0.0d0, 0.0d0)
    ys_local_operator = 0.0d0
    ys_lower_ghost_rhs = (0.0d0, 0.0d0)
    ys_lower_boundary_rhs = (0.0d0, 0.0d0)
    ys_upper_boundary_rhs = (0.0d0, 0.0d0)
    ys_upper_ghost_rhs = (0.0d0, 0.0d0)
    ys_lower_ghost_row = 0.0d0
    ys_lower_boundary_row = 0.0d0
    ys_upper_boundary_row = 0.0d0
    ys_upper_ghost_row = 0.0d0
    ys_boundary_lower_rhs0 = (0.0d0, 0.0d0)
    ys_boundary_upper_rhsn = (0.0d0, 0.0d0)
    ys_boundary_lower_eq = 0.0d0
    ys_boundary_upper_eq = 0.0d0
    ys_interior_lu = 0.0d0
    ys_interior_response_columns = (0.0d0, 0.0d0)
    ys_reduced_rows_send = (0.0d0, 0.0d0)
    ys_left_interface_values = (0.0d0, 0.0d0)
    ys_right_interface_values = (0.0d0, 0.0d0)
    ys_reduced_rows_recv = (0.0d0, 0.0d0)
    ys_reduced_matrix_lu = (0.0d0, 0.0d0)
    ys_reduced_rhs = (0.0d0, 0.0d0)
  end subroutine ys_prepare_ghost_field_workspace

  subroutine ys_release_ghost_field_workspace()
    implicit none

    if (.not. allocated(ys_local_rhs)) return

    !$omp target exit data map(delete: ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
    !$omp& ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ys_boundary_lower_rhs0, ys_boundary_upper_rhsn, &
    !$omp& ys_boundary_lower_eq, ys_boundary_upper_eq, ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, &
    !$omp& ys_right_interface_values, ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
  deallocate (ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs)
    deallocate (ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row)
    deallocate (ys_boundary_lower_rhs0, ys_boundary_upper_rhsn, ys_boundary_lower_eq, ys_boundary_upper_eq)
deallocate (ys_interior_lu, ys_interior_response_columns, ys_reduced_rows_send, ys_left_interface_values, ys_right_interface_values)
    deallocate (ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs)
    call ys_release_gpsv_workspace()

    ys_workspace_ny = -1
    ys_workspace_nz = -1
    ys_workspace_nx = -1
    ys_workspace_nlines = 0
    ys_workspace_active_n = 0
    ys_workspace_npy = -1
  end subroutine ys_release_ghost_field_workspace

#ifdef HAVE_CUDA
  subroutine ys_check_cusparse(status, where)
    integer(C_INT), intent(in) :: status
    character(*), intent(in) :: where

    if (status /= CUSPARSE_STATUS_SUCCESS) then
      print *, "cuSPARSE error in ", trim(where), ": status=", status
      error stop
    end if
  end subroutine ys_check_cusparse
#endif

  subroutine ys_release_gpsv_workspace()
    implicit none
#ifdef HAVE_CUDA
    integer(C_INT) :: status
#endif

    if (allocated(ys_gpsv_ds)) then
      !$omp target exit data map(delete: ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x)
      deallocate (ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x)
    end if
#ifdef HAVE_CUDA
    if (allocated(ys_gpsv_buffer)) then
      !$omp target exit data map(delete: ys_gpsv_buffer)
      deallocate (ys_gpsv_buffer)
    end if
    if (ys_gpsv_handle_created) then
      status = cusparseDestroy(ys_gpsv_handle)
      call ys_check_cusparse(status, "cusparseDestroy")
      ys_gpsv_handle_created = .false.
    end if
#endif
    ys_gpsv_n = -1
    ys_gpsv_batch = -1
#ifdef HAVE_CUDA
    ys_gpsv_buffer_size = 0_8
#endif
  end subroutine ys_release_gpsv_workspace

  subroutine ys_prepare_gpsv_workspace(n, batch_count)
    implicit none
    integer(C_INT), intent(in) :: n, batch_count
#ifdef HAVE_CUDA
    integer(C_INT) :: status
    integer(8) :: buffer_size
#endif

    if (ys_gpsv_n == n .and. ys_gpsv_batch == batch_count .and. allocated(ys_gpsv_ds)) return

    call ys_release_gpsv_workspace()

#ifdef HAVE_CUDA
    status = cusparseCreate(ys_gpsv_handle)
    call ys_check_cusparse(status, "cusparseCreate")
    ys_gpsv_handle_created = .true.
#endif

    allocate (ys_gpsv_ds(n*batch_count), ys_gpsv_dl(n*batch_count), ys_gpsv_d(n*batch_count), &
              ys_gpsv_du(n*batch_count), ys_gpsv_dw(n*batch_count), ys_gpsv_x(n*batch_count))
    !$omp target enter data map(alloc: ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x)
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x)
    status = cusparseZgpsvInterleavedBatch_bufferSize(ys_gpsv_handle, 0_C_INT, n, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, &
                                                      ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, batch_count, buffer_size)
    !$omp end target data
    call ys_check_cusparse(status, "cusparseZgpsvInterleavedBatch_bufferSize")
    ys_gpsv_buffer_size = buffer_size
    allocate (ys_gpsv_buffer(max(1, int(buffer_size))))
    !$omp target enter data map(alloc: ys_gpsv_buffer)
#endif

    ys_gpsv_n = n
    ys_gpsv_batch = batch_count
  end subroutine ys_prepare_gpsv_workspace

  subroutine ys_solve_packed_pentadiagonal(n, batch_count, label)
    implicit none
    integer(C_INT), intent(in) :: n, batch_count
    character(*), intent(in) :: label
#ifdef HAVE_CUDA
    integer(C_INT) :: status
#else
    integer(C_INT) :: iline
#endif

    call roctxPush(label)
#ifdef HAVE_CUDA
    !$omp target data use_device_addr(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, ys_gpsv_buffer)
    status = cusparseZgpsvInterleavedBatch(ys_gpsv_handle, 0_C_INT, n, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, &
                                           ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, batch_count, ys_gpsv_buffer)
    !$omp end target data
    call ys_check_cusparse(status, "cusparseZgpsvInterleavedBatch "//trim(label))
#else
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, n, batch_count) &
    !$omp private(iline)
    do iline = 1, batch_count
      call ys_factor_penta_interleaved(ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, batch_count, iline, n)
      call ys_solve_factored_penta_interleaved(ys_gpsv_x, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, &
                                               batch_count, iline, n)
    end do
    !$omp end target teams distribute parallel do
#endif
    call roctxPop(label)
  end subroutine ys_solve_packed_pentadiagonal

  !$omp declare target(ys_compact_physical_boundary_equations)
  subroutine ys_compact_physical_boundary_equations(lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, lower_boundary_value, &
                                                lower_ghost_value, upper_boundary_value, upper_ghost_value, lower_rhs0, lower_eq0, &
                                                    upper_rhsn, upper_eqn)
    implicit none
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: lower_boundary_value, lower_ghost_value, upper_boundary_value, upper_ghost_value
    complex(C_DOUBLE_COMPLEX), intent(out) :: lower_rhs0, upper_rhsn
    real(C_DOUBLE), intent(out) :: lower_eq0(-2:2), upper_eqn(-2:2)

    lower_rhs0 = lower_boundary_value - lower_ghost_value*lower_bc(-2)/lower_ghost_bc(-2)
    lower_eq0 = lower_bc - lower_ghost_bc*lower_bc(-2)/lower_ghost_bc(-2)
    lower_eq0(-2) = 0.0d0
    upper_rhsn = upper_boundary_value - upper_ghost_value*upper_bc(2)/upper_ghost_bc(2)
    upper_eqn = upper_bc - upper_ghost_bc*upper_bc(2)/upper_ghost_bc(2)
    upper_eqn(2) = 0.0d0
  end subroutine ys_compact_physical_boundary_equations

  !$omp declare target(ys_eliminate_physical_boundary_row)
  subroutine ys_eliminate_physical_boundary_row(local_idx, active_n, has_lower_boundary, has_upper_boundary, rhs_value, row_coeffs, &
                                                lower_ghost_value, lower_rhs0, lower_ghost_row, lower_eq0, upper_ghost_value, &
                                                upper_rhsn, upper_ghost_row, upper_eqn)
    implicit none
    integer(C_INT), intent(in) :: local_idx, active_n
    logical, intent(in) :: has_lower_boundary, has_upper_boundary
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs_value
    real(C_DOUBLE), intent(inout) :: row_coeffs(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: lower_ghost_value, lower_rhs0, upper_ghost_value, upper_rhsn
    real(C_DOUBLE), intent(in) :: lower_ghost_row(-2:2), lower_eq0(-2:2), upper_ghost_row(-2:2), upper_eqn(-2:2)
    real(C_DOUBLE) :: fac

    if (has_lower_boundary .and. local_idx == 0) then
      fac = row_coeffs(-2)/lower_ghost_row(-2)
      rhs_value = rhs_value - lower_ghost_value*fac
      row_coeffs = row_coeffs - lower_ghost_row*fac
      row_coeffs(-2) = 0.0d0
      fac = row_coeffs(-1)/lower_eq0(-1)
      rhs_value = rhs_value - lower_rhs0*fac
      row_coeffs = row_coeffs - lower_eq0*fac
      row_coeffs(-1) = 0.0d0
    else if (has_lower_boundary .and. local_idx == 1) then
      fac = row_coeffs(-2)/lower_eq0(-1)
      rhs_value = rhs_value - lower_rhs0*fac
      row_coeffs(-2:1) = row_coeffs(-2:1) - lower_eq0(-1:2)*fac
      row_coeffs(-2) = 0.0d0
    end if

    if (has_upper_boundary .and. local_idx == active_n - 2) then
      fac = row_coeffs(2)/upper_eqn(1)
      rhs_value = rhs_value - upper_rhsn*fac
      row_coeffs(-1:2) = row_coeffs(-1:2) - upper_eqn(-2:1)*fac
      row_coeffs(2) = 0.0d0
    else if (has_upper_boundary .and. local_idx == active_n - 1) then
      fac = row_coeffs(2)/upper_ghost_row(2)
      rhs_value = rhs_value - upper_ghost_value*fac
      row_coeffs = row_coeffs - upper_ghost_row*fac
      row_coeffs(2) = 0.0d0
      fac = row_coeffs(1)/upper_eqn(1)
      rhs_value = rhs_value - upper_rhsn*fac
      row_coeffs = row_coeffs - upper_eqn*fac
      row_coeffs(1) = 0.0d0
    end if
  end subroutine ys_eliminate_physical_boundary_row

subroutine ys_solve_full_y_lines_packed(slab_values, der_values, k2_values, bc0_values, bcn_values, ny_value, nz_value, nx0_value, &
                                          ni_value, component_index, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                          lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, &
                                          lambda_coeff, diffusion_coeff, first_line, line_count, solve_label)
    implicit none
    integer(C_INT), intent(in) :: ny_value, nz_value, nx0_value
    integer(C_INT), intent(in) :: component_index, lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index
    integer(C_INT), intent(in) :: first_line, line_count
    real(C_DOUBLE), intent(in) :: ni_value, lambda_coeff, diffusion_coeff
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    real(C_DOUBLE), intent(in) :: der_values(1:ny_value - 1, 0:3, -2:2)
    real(C_DOUBLE), intent(in) :: k2_values(-nz_value:, nx0_value:)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: slab_values(:, :)
    complex(C_DOUBLE_COMPLEX), intent(in) :: bc0_values(-nz_value:, nx0_value:, :)
    complex(C_DOUBLE_COMPLEX), intent(in) :: bcn_values(-nz_value:, nx0_value:, :)
    character(len=*), intent(in) :: solve_label
    complex(C_DOUBLE_COMPLEX) :: rhs_value, lower_ghost_value, lower_boundary_value, upper_boundary_value, upper_ghost_value
    complex(C_DOUBLE_COMPLEX) :: lower_rhs0, upper_rhsn
    real(C_DOUBLE) :: lower_eq0(-2:2), upper_eqn(-2:2), row_coeffs(-2:2)
    real(C_DOUBLE) :: kk, c_m2, c_m1, c_0, c_p1, c_p2
    integer(C_INT) :: nlines_z, ilocal, iline, ix, iz, iy, p

    if (line_count <= 0) return

    nlines_z = 2*nz_value + 1
    call ys_prepare_gpsv_workspace(ny_value - 1, line_count)

    call roctxPush("compact full-y pack")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(slab_values, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, der_values, k2_values, ni_value, &
    !$omp& bc0_values, bcn_values, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, lambda_coeff, diffusion_coeff, component_index, &
    !$omp& lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, first_line, line_count, nlines_z, nz_value, nx0_value, ny_value) &
    !$omp private(ilocal, iy, iline, ix, iz, p, kk, c_m2, c_m1, c_0, c_p1, c_p2, rhs_value, lower_ghost_value, &
    !$omp& lower_boundary_value, upper_boundary_value, upper_ghost_value, lower_rhs0, upper_rhsn, lower_eq0, upper_eqn, row_coeffs)
    do ilocal = 1, line_count
      do iy = 1, ny_value - 1
        iline = first_line + ilocal - 1
        ix = (iline - 1)/nlines_z + nx0_value
        iz = mod(iline - 1, nlines_z) - nz_value

        lower_ghost_value = (0.0d0, 0.0d0)
        if (lower_ghost_rhs_index > 0_C_INT) lower_ghost_value = bc0_values(iz, ix, lower_ghost_rhs_index)
        if (lower_ghost_rhs_index < 0_C_INT) lower_ghost_value = bcn_values(iz, ix, -lower_ghost_rhs_index)
        lower_boundary_value = (0.0d0, 0.0d0)
        if (lower_rhs_index > 0_C_INT) lower_boundary_value = bc0_values(iz, ix, lower_rhs_index)
        if (lower_rhs_index < 0_C_INT) lower_boundary_value = bcn_values(iz, ix, -lower_rhs_index)
        upper_boundary_value = (0.0d0, 0.0d0)
        if (upper_rhs_index > 0_C_INT) upper_boundary_value = bc0_values(iz, ix, upper_rhs_index)
        if (upper_rhs_index < 0_C_INT) upper_boundary_value = bcn_values(iz, ix, -upper_rhs_index)
        upper_ghost_value = (0.0d0, 0.0d0)
        if (upper_ghost_rhs_index > 0_C_INT) upper_ghost_value = bc0_values(iz, ix, upper_ghost_rhs_index)
        if (upper_ghost_rhs_index < 0_C_INT) upper_ghost_value = bcn_values(iz, ix, -upper_ghost_rhs_index)

        lower_rhs0 = lower_boundary_value - lower_ghost_value*lower_bc(-2)/lower_ghost_bc(-2)
        lower_eq0 = lower_bc - lower_ghost_bc*lower_bc(-2)/lower_ghost_bc(-2)
        lower_eq0(-2) = 0.0d0
        upper_rhsn = upper_boundary_value - upper_ghost_value*upper_bc(2)/upper_ghost_bc(2)
        upper_eqn = upper_bc - upper_ghost_bc*upper_bc(2)/upper_ghost_bc(2)
        upper_eqn(2) = 0.0d0

        kk = k2_values(iz, ix)
        rhs_value = slab_values(iy + 2, ilocal)
        if (component_index == 2_C_INT) then
          c_m2 = lambda_coeff*(der_values(iy, 2, -2) - kk*der_values(iy, 0, -2)) - &
                 ni_value*(der_values(iy, 3, -2) - 2.0d0*kk*der_values(iy, 2, -2) + kk*kk*der_values(iy, 0, -2))
          c_m1 = lambda_coeff*(der_values(iy, 2, -1) - kk*der_values(iy, 0, -1)) - &
                 ni_value*(der_values(iy, 3, -1) - 2.0d0*kk*der_values(iy, 2, -1) + kk*kk*der_values(iy, 0, -1))
          c_0 = lambda_coeff*(der_values(iy, 2, 0) - kk*der_values(iy, 0, 0)) - &
                ni_value*(der_values(iy, 3, 0) - 2.0d0*kk*der_values(iy, 2, 0) + kk*kk*der_values(iy, 0, 0))
          c_p1 = lambda_coeff*(der_values(iy, 2, 1) - kk*der_values(iy, 0, 1)) - &
                 ni_value*(der_values(iy, 3, 1) - 2.0d0*kk*der_values(iy, 2, 1) + kk*kk*der_values(iy, 0, 1))
          c_p2 = lambda_coeff*(der_values(iy, 2, 2) - kk*der_values(iy, 0, 2)) - &
                 ni_value*(der_values(iy, 3, 2) - 2.0d0*kk*der_values(iy, 2, 2) + kk*kk*der_values(iy, 0, 2))
        else
          c_m2 = lambda_coeff*der_values(iy, 0, -2) - diffusion_coeff*ni_value*(der_values(iy, 2, -2) - kk*der_values(iy, 0, -2))
          c_m1 = lambda_coeff*der_values(iy, 0, -1) - diffusion_coeff*ni_value*(der_values(iy, 2, -1) - kk*der_values(iy, 0, -1))
          c_0 = lambda_coeff*der_values(iy, 0, 0) - diffusion_coeff*ni_value*(der_values(iy, 2, 0) - kk*der_values(iy, 0, 0))
          c_p1 = lambda_coeff*der_values(iy, 0, 1) - diffusion_coeff*ni_value*(der_values(iy, 2, 1) - kk*der_values(iy, 0, 1))
          c_p2 = lambda_coeff*der_values(iy, 0, 2) - diffusion_coeff*ni_value*(der_values(iy, 2, 2) - kk*der_values(iy, 0, 2))
        end if

        row_coeffs = (/c_m2, c_m1, c_0, c_p1, c_p2/)
        call ys_eliminate_physical_boundary_row(iy - 1, ny_value - 1, .true., .true., rhs_value, row_coeffs, lower_ghost_value, &
                                                lower_rhs0, lower_ghost_bc, lower_eq0, upper_ghost_value, upper_rhsn, &
                                                upper_ghost_bc, upper_eqn)
        c_m2 = row_coeffs(-2)
        c_m1 = row_coeffs(-1)
        c_0 = row_coeffs(0)
        c_p1 = row_coeffs(1)
        c_p2 = row_coeffs(2)

        p = (iy - 1)*line_count + ilocal
        ys_gpsv_x(p) = rhs_value
        ys_gpsv_ds(p) = cmplx(c_m2, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_dl(p) = cmplx(c_m1, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_d(p) = cmplx(c_0, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_du(p) = cmplx(c_p1, 0.0d0, kind=C_DOUBLE)
        ys_gpsv_dw(p) = cmplx(c_p2, 0.0d0, kind=C_DOUBLE)
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("compact full-y pack")

    call ys_solve_packed_pentadiagonal(ny_value - 1, line_count, solve_label)

    call roctxPush("compact full-y unpack")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(slab_values, ys_gpsv_x, bc0_values, bcn_values, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
    !$omp& lower_rhs_index, lower_ghost_rhs_index, upper_rhs_index, upper_ghost_rhs_index, first_line, line_count, nlines_z, nz_value, nx0_value, ny_value) &
    !$omp private(ilocal, iy, iline, ix, iz, p, lower_ghost_value, lower_boundary_value, upper_boundary_value, upper_ghost_value, &
    !$omp& lower_rhs0, upper_rhsn, lower_eq0, upper_eqn)
    do ilocal = 1, line_count
      iline = first_line + ilocal - 1
      ix = (iline - 1)/nlines_z + nx0_value
      iz = mod(iline - 1, nlines_z) - nz_value
      do iy = 1, ny_value - 1
        p = (iy - 1)*line_count + ilocal
        slab_values(iy + 2, ilocal) = ys_gpsv_x(p)
      end do

      lower_ghost_value = (0.0d0, 0.0d0)
      if (lower_ghost_rhs_index > 0_C_INT) lower_ghost_value = bc0_values(iz, ix, lower_ghost_rhs_index)
      if (lower_ghost_rhs_index < 0_C_INT) lower_ghost_value = bcn_values(iz, ix, -lower_ghost_rhs_index)
      lower_boundary_value = (0.0d0, 0.0d0)
      if (lower_rhs_index > 0_C_INT) lower_boundary_value = bc0_values(iz, ix, lower_rhs_index)
      if (lower_rhs_index < 0_C_INT) lower_boundary_value = bcn_values(iz, ix, -lower_rhs_index)
      upper_boundary_value = (0.0d0, 0.0d0)
      if (upper_rhs_index > 0_C_INT) upper_boundary_value = bc0_values(iz, ix, upper_rhs_index)
      if (upper_rhs_index < 0_C_INT) upper_boundary_value = bcn_values(iz, ix, -upper_rhs_index)
      upper_ghost_value = (0.0d0, 0.0d0)
      if (upper_ghost_rhs_index > 0_C_INT) upper_ghost_value = bc0_values(iz, ix, upper_ghost_rhs_index)
      if (upper_ghost_rhs_index < 0_C_INT) upper_ghost_value = bcn_values(iz, ix, -upper_ghost_rhs_index)
      call ys_compact_physical_boundary_equations(lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, lower_boundary_value, &
                                                lower_ghost_value, upper_boundary_value, upper_ghost_value, lower_rhs0, lower_eq0, &
                                                  upper_rhsn, upper_eqn)

      slab_values(2, ilocal) = (lower_rhs0 - lower_eq0(0)*slab_values(3, ilocal) - lower_eq0(1)*slab_values(4, ilocal) - &
                                lower_eq0(2)*slab_values(5, ilocal))/lower_eq0(-1)
      slab_values(1, ilocal) = (lower_ghost_value - lower_ghost_bc(-1)*slab_values(2, ilocal) - &
                                lower_ghost_bc(0)*slab_values(3, ilocal) - lower_ghost_bc(1)*slab_values(4, ilocal) - &
                                lower_ghost_bc(2)*slab_values(5, ilocal))/lower_ghost_bc(-2)
      slab_values(ny_value + 2, ilocal) = (upper_rhsn - upper_eqn(-2)*slab_values(ny_value - 1, ilocal) - &
                                           upper_eqn(-1)*slab_values(ny_value, ilocal) - &
                                           upper_eqn(0)*slab_values(ny_value + 1, ilocal))/upper_eqn(1)
      slab_values(ny_value + 3, ilocal) = (upper_ghost_value - upper_ghost_bc(-2)*slab_values(ny_value - 1, ilocal) - &
                                           upper_ghost_bc(-1)*slab_values(ny_value, ilocal) - &
                                           upper_ghost_bc(0)*slab_values(ny_value + 1, ilocal) - &
                                           upper_ghost_bc(1)*slab_values(ny_value + 2, ilocal))/upper_ghost_bc(2)
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("compact full-y unpack")
  end subroutine ys_solve_full_y_lines_packed

  subroutine ys_lu5decomp(a)
    real(C_DOUBLE), intent(inout) :: a(0:, -2:)
    integer(C_INT) :: hi1, hi2
    real(C_DOUBLE) :: piv
    integer :: i, k, j

    hi1 = size(a, 1) - 1
    hi2 = size(a, 2) - 3
    a(hi1 - 2, 1:2) = 0
    a(hi1 - 3, 2) = 0
    do i = hi1 - hi2, 0, -1
      do k = hi2, 1, -1
        piv = a(i, k)
        do j = -1, -2, -1
          a(i, j + k) = a(i, j + k) - piv*a(i + k, j)
        end do
      end do
      piv = 1.0d0/a(i, 0)
      a(i, 0) = piv
      a(i, -2:-1) = a(i, -2:-1)*piv
    end do
    a(0, -2:-1) = 0
    a(1, -2) = 0
  end subroutine ys_lu5decomp

  subroutine ys_leftlu5div(x, a)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(-2:)
    real(C_DOUBLE), intent(in) :: a(0:, -2:)
    integer(C_INT) :: hi1, hi2, i

    hi1 = size(a, 1) - 1
    hi2 = size(a, 2) - 3

    do i = hi1 - hi2, 0, -1
      x(i) = x(i) - (a(i, 1)*x(i + 1) + a(i, 2)*x(i + 2))
      x(i) = x(i)*a(i, 0)
    end do

    do i = 0, hi1
      x(i) = x(i) - (a(i, -2)*x(i - 2) + a(i, -1)*x(i - 1))
    end do
  end subroutine ys_leftlu5div

  subroutine ys_solve_compact_derivative(f0, f1, der, d0mat, d140, d14m1, d14n, d14np1, ny, ny0, nyN)
    integer(C_INT), intent(in) :: ny, ny0, nyN
    complex(C_DOUBLE_COMPLEX), intent(in) :: f0(-1:ny + 1)
    complex(C_DOUBLE_COMPLEX), intent(out) :: f1(-1:ny + 1)
    real(C_DOUBLE), intent(in) :: der(ny0:nyN, 0:3, -2:2)
    real(C_DOUBLE), intent(in) :: d0mat(ny0:nyN + 2, -2:2)
    real(C_DOUBLE), intent(in) :: d140(-2:2), d14m1(-2:2), d14n(-2:2), d14np1(-2:2)
    integer(C_INT) :: iy

    f1(0) = sum(d140(-2:2)*f0(-1:3))
    f1(-1) = sum(d14m1(-2:2)*f0(-1:3))
    f1(ny) = sum(d14n(-2:2)*f0(ny - 3:ny + 1))
    f1(ny + 1) = sum(d14np1(-2:2)*f0(ny - 3:ny + 1))
    do iy = ny0, nyN
      f1(iy) = sum(der(iy, 1, -2:2)*f0(iy - 2:iy + 2))
    end do
    f1(1) = f1(1) - (der(1, 0, -1)*f1(0) + der(1, 0, -2)*f1(-1))
    f1(2) = f1(2) - der(2, 0, -2)*f1(0)
    f1(ny - 1) = f1(ny - 1) - (der(ny - 1, 0, 1)*f1(ny) + der(ny - 1, 0, 2)*f1(ny + 1))
    f1(ny - 2) = f1(ny - 2) - der(ny - 2, 0, 2)*f1(ny)
    call ys_leftlu5div(f1, d0mat)
  end subroutine ys_solve_compact_derivative
  subroutine ys_solve_compact_system(x, a, lower_bc, lower_ghost_bc, upper_bc, upper_ghost_bc, &
                                     rhs_lower, rhs_lower_ghost, rhs_upper, rhs_upper_ghost, ny, ny0, nyN)
    integer(C_INT), intent(in) :: ny, ny0, nyN
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(-1:ny + 1)
    real(C_DOUBLE), intent(inout) :: a(ny0:nyN + 2, -2:2)
    real(C_DOUBLE), intent(in) :: lower_bc(-2:2), lower_ghost_bc(-2:2), upper_bc(-2:2), upper_ghost_bc(-2:2)
    complex(C_DOUBLE_COMPLEX), intent(in) :: rhs_lower, rhs_lower_ghost, rhs_upper, rhs_upper_ghost
    real(C_DOUBLE) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)

    eqm1 = lower_ghost_bc
    eq0 = lower_bc
    eqn = upper_bc
    eqnp1 = upper_ghost_bc
    x(-1) = rhs_lower_ghost
    x(0) = rhs_lower
    x(ny) = rhs_upper
    x(ny + 1) = rhs_upper_ghost
    call ys_solve_ghost_system(x, a, eqm1, eq0, eqn, eqnp1, ny)
  end subroutine ys_solve_compact_system

  subroutine ys_solve_ghost_system(x, a, eqm1, eq0, eqn, eqnp1, ny)
    integer(C_INT), intent(in) :: ny
    complex(C_DOUBLE_COMPLEX), intent(inout) :: x(-1:ny + 1)
    real(C_DOUBLE), intent(inout) :: a(1:ny + 1, -2:2)
    real(C_DOUBLE), intent(inout) :: eqm1(-2:2), eq0(-2:2), eqn(-2:2), eqnp1(-2:2)

    x(0) = x(0) - x(-1)*eq0(-2)/eqm1(-2)
    eq0(-2:2) = eq0(-2:2) - eqm1(-2:2)*eq0(-2)/eqm1(-2)
    eq0(-2) = 0.0d0

    x(1) = x(1) - x(-1)*a(1, -2)/eqm1(-2)
    a(1, -2:2) = a(1, -2:2) - eqm1(-2:2)*a(1, -2)/eqm1(-2)
    a(1, -2) = 0.0d0

    x(1) = x(1) - x(0)*a(1, -1)/eq0(-1)
    a(1, -2:2) = a(1, -2:2) - eq0(-2:2)*a(1, -1)/eq0(-1)
    a(1, -1) = 0.0d0

    x(2) = x(2) - x(0)*a(2, -2)/eq0(-1)
    a(2, -2:1) = a(2, -2:1) - eq0(-1:2)*a(2, -2)/eq0(-1)
    a(2, -2) = 0.0d0

    x(ny) = x(ny) - x(ny + 1)*eqn(2)/eqnp1(2)
    eqn(-2:2) = eqn(-2:2) - eqnp1(-2:2)*eqn(2)/eqnp1(2)
    eqn(2) = 0.0d0

    x(ny - 1) = x(ny - 1) - x(ny + 1)*a(ny - 1, 2)/eqnp1(2)
    a(ny - 1, -2:2) = a(ny - 1, -2:2) - eqnp1(-2:2)*a(ny - 1, 2)/eqnp1(2)
    a(ny - 1, 2) = 0.0d0

    x(ny - 1) = x(ny - 1) - x(ny)*a(ny - 1, 1)/eqn(1)
    a(ny - 1, -2:2) = a(ny - 1, -2:2) - eqn(-2:2)*a(ny - 1, 1)/eqn(1)
    a(ny - 1, 1) = 0.0d0

    x(ny - 2) = x(ny - 2) - x(ny)*a(ny - 2, 2)/eqn(1)
    a(ny - 2, -1:2) = a(ny - 2, -1:2) - eqn(-2:1)*a(ny - 2, 2)/eqn(1)
    a(ny - 2, 2) = 0.0d0

    call ys_lu5decomp(a)
    call ys_leftlu5div(x, a)

    x(0) = (x(0) - sum(eq0(0:2)*x(1:3)))/eq0(-1)
    x(-1) = (x(-1) - sum(eqm1(-1:2)*x(0:3)))/eqm1(-2)
    x(ny) = (x(ny) - sum(eqn(-2:0)*x(ny - 3:ny - 1)))/eqn(1)
    x(ny + 1) = (x(ny + 1) - sum(eqnp1(-2:1)*x(ny - 3:ny)))/eqnp1(2)
  end subroutine ys_solve_ghost_system

  subroutine ys_solve_ghost_field_reduced(dst, ny, nz)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)

    call ys_solve_assembled_ghost_field(dst, ny, nz, YS_ENDPOINT_RESPONSE_FULL)
  end subroutine ys_solve_ghost_field_reduced

  subroutine ys_endpoint_context(dst, ny, nz, row_start, row_end, active_n, nlines, nlines_z, dst_row_base, &
                                 has_lower_boundary, has_upper_boundary, has_padded_dst)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(in) :: dst(:, :, :)
    integer(C_INT), intent(out) :: row_start, row_end, active_n, nlines, nlines_z, dst_row_base
    logical, intent(out) :: has_lower_boundary, has_upper_boundary, has_padded_dst

    row_start = ny0
    row_end = nyN
    active_n = row_end - row_start + 1
    nlines = ys_workspace_nlines
    nlines_z = 2*nz + 1
    has_lower_boundary = (row_start == 1)
    has_upper_boundary = (row_end == ny - 1)
    has_padded_dst = (size(dst, 1) == active_n + 4)
    dst_row_base = 1
    if (has_padded_dst) dst_row_base = 3
  end subroutine ys_endpoint_context

  subroutine ys_padded_inner_indices(active_n, dst_row_base, lower_inner0, lower_inner2, upper_inner0, upper_inner2, upper_inner3)
    implicit none
    integer(C_INT), intent(in) :: active_n, dst_row_base
    integer(C_INT), intent(out) :: lower_inner0, lower_inner2, upper_inner0, upper_inner2, upper_inner3

    lower_inner0 = dst_row_base
    lower_inner2 = dst_row_base + 2
    upper_inner0 = active_n + dst_row_base - 3
    upper_inner2 = upper_inner0 + 2
    upper_inner3 = upper_inner0 + 3
  end subroutine ys_padded_inner_indices

  subroutine ys_solve_assembled_ghost_field(dst, ny, nz, response_mode)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, response_mode
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)
    integer(C_INT) :: row_start, row_end, active_n, nlines, nlines_z, dst_row_base
    logical :: has_lower_boundary, has_upper_boundary, has_padded_dst

    call ys_endpoint_context(dst, ny, nz, row_start, row_end, active_n, nlines, nlines_z, dst_row_base, &
                             has_lower_boundary, has_upper_boundary, has_padded_dst)
    if (.not. has_padded_dst .and. size(dst, 1) /= active_n) then
      error stop "ys_solve_assembled_ghost_field expected either local-only or ghost-padded dst"
    end if
    if (.not. allocated(ys_local_rhs)) error stop "ys_prepare_ghost_field_workspace must be called before assembled ghost solve"
    if (ys_workspace_ny /= ny .or. ys_workspace_nz /= nz .or. ys_workspace_nx /= size(dst, 3) .or. &
        ys_workspace_active_n /= active_n) then
      error stop "ys_solve_assembled_ghost_field workspace does not match requested local solve dimensions"
    end if

    call ys_eliminate_assembled_boundaries(row_start, active_n, nlines, has_lower_boundary, has_upper_boundary)

    if (npy_grid == 1) then
      call ys_pack_assembled_direct(row_start, active_n, nlines)
      call ys_solve_ghost_field_single_rank_packed(dst, ny, nz)
      call ys_reconstruct_single_rank_boundaries(dst, ny, nz)
      return
    end if

    call roctxPush("ys_endpoint_schur")
    call ys_solve_endpoint_schur(dst, ny, nz, has_lower_boundary, has_upper_boundary, has_padded_dst, &
                                 row_start, row_end, active_n, nlines, nlines_z, dst_row_base, response_mode)
    call roctxPop("ys_endpoint_schur")
  end subroutine ys_solve_assembled_ghost_field

  subroutine ys_eliminate_assembled_boundaries(row_start, active_n, nlines, has_lower_boundary, has_upper_boundary)
    implicit none
    integer(C_INT), intent(in) :: row_start, active_n, nlines
    logical, intent(in) :: has_lower_boundary, has_upper_boundary
    complex(C_DOUBLE_COMPLEX) :: lower_rhs0, upper_rhsn, rhs_value
    real(C_DOUBLE) :: row_coeffs(-2:2), lower_eq0(-2:2), upper_eqn(-2:2)
    integer(C_INT) :: iline, row, local_idx

    if (.not. has_lower_boundary .and. .not. has_upper_boundary) return

    call roctxPush("ys_eliminate_assembled_boundaries")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
    !$omp& ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ys_boundary_lower_rhs0, ys_boundary_upper_rhsn, &
    !$omp& ys_boundary_lower_eq, ys_boundary_upper_eq, row_start, active_n, nlines, has_lower_boundary, has_upper_boundary) &
    !$omp private(iline, lower_rhs0, upper_rhsn, rhs_value, row_coeffs, lower_eq0, upper_eqn, row, local_idx)
    do iline = 1, nlines
      lower_rhs0 = (0.0d0, 0.0d0)
      lower_eq0 = 0.0d0
      if (has_lower_boundary) then
        lower_rhs0 = ys_lower_boundary_rhs(iline) - &
                     ys_lower_ghost_rhs(iline)*ys_lower_boundary_row(-2, iline)/ys_lower_ghost_row(-2, iline)
        lower_eq0 = ys_lower_boundary_row(:, iline) - &
                    ys_lower_ghost_row(:, iline)*ys_lower_boundary_row(-2, iline)/ys_lower_ghost_row(-2, iline)
        lower_eq0(-2) = 0.0d0
        ys_boundary_lower_rhs0(iline) = lower_rhs0
        ys_boundary_lower_eq(:, iline) = lower_eq0(-1:2)

        do local_idx = 0, 1
          row = row_start + local_idx
          rhs_value = ys_local_rhs(row, iline)
          row_coeffs = ys_local_operator(row, -2:2, iline)
          call ys_eliminate_physical_boundary_row(local_idx, active_n, .true., .false., rhs_value, row_coeffs, ys_lower_ghost_rhs(iline), &
                                                  lower_rhs0, ys_lower_ghost_row(:, iline), lower_eq0, (0.0d0, 0.0d0), &
                                                  (0.0d0, 0.0d0), ys_upper_ghost_row(:, iline), upper_eqn)
          ys_local_rhs(row, iline) = rhs_value
          ys_local_operator(row, -2:2, iline) = row_coeffs
        end do
      end if

      upper_rhsn = (0.0d0, 0.0d0)
      upper_eqn = 0.0d0
      if (has_upper_boundary) then
        upper_rhsn = ys_upper_boundary_rhs(iline) - &
                     ys_upper_ghost_rhs(iline)*ys_upper_boundary_row(2, iline)/ys_upper_ghost_row(2, iline)
        upper_eqn = ys_upper_boundary_row(:, iline) - &
                    ys_upper_ghost_row(:, iline)*ys_upper_boundary_row(2, iline)/ys_upper_ghost_row(2, iline)
        upper_eqn(2) = 0.0d0
        ys_boundary_upper_rhsn(iline) = upper_rhsn
        ys_boundary_upper_eq(:, iline) = upper_eqn(-2:1)

        do local_idx = active_n - 2, active_n - 1
          row = row_start + local_idx
          rhs_value = ys_local_rhs(row, iline)
          row_coeffs = ys_local_operator(row, -2:2, iline)
          call ys_eliminate_physical_boundary_row(local_idx, active_n, .false., .true., rhs_value, row_coeffs, (0.0d0, 0.0d0), &
                                                  (0.0d0, 0.0d0), ys_lower_ghost_row(:, iline), lower_eq0, ys_upper_ghost_rhs(iline), &
                                                  upper_rhsn, ys_upper_ghost_row(:, iline), upper_eqn)
          ys_local_rhs(row, iline) = rhs_value
          ys_local_operator(row, -2:2, iline) = row_coeffs
        end do
      end if
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_eliminate_assembled_boundaries")
  end subroutine ys_eliminate_assembled_boundaries

  subroutine ys_solve_ghost_field_reduced_const_operator(dst, ny, nz)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)

    call ys_solve_assembled_ghost_field(dst, ny, nz, YS_ENDPOINT_RESPONSE_CONST)
  end subroutine ys_solve_ghost_field_reduced_const_operator

  subroutine ys_solve_ghost_field_reduced_symmetric_operator(dst, ny, nz)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)

    call ys_solve_assembled_ghost_field(dst, ny, nz, YS_ENDPOINT_RESPONSE_SYMMETRIC)
  end subroutine ys_solve_ghost_field_reduced_symmetric_operator

  subroutine ys_solve_endpoint_schur(dst, ny, nz, has_lower_boundary, has_upper_boundary, has_padded_dst, &
                                     row_start, row_end, active_n, nlines, nlines_z, dst_row_base, response_mode)
    implicit none
    integer(C_INT), intent(in) :: ny, nz, row_start, row_end, active_n, nlines, nlines_z, dst_row_base, response_mode
    logical, intent(in) :: has_lower_boundary, has_upper_boundary, has_padded_dst
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)
    complex(C_DOUBLE_COMPLEX) :: lower_rhs0, upper_rhsn, rhs_value
    complex(C_DOUBLE_COMPLEX) :: c0, c1, y00, y01, y10, y11, g0, g1
    complex(C_DOUBLE_COMPLEX) :: s00, s01, s10, s11, t00, t01, t10, t11, det
    complex(C_DOUBLE_COMPLEX) :: vleft1, vleft2, vright1, vright2
    complex(C_DOUBLE_COMPLEX) :: s4(4, 4), rhs4(4, 5), pivot4, factor4
    real(C_DOUBLE) :: row_coeffs(-2:2), lower_eq0(-2:2), upper_eqn(-2:2), coeff
    integer(C_INT) :: nI, nresp, batch_count, exposed_n, interior_base
    integer(C_INT) :: sys, iline, ref_iline, resp, resp_index
    integer(C_INT) :: local_i, local_idx, row, col, coupled_row, p, j, offset
    integer(C_INT) :: exposed_slot, response_slot, rhs_col, iface, k, m
    integer(C_INT) :: ix, iz, abs_iz, ix_local
    integer(C_INT) :: lower_inner0, lower_inner2, upper_inner0, upper_inner2, upper_inner3
    logical :: is_actual
    logical :: has_left_interface, has_right_interface

    has_left_interface = (ipy > 0)
    has_right_interface = (ipy < npy_grid - 1)
    exposed_n = 0
    if (has_left_interface) exposed_n = exposed_n + 2
    if (has_right_interface) exposed_n = exposed_n + 2
    interior_base = 0
    if (has_left_interface) interior_base = 2
    nI = active_n - exposed_n
    if (nI <= 0) error stop "endpoint Schur solve needs at least one interior row"
    select case (response_mode)
    case (YS_ENDPOINT_RESPONSE_FULL)
      nresp = nlines
    case (YS_ENDPOINT_RESPONSE_CONST)
      nresp = 1_C_INT
    case (YS_ENDPOINT_RESPONSE_SYMMETRIC)
      nresp = (nxN - nx0 + 1)*(nz + 1)
    case default
      error stop "unknown endpoint Schur response mode"
    end select
    batch_count = nlines + exposed_n*nresp

    if (.not. allocated(ys_local_rhs)) error stop "ys_prepare_ghost_field_workspace must be called before endpoint Schur solve"
    call ys_prepare_gpsv_workspace(nI, batch_count)

    call roctxPush("ys_endpoint_pack_plus_response")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_local_rhs, ys_local_operator, ys_lower_ghost_rhs, ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_upper_ghost_rhs, &
    !$omp& ys_lower_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, ys_upper_ghost_row, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, &
    !$omp& ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, row_start, active_n, nI, nlines, nlines_z, nresp, batch_count, nz, nx0, &
    !$omp& has_lower_boundary, has_upper_boundary, response_mode, ipy, has_left_interface, has_right_interface, exposed_n, interior_base) &
    !$omp private(sys, iline, ref_iline, resp, local_i, local_idx, row, col, coupled_row, p, j, offset, is_actual, ix_local, abs_iz, &
    !$omp& rhs_value, row_coeffs, coeff, exposed_slot, response_slot)
    do local_i = 0, nI - 1
      do sys = 1, batch_count
        is_actual = (sys <= nlines)
        if (is_actual) then
          iline = sys
          ref_iline = iline
        else
          resp = mod(sys - nlines - 1, nresp) + 1
          select case (response_mode)
          case (YS_ENDPOINT_RESPONSE_FULL)
            ref_iline = resp
          case (YS_ENDPOINT_RESPONSE_CONST)
            ref_iline = 1_C_INT
          case default
            ix_local = (resp - 1)/(nz + 1)
            abs_iz = mod(resp - 1, nz + 1)
            ref_iline = ix_local*nlines_z + (abs_iz + nz + 1)
          end select
          iline = ref_iline
        end if
        local_idx = interior_base + local_i
        row = row_start + local_idx
        rhs_value = (0.0d0, 0.0d0)
        if (is_actual) rhs_value = ys_local_rhs(row, iline)
        row_coeffs = ys_local_operator(row, -2:2, ref_iline)

        p = local_i*batch_count + sys
        ys_gpsv_ds(p) = (0.0d0, 0.0d0)
        ys_gpsv_dl(p) = (0.0d0, 0.0d0)
        ys_gpsv_d(p) = (0.0d0, 0.0d0)
        ys_gpsv_du(p) = (0.0d0, 0.0d0)
        ys_gpsv_dw(p) = (0.0d0, 0.0d0)

        do col = -2, 2
          coeff = row_coeffs(col)
          if (coeff == 0.0d0) cycle
          coupled_row = local_idx + col
          if (coupled_row >= interior_base .and. coupled_row <= interior_base + nI - 1) then
            j = coupled_row - interior_base
            offset = j - local_i
          else
            exposed_slot = 0
            if (has_left_interface) then
              if (coupled_row == 0) exposed_slot = 1
              if (coupled_row == 1) exposed_slot = 2
            end if
            if (has_right_interface) then
              if (coupled_row == active_n - 2) exposed_slot = exposed_n - 1
              if (coupled_row == active_n - 1) exposed_slot = exposed_n
            end if
            if (exposed_slot > 0) then
              if (sys > nlines) then
                response_slot = (sys - nlines - 1)/nresp + 1
                if (response_slot == exposed_slot) rhs_value = cmplx(coeff, 0.0d0, kind=C_DOUBLE)
              end if
            end if
            cycle
          end if

          select case (offset)
          case (-2)
            ys_gpsv_ds(p) = cmplx(coeff, 0.0d0, kind=C_DOUBLE)
          case (-1)
            ys_gpsv_dl(p) = cmplx(coeff, 0.0d0, kind=C_DOUBLE)
          case (0)
            ys_gpsv_d(p) = cmplx(coeff, 0.0d0, kind=C_DOUBLE)
          case (1)
            ys_gpsv_du(p) = cmplx(coeff, 0.0d0, kind=C_DOUBLE)
          case (2)
            ys_gpsv_dw(p) = cmplx(coeff, 0.0d0, kind=C_DOUBLE)
          end select
        end do
        ys_gpsv_x(p) = rhs_value
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_endpoint_pack_plus_response")

    call ys_solve_packed_pentadiagonal(nI, batch_count, "ys_endpoint_gpsv_plus_response")

    call roctxPush("ys_endpoint_pack_schur")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_local_rhs, ys_local_operator, ys_gpsv_x, ys_reduced_rows_send, nI, nlines, nlines_z, nresp, batch_count, &
    !$omp& nz, nx0, row_start, active_n, response_mode, ipy, has_left_interface, has_right_interface, exposed_n, interior_base) &
    !$omp private(iline, row, ix, iz, abs_iz, resp_index, c0, c1, y00, y01, y10, y11, g0, g1, s00, s01, s10, s11, &
    !$omp& t00, t01, t10, t11, det, row_coeffs, s4, rhs4, pivot4, factor4, local_idx, col, coeff, coupled_row, j, &
    !$omp& exposed_slot, iface, k, m, rhs_col)
    do iline = 1, nlines
      ix = (iline - 1)/nlines_z + nx0
      iz = mod(iline - 1, nlines_z) - nz
      select case (response_mode)
      case (YS_ENDPOINT_RESPONSE_FULL)
        resp_index = iline
      case (YS_ENDPOINT_RESPONSE_CONST)
        resp_index = 1_C_INT
      case default
        abs_iz = abs(iz)
        resp_index = (ix - nx0)*(nz + 1) + abs_iz + 1
      end select
      ys_reduced_rows_send(:, iline) = (0.0d0, 0.0d0)
      if (.not. has_left_interface) then
        row = row_start + nI
        row_coeffs = ys_local_operator(row, -2:2, iline)
        c0 = ys_gpsv_x((nI - 2)*batch_count + iline)
        c1 = ys_gpsv_x((nI - 1)*batch_count + iline)
        y00 = ys_gpsv_x((nI - 2)*batch_count + nlines + resp_index)
        y10 = ys_gpsv_x((nI - 1)*batch_count + nlines + resp_index)
        y01 = ys_gpsv_x((nI - 2)*batch_count + nlines + nresp + resp_index)
        y11 = ys_gpsv_x((nI - 1)*batch_count + nlines + nresp + resp_index)
        g0 = ys_local_rhs(row, iline) - row_coeffs(-2)*c0 - row_coeffs(-1)*c1
        s00 = cmplx(row_coeffs(0), 0.0d0, kind=C_DOUBLE) - row_coeffs(-2)*y00 - row_coeffs(-1)*y10
        s01 = cmplx(row_coeffs(1), 0.0d0, kind=C_DOUBLE) - row_coeffs(-2)*y01 - row_coeffs(-1)*y11
        t00 = cmplx(row_coeffs(2), 0.0d0, kind=C_DOUBLE)
        t01 = (0.0d0, 0.0d0)

        row = row_start + nI + 1
        row_coeffs = ys_local_operator(row, -2:2, iline)
        c1 = ys_gpsv_x((nI - 1)*batch_count + iline)
        y10 = ys_gpsv_x((nI - 1)*batch_count + nlines + resp_index)
        y11 = ys_gpsv_x((nI - 1)*batch_count + nlines + nresp + resp_index)
        g1 = ys_local_rhs(row, iline) - row_coeffs(-2)*c1
        s10 = cmplx(row_coeffs(-1), 0.0d0, kind=C_DOUBLE) - row_coeffs(-2)*y10
        s11 = cmplx(row_coeffs(0), 0.0d0, kind=C_DOUBLE) - row_coeffs(-2)*y11
        t10 = cmplx(row_coeffs(1), 0.0d0, kind=C_DOUBLE)
        t11 = cmplx(row_coeffs(2), 0.0d0, kind=C_DOUBLE)

        det = s00*s11 - s01*s10
        ys_reduced_rows_send(3, iline) = (s11*g0 - s01*g1)/det
        ys_reduced_rows_send(4, iline) = (-s10*g0 + s00*g1)/det
        ys_reduced_rows_send(15, iline) = -((s11*t00 - s01*t10)/det)
        ys_reduced_rows_send(16, iline) = -((-s10*t00 + s00*t10)/det)
        ys_reduced_rows_send(19, iline) = -((s11*t01 - s01*t11)/det)
        ys_reduced_rows_send(20, iline) = -((-s10*t01 + s00*t11)/det)
      else if (.not. has_right_interface) then
        row = row_start
        row_coeffs = ys_local_operator(row, -2:2, iline)
        c0 = ys_gpsv_x(iline)
        y00 = ys_gpsv_x(nlines + resp_index)
        y01 = ys_gpsv_x(nlines + nresp + resp_index)
        g0 = ys_local_rhs(row, iline) - row_coeffs(2)*c0
        s00 = cmplx(row_coeffs(0), 0.0d0, kind=C_DOUBLE) - row_coeffs(2)*y00
        s01 = cmplx(row_coeffs(1), 0.0d0, kind=C_DOUBLE) - row_coeffs(2)*y01
        t00 = cmplx(row_coeffs(-2), 0.0d0, kind=C_DOUBLE)
        t01 = cmplx(row_coeffs(-1), 0.0d0, kind=C_DOUBLE)

        row = row_start + 1
        row_coeffs = ys_local_operator(row, -2:2, iline)
        c0 = ys_gpsv_x(iline)
        c1 = ys_gpsv_x(batch_count + iline)
        y00 = ys_gpsv_x(nlines + resp_index)
        y10 = ys_gpsv_x(batch_count + nlines + resp_index)
        y01 = ys_gpsv_x(nlines + nresp + resp_index)
        y11 = ys_gpsv_x(batch_count + nlines + nresp + resp_index)
        g1 = ys_local_rhs(row, iline) - row_coeffs(1)*c0 - row_coeffs(2)*c1
        s10 = cmplx(row_coeffs(-1), 0.0d0, kind=C_DOUBLE) - row_coeffs(1)*y00 - row_coeffs(2)*y10
        s11 = cmplx(row_coeffs(0), 0.0d0, kind=C_DOUBLE) - row_coeffs(1)*y01 - row_coeffs(2)*y11
        t10 = (0.0d0, 0.0d0)
        t11 = cmplx(row_coeffs(-2), 0.0d0, kind=C_DOUBLE)

        det = s00*s11 - s01*s10
        ys_reduced_rows_send(1, iline) = (s11*g0 - s01*g1)/det
        ys_reduced_rows_send(2, iline) = (-s10*g0 + s00*g1)/det
        ys_reduced_rows_send(5, iline) = -((s11*t00 - s01*t10)/det)
        ys_reduced_rows_send(6, iline) = -((-s10*t00 + s00*t10)/det)
        ys_reduced_rows_send(9, iline) = -((s11*t01 - s01*t11)/det)
        ys_reduced_rows_send(10, iline) = -((-s10*t01 + s00*t11)/det)
      else
        do k = 1, 4
          do m = 1, 4
            s4(k, m) = (0.0d0, 0.0d0)
          end do
          do rhs_col = 1, 5
            rhs4(k, rhs_col) = (0.0d0, 0.0d0)
          end do
        end do
        do iface = 1, 4
          select case (iface)
          case (1)
            local_idx = 0
          case (2)
            local_idx = 1
          case (3)
            local_idx = active_n - 2
          case default
            local_idx = active_n - 1
          end select
          row = row_start + local_idx
          row_coeffs = ys_local_operator(row, -2:2, iline)
          rhs4(iface, 1) = ys_local_rhs(row, iline)
          do col = -2, 2
            coeff = row_coeffs(col)
            if (coeff == 0.0d0) cycle
            coupled_row = local_idx + col
            if (coupled_row >= interior_base .and. coupled_row <= interior_base + nI - 1) then
              j = coupled_row - interior_base
              rhs4(iface, 1) = rhs4(iface, 1) - coeff*ys_gpsv_x(j*batch_count + iline)
              do exposed_slot = 1, exposed_n
                s4(iface, exposed_slot) = s4(iface, exposed_slot) - &
                                          coeff*ys_gpsv_x(j*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)
              end do
            else if (coupled_row == -2) then
              rhs4(iface, 2) = rhs4(iface, 2) + coeff
            else if (coupled_row == -1) then
              rhs4(iface, 3) = rhs4(iface, 3) + coeff
            else if (coupled_row == active_n) then
              rhs4(iface, 4) = rhs4(iface, 4) + coeff
            else if (coupled_row == active_n + 1) then
              rhs4(iface, 5) = rhs4(iface, 5) + coeff
            else
              exposed_slot = 0
              if (coupled_row == 0) exposed_slot = 1
              if (coupled_row == 1) exposed_slot = 2
              if (coupled_row == active_n - 2) exposed_slot = 3
              if (coupled_row == active_n - 1) exposed_slot = 4
              if (exposed_slot > 0) s4(iface, exposed_slot) = s4(iface, exposed_slot) + coeff
            end if
          end do
        end do

        do k = 1, 4
          pivot4 = s4(k, k)
          do m = k + 1, 4
            factor4 = s4(m, k)/pivot4
            s4(m, k) = factor4
            do j = k + 1, 4
              s4(m, j) = s4(m, j) - factor4*s4(k, j)
            end do
            do rhs_col = 1, 5
              rhs4(m, rhs_col) = rhs4(m, rhs_col) - factor4*rhs4(k, rhs_col)
            end do
          end do
        end do
        do rhs_col = 1, 5
          do k = 4, 1, -1
            do j = k + 1, 4
              rhs4(k, rhs_col) = rhs4(k, rhs_col) - s4(k, j)*rhs4(j, rhs_col)
            end do
            rhs4(k, rhs_col) = rhs4(k, rhs_col)/s4(k, k)
          end do
        end do

        do iface = 1, 4
          ys_reduced_rows_send(iface, iline) = rhs4(iface, 1)
          ys_reduced_rows_send(4 + iface, iline) = -rhs4(iface, 2)
          ys_reduced_rows_send(8 + iface, iline) = -rhs4(iface, 3)
          ys_reduced_rows_send(12 + iface, iline) = -rhs4(iface, 4)
          ys_reduced_rows_send(16 + iface, iline) = -rhs4(iface, 5)
        end do
      end if
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_endpoint_pack_schur")

    call ys_solve_reduced_interfaces()

    call ys_padded_inner_indices(active_n, dst_row_base, lower_inner0, lower_inner2, upper_inner0, upper_inner2, upper_inner3)

    call roctxPush("ys_endpoint_reconstruct")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(dst, ys_gpsv_x, ys_reduced_rhs, ys_left_interface_values, ys_right_interface_values, ys_lower_ghost_rhs, ys_upper_ghost_rhs, &
    !$omp& ys_lower_boundary_rhs, ys_upper_boundary_rhs, ys_lower_ghost_row, ys_upper_ghost_row, ys_lower_boundary_row, ys_upper_boundary_row, &
    !$omp& nlines, nlines_z, nresp, nx0, nz, active_n, dst_row_base, has_lower_boundary, has_upper_boundary, lower_inner0, lower_inner2, &
    !$omp& upper_inner0, upper_inner2, upper_inner3, has_padded_dst, nI, batch_count, response_mode, ipy, has_left_interface, has_right_interface, interior_base) &
    !$omp private(iline, local_i, local_idx, p, ix, iz, abs_iz, resp_index, vleft1, vleft2, vright1, vright2, lower_rhs0, upper_rhsn, lower_eq0, upper_eqn, exposed_slot)
    do iline = 1, nlines
      ix = (iline - 1)/nlines_z + nx0
      iz = mod(iline - 1, nlines_z) - nz
      select case (response_mode)
      case (YS_ENDPOINT_RESPONSE_FULL)
        resp_index = iline
      case (YS_ENDPOINT_RESPONSE_CONST)
        resp_index = 1_C_INT
      case default
        abs_iz = abs(iz)
        resp_index = (ix - nx0)*(nz + 1) + abs_iz + 1
      end select
      vleft1 = (0.0d0, 0.0d0)
      vleft2 = (0.0d0, 0.0d0)
      vright1 = (0.0d0, 0.0d0)
      vright2 = (0.0d0, 0.0d0)
      if (has_left_interface) then
        vleft1 = ys_reduced_rhs(4*ipy + 1, iline)
        vleft2 = ys_reduced_rhs(4*ipy + 2, iline)
        dst(dst_row_base, iz + nz + 1, ix - nx0 + 1) = vleft1
        dst(dst_row_base + 1, iz + nz + 1, ix - nx0 + 1) = vleft2
      end if
      if (has_right_interface) then
        vright1 = ys_reduced_rhs(4*ipy + 3, iline)
        vright2 = ys_reduced_rhs(4*ipy + 4, iline)
        dst(active_n - 2 + dst_row_base, iz + nz + 1, ix - nx0 + 1) = vright1
        dst(active_n - 1 + dst_row_base, iz + nz + 1, ix - nx0 + 1) = vright2
      end if
      do local_i = 0, nI - 1
        p = local_i*batch_count + iline
        local_idx = interior_base + local_i
        dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = ys_gpsv_x(p)
        exposed_slot = 0
        if (has_left_interface) then
          exposed_slot = exposed_slot + 1
          dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = &
            dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) - &
            ys_gpsv_x(local_i*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)*vleft1
          exposed_slot = exposed_slot + 1
          dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = &
            dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) - &
            ys_gpsv_x(local_i*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)*vleft2
        end if
        if (has_right_interface) then
          exposed_slot = exposed_slot + 1
          dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = &
            dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) - &
            ys_gpsv_x(local_i*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)*vright1
          exposed_slot = exposed_slot + 1
          dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) = &
            dst(local_idx + dst_row_base, iz + nz + 1, ix - nx0 + 1) - &
            ys_gpsv_x(local_i*batch_count + nlines + (exposed_slot - 1)*nresp + resp_index)*vright2
        end if
      end do
      if (has_padded_dst) then
        if (has_lower_boundary) then
          lower_rhs0 = ys_lower_boundary_rhs(iline) - &
                       ys_lower_ghost_rhs(iline)*ys_lower_boundary_row(-2, iline)/ys_lower_ghost_row(-2, iline)
          lower_eq0 = ys_lower_boundary_row(:, iline) - &
                      ys_lower_ghost_row(:, iline)*ys_lower_boundary_row(-2, iline)/ys_lower_ghost_row(-2, iline)
          lower_eq0(-2) = 0.0d0
          dst(2, iz + nz + 1, ix - nx0 + 1) = (lower_rhs0 - &
                                               lower_eq0(0)*dst(lower_inner0, iz + nz + 1, ix - nx0 + 1) - &
                                               lower_eq0(1)*dst(lower_inner0 + 1, iz + nz + 1, ix - nx0 + 1) - &
                                               lower_eq0(2)*dst(lower_inner2, iz + nz + 1, ix - nx0 + 1))/lower_eq0(-1)
          dst(1, iz + nz + 1, ix - nx0 + 1) = (ys_lower_ghost_rhs(iline) - &
                                               ys_lower_ghost_row(-1, iline)*dst(2, iz + nz + 1, ix - nx0 + 1) - &
                                               ys_lower_ghost_row(0, iline)*dst(3, iz + nz + 1, ix - nx0 + 1) - &
                                               ys_lower_ghost_row(1, iline)*dst(4, iz + nz + 1, ix - nx0 + 1) - &
                                       ys_lower_ghost_row(2, iline)*dst(5, iz + nz + 1, ix - nx0 + 1))/ys_lower_ghost_row(-2, iline)
        else
          dst(1, iz + nz + 1, ix - nx0 + 1) = ys_left_interface_values(1, iline)
          dst(2, iz + nz + 1, ix - nx0 + 1) = ys_left_interface_values(2, iline)
        end if
        if (has_upper_boundary) then
          upper_rhsn = ys_upper_boundary_rhs(iline) - &
                       ys_upper_ghost_rhs(iline)*ys_upper_boundary_row(2, iline)/ys_upper_ghost_row(2, iline)
          upper_eqn = ys_upper_boundary_row(:, iline) - &
                      ys_upper_ghost_row(:, iline)*ys_upper_boundary_row(2, iline)/ys_upper_ghost_row(2, iline)
          upper_eqn(2) = 0.0d0
          dst(active_n + 3, iz + nz + 1, ix - nx0 + 1) = (upper_rhsn - &
                                                          upper_eqn(-2)*dst(upper_inner0, iz + nz + 1, ix - nx0 + 1) - &
                                                          upper_eqn(-1)*dst(upper_inner0 + 1, iz + nz + 1, ix - nx0 + 1) - &
                                                          upper_eqn(0)*dst(upper_inner2, iz + nz + 1, ix - nx0 + 1))/upper_eqn(1)
          dst(active_n + 4, iz + nz + 1, ix - nx0 + 1) = (ys_upper_ghost_rhs(iline) - &
                                                      ys_upper_ghost_row(-2, iline)*dst(upper_inner0, iz + nz + 1, ix - nx0 + 1) - &
                                                  ys_upper_ghost_row(-1, iline)*dst(upper_inner0 + 1, iz + nz + 1, ix - nx0 + 1) - &
                                                       ys_upper_ghost_row(0, iline)*dst(upper_inner2, iz + nz + 1, ix - nx0 + 1) - &
                             ys_upper_ghost_row(1, iline)*dst(upper_inner3, iz + nz + 1, ix - nx0 + 1))/ys_upper_ghost_row(2, iline)
        else
          dst(active_n + 3, iz + nz + 1, ix - nx0 + 1) = ys_right_interface_values(1, iline)
          dst(active_n + 4, iz + nz + 1, ix - nx0 + 1) = ys_right_interface_values(2, iline)
        end if
      end if
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_endpoint_reconstruct")
  end subroutine ys_solve_endpoint_schur

  subroutine ys_pack_assembled_direct(row_start, active_n, nlines)
    implicit none
    integer(C_INT), intent(in) :: row_start, active_n, nlines
    integer(C_INT) :: iline, p, row, local_idx

    call roctxPush("ys_direct_pack")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_local_rhs, ys_local_operator, ys_gpsv_ds, ys_gpsv_dl, ys_gpsv_d, ys_gpsv_du, ys_gpsv_dw, ys_gpsv_x, row_start, active_n, nlines) &
    !$omp private(iline, p, row, local_idx)
    do iline = 1, nlines
      do local_idx = 0, active_n - 1
        row = row_start + local_idx
        p = local_idx*nlines + iline
        ys_gpsv_x(p) = ys_local_rhs(row, iline)
        ys_gpsv_ds(p) = cmplx(ys_local_operator(row, -2, iline), 0.0d0, kind=C_DOUBLE)
        ys_gpsv_dl(p) = cmplx(ys_local_operator(row, -1, iline), 0.0d0, kind=C_DOUBLE)
        ys_gpsv_d(p) = cmplx(ys_local_operator(row, 0, iline), 0.0d0, kind=C_DOUBLE)
        ys_gpsv_du(p) = cmplx(ys_local_operator(row, 1, iline), 0.0d0, kind=C_DOUBLE)
        ys_gpsv_dw(p) = cmplx(ys_local_operator(row, 2, iline), 0.0d0, kind=C_DOUBLE)
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_direct_pack")
  end subroutine ys_pack_assembled_direct

  subroutine ys_reconstruct_single_rank_boundaries(dst, ny, nz)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)
    integer(C_INT) :: ix, iz, iline, nlines, nlines_z, active_n
    integer(C_INT) :: dst_row_base, lower_inner0, lower_inner2, upper_inner0, upper_inner2, upper_inner3
    logical :: has_padded_dst

    active_n = nyN - ny0 + 1
    has_padded_dst = (size(dst, 1) == active_n + 4)
    if (.not. has_padded_dst) return

    dst_row_base = 3
    nlines_z = 2*nz + 1
    nlines = ys_workspace_nlines
    lower_inner0 = dst_row_base
    lower_inner2 = dst_row_base + 2
    upper_inner0 = active_n + dst_row_base - 3
    upper_inner2 = upper_inner0 + 2
    upper_inner3 = upper_inner0 + 3

    call roctxPush("ys_single_rank_reconstruct_boundaries")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(dst, ys_lower_ghost_rhs, ys_upper_ghost_rhs, ys_lower_ghost_row, ys_upper_ghost_row, ys_boundary_lower_rhs0, &
    !$omp& ys_boundary_upper_rhsn, ys_boundary_lower_eq, ys_boundary_upper_eq, nlines, nlines_z, nx0, nz, active_n, lower_inner0, &
    !$omp& lower_inner2, upper_inner0, upper_inner2, upper_inner3) &
    !$omp private(iline, ix, iz)
    do iline = 1, nlines
      ix = (iline - 1)/nlines_z + nx0
      iz = mod(iline - 1, nlines_z) - nz

      dst(2, iz + nz + 1, ix - nx0 + 1) = (ys_boundary_lower_rhs0(iline) - &
                                           ys_boundary_lower_eq(0, iline)*dst(lower_inner0, iz + nz + 1, ix - nx0 + 1) - &
                                           ys_boundary_lower_eq(1, iline)*dst(lower_inner0 + 1, iz + nz + 1, ix - nx0 + 1) - &
                        ys_boundary_lower_eq(2, iline)*dst(lower_inner2, iz + nz + 1, ix - nx0 + 1))/ys_boundary_lower_eq(-1, iline)
      dst(1, iz + nz + 1, ix - nx0 + 1) = (ys_lower_ghost_rhs(iline) - &
                                           ys_lower_ghost_row(-1, iline)*dst(2, iz + nz + 1, ix - nx0 + 1) - &
                                           ys_lower_ghost_row(0, iline)*dst(3, iz + nz + 1, ix - nx0 + 1) - &
                                           ys_lower_ghost_row(1, iline)*dst(4, iz + nz + 1, ix - nx0 + 1) - &
                                       ys_lower_ghost_row(2, iline)*dst(5, iz + nz + 1, ix - nx0 + 1))/ys_lower_ghost_row(-2, iline)

      dst(active_n + 3, iz + nz + 1, ix - nx0 + 1) = (ys_boundary_upper_rhsn(iline) - &
                                                    ys_boundary_upper_eq(-2, iline)*dst(upper_inner0, iz + nz + 1, ix - nx0 + 1) - &
                                                ys_boundary_upper_eq(-1, iline)*dst(upper_inner0 + 1, iz + nz + 1, ix - nx0 + 1) - &
                         ys_boundary_upper_eq(0, iline)*dst(upper_inner2, iz + nz + 1, ix - nx0 + 1))/ys_boundary_upper_eq(1, iline)
      dst(active_n + 4, iz + nz + 1, ix - nx0 + 1) = (ys_upper_ghost_rhs(iline) - &
                                                      ys_upper_ghost_row(-2, iline)*dst(upper_inner0, iz + nz + 1, ix - nx0 + 1) - &
                                                  ys_upper_ghost_row(-1, iline)*dst(upper_inner0 + 1, iz + nz + 1, ix - nx0 + 1) - &
                                                      ys_upper_ghost_row(0, iline)*dst(upper_inner2, iz + nz + 1, ix - nx0 + 1) - &
                             ys_upper_ghost_row(1, iline)*dst(upper_inner3, iz + nz + 1, ix - nx0 + 1))/ys_upper_ghost_row(2, iline)
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_single_rank_reconstruct_boundaries")
  end subroutine ys_reconstruct_single_rank_boundaries

  subroutine ys_solve_ghost_field_single_rank_packed(dst, ny, nz)
    implicit none
    integer(C_INT), intent(in) :: ny, nz
    complex(C_DOUBLE_COMPLEX), intent(inout) :: dst(:, :, :)
    integer(C_INT) :: ix, iz, iline, p, nlines, nlines_z, row_start, row_end, active_n, row
    integer(C_INT) :: dst_row_base
    logical :: has_padded_dst

    row_start = ny0
    row_end = nyN
    active_n = row_end - row_start + 1
    has_padded_dst = (size(dst, 1) == active_n + 4)
    dst_row_base = 1
    if (has_padded_dst) dst_row_base = 3
    nlines_z = 2*nz + 1
    nlines = ys_workspace_nlines
    call ys_solve_packed_pentadiagonal(active_n, nlines, "ys_single_rank_gpsv")

    call roctxPush("ys_single_rank_unpack")
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(dst, ys_gpsv_x, nlines, nlines_z, nx0, nz, active_n, dst_row_base) &
    !$omp private(iline, p, ix, iz, row)
    do iline = 1, nlines
      do row = 0, active_n - 1
        ix = (iline - 1)/nlines_z + nx0
        iz = mod(iline - 1, nlines_z) - nz
        p = row*nlines + iline
        dst(row + dst_row_base, iz + nz + 1, ix - nx0 + 1) = ys_gpsv_x(p)
      end do
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_single_rank_unpack")

  end subroutine ys_solve_ghost_field_single_rank_packed

  subroutine ys_solve_reduced_interfaces()
    integer(C_INT), parameter :: bw = 5
    complex(C_DOUBLE_COMPLEX) :: solve_piv, solve_factor, packed_remote(20)
    integer(C_INT) :: nlines, niface, iline, iblock, row0, i, j, t

    nlines = size(ys_reduced_rows_send, 2)
    niface = 4*npy_grid

#ifdef HAVE_MPI
    call roctxPush("MPI_Allgather reduced_y_interfaces")
#ifndef HAVE_HIP
    !$omp target data use_device_ptr(ys_reduced_rows_send, ys_reduced_rows_recv)
#endif
    call MPI_Allgather(ys_reduced_rows_send, 20*nlines, MPI_DOUBLE_COMPLEX, ys_reduced_rows_recv, 20*nlines, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr)
#ifndef HAVE_HIP
    !$omp end target data
#endif
    call roctxPop("MPI_Allgather reduced_y_interfaces")
#else
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(ys_reduced_rows_send, ys_reduced_rows_recv, nlines) private(iline, i)
    do iline = 1, nlines
      do i = 1, 20
        ys_reduced_rows_recv(i, iline, 1) = ys_reduced_rows_send(i, iline)
      end do
    end do
    !$omp end target teams distribute parallel do
#endif

    call roctxPush("ys_reduced_interfaces_solve")
    !$omp target teams distribute parallel do default(none) &
    !$omp shared(ys_reduced_rows_recv, ys_reduced_matrix_lu, ys_reduced_rhs, ys_left_interface_values, ys_right_interface_values, nlines, niface, npy_grid, ipy) &
    !$omp private(iline, iblock, row0, packed_remote, i, j, t, solve_piv, solve_factor)
    do iline = 1, nlines
      ys_reduced_matrix_lu(:, :, iline) = (0.0d0, 0.0d0)
      ys_reduced_rhs(:, iline) = (0.0d0, 0.0d0)
      do iblock = 0, npy_grid - 1
        row0 = 4*iblock
        packed_remote = ys_reduced_rows_recv(:, iline, iblock + 1)
        ys_reduced_matrix_lu(row0 + 1, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 2, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 3, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_matrix_lu(row0 + 4, bw + 1, iline) = (1.0d0, 0.0d0)
        ys_reduced_rhs(row0 + 1:row0 + 4, iline) = packed_remote(1:4)
        if (iblock > 0) then
          ys_reduced_matrix_lu(row0 + 1, bw - 1, iline) = -packed_remote(5)
          ys_reduced_matrix_lu(row0 + 2, bw - 2, iline) = -packed_remote(6)
          ys_reduced_matrix_lu(row0 + 3, bw - 3, iline) = -packed_remote(7)
          ys_reduced_matrix_lu(row0 + 4, bw - 4, iline) = -packed_remote(8)
          ys_reduced_matrix_lu(row0 + 1, bw, iline) = -packed_remote(9)
          ys_reduced_matrix_lu(row0 + 2, bw - 1, iline) = -packed_remote(10)
          ys_reduced_matrix_lu(row0 + 3, bw - 2, iline) = -packed_remote(11)
          ys_reduced_matrix_lu(row0 + 4, bw - 3, iline) = -packed_remote(12)
        end if
        if (iblock < npy_grid - 1) then
          ys_reduced_matrix_lu(row0 + 1, bw + 5, iline) = -packed_remote(13)
          ys_reduced_matrix_lu(row0 + 2, bw + 4, iline) = -packed_remote(14)
          ys_reduced_matrix_lu(row0 + 3, bw + 3, iline) = -packed_remote(15)
          ys_reduced_matrix_lu(row0 + 4, bw + 2, iline) = -packed_remote(16)
          ys_reduced_matrix_lu(row0 + 1, bw + 6, iline) = -packed_remote(17)
          ys_reduced_matrix_lu(row0 + 2, bw + 5, iline) = -packed_remote(18)
          ys_reduced_matrix_lu(row0 + 3, bw + 4, iline) = -packed_remote(19)
          ys_reduced_matrix_lu(row0 + 4, bw + 3, iline) = -packed_remote(20)
        end if
      end do

      call ys_factor_banded_complex(ys_reduced_matrix_lu(:, :, iline))
      call ys_solve_factored_banded_complex(ys_reduced_rhs(:, iline), ys_reduced_matrix_lu(:, :, iline))

      ys_left_interface_values(:, iline) = (0.0d0, 0.0d0)
      ys_right_interface_values(:, iline) = (0.0d0, 0.0d0)
      row0 = 4*ipy
      if (ipy > 0) then
        ys_left_interface_values(:, iline) = ys_reduced_rhs(row0 - 1:row0, iline)
      end if
      if (ipy < npy_grid - 1) then
        ys_right_interface_values(:, iline) = ys_reduced_rhs(row0 + 5:row0 + 6, iline)
      end if
    end do
    !$omp end target teams distribute parallel do
    call roctxPop("ys_reduced_interfaces_solve")
  end subroutine ys_solve_reduced_interfaces

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_factor_banded_complex)
#endif
  subroutine ys_factor_banded_complex(a)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: a(:, :)
    integer(C_INT), parameter :: bw = 5
    integer(C_INT) :: n, i, j, t
    complex(C_DOUBLE_COMPLEX) :: piv, factor

    n = size(a, 1)
    do i = 1, n
      piv = a(i, bw + 1)
      do j = 1, min(bw, n - i)
        factor = a(i + j, bw + 1 - j)/piv
        a(i + j, bw + 1 - j) = factor
        do t = 1, min(bw, n - i)
          if (t - j > bw) cycle
          a(i + j, bw + 1 + t - j) = a(i + j, bw + 1 + t - j) - factor*a(i, bw + 1 + t)
        end do
      end do
    end do
  end subroutine ys_factor_banded_complex

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_solve_factored_banded_complex)
#endif
  subroutine ys_solve_factored_banded_complex(rhs, a)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: a(:, :)
    integer(C_INT), parameter :: bw = 5
    integer(C_INT) :: n, i, j

    n = size(a, 1)
    do i = 1, n
      do j = max(1_C_INT, i - bw), i - 1
        rhs(i) = rhs(i) - a(i, bw + 1 + j - i)*rhs(j)
      end do
    end do

    do i = n, 1, -1
      do j = i + 1, min(n, i + bw)
        rhs(i) = rhs(i) - a(i, bw + 1 + j - i)*rhs(j)
      end do
      rhs(i) = rhs(i)/a(i, bw + 1)
    end do
  end subroutine ys_solve_factored_banded_complex

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_factor_penta)
#endif
  subroutine ys_factor_penta(a)
    real(C_DOUBLE), intent(inout) :: a(0:, -2:)
    integer(C_INT) :: n, i
    real(C_DOUBLE) :: piv, factor

    ! Generic pentadiagonal LU on the local reduced block. The wall/ghost rows
    ! are already eliminated outside this factorization, so this is a pure
    ! interior solve rather than the boundary-aware variant used elsewhere.
    n = size(a, 1)
    do i = 0, n - 1
      piv = a(i, 0)
      a(i, 0) = 1.0d0/piv

      if (i + 1 < n) then
        factor = a(i + 1, -1)*a(i, 0)
        a(i + 1, -1) = factor
        a(i + 1, 0) = a(i + 1, 0) - factor*a(i, 1)
        if (i + 2 < n) a(i + 1, 1) = a(i + 1, 1) - factor*a(i, 2)
      end if

      if (i + 2 < n) then
        factor = a(i + 2, -2)*a(i, 0)
        a(i + 2, -2) = factor
        a(i + 2, -1) = a(i + 2, -1) - factor*a(i, 1)
        a(i + 2, 0) = a(i + 2, 0) - factor*a(i, 2)
      end if
    end do
  end subroutine ys_factor_penta

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_solve_factored_penta_one)
#endif
  subroutine ys_solve_factored_penta_one(rhs, a)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(0:)
    real(C_DOUBLE), intent(in) :: a(0:, -2:)
    integer(C_INT) :: n, i

    n = size(a, 1)

    do i = 0, n - 1
      if (i >= 2) rhs(i) = rhs(i) - a(i, -2)*rhs(i - 2)
      if (i >= 1) rhs(i) = rhs(i) - a(i, -1)*rhs(i - 1)
    end do

    do i = n - 1, 0, -1
      if (i + 1 < n) rhs(i) = rhs(i) - a(i, 1)*rhs(i + 1)
      if (i + 2 < n) rhs(i) = rhs(i) - a(i, 2)*rhs(i + 2)
      rhs(i) = rhs(i)*a(i, 0)
    end do
  end subroutine ys_solve_factored_penta_one

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_factor_penta_interleaved)
#endif
  subroutine ys_factor_penta_interleaved(ds, dl, d, du, dw, stride, first, n)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i, p, p1, p2
    complex(C_DOUBLE_COMPLEX) :: factor

    do i = 0, n - 1
      p = first + i*stride
      d(p) = 1.0d0/d(p)

      if (i + 1 < n) then
        p1 = first + (i + 1)*stride
        factor = dl(p1)*d(p)
        dl(p1) = factor
        d(p1) = d(p1) - factor*du(p)
        if (i + 2 < n) du(p1) = du(p1) - factor*dw(p)
      end if

      if (i + 2 < n) then
        p2 = first + (i + 2)*stride
        factor = ds(p2)*d(p)
        ds(p2) = factor
        dl(p2) = dl(p2) - factor*du(p)
        d(p2) = d(p2) - factor*dw(p)
      end if
    end do
  end subroutine ys_factor_penta_interleaved

#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  !$omp declare target(ys_solve_factored_penta_interleaved)
#endif
  subroutine ys_solve_factored_penta_interleaved(rhs, ds, dl, d, du, dw, stride, first, n)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: rhs(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: ds(:), dl(:), d(:), du(:), dw(:)
    integer(C_INT), intent(in) :: stride, first, n
    integer(C_INT) :: i, p

    do i = 0, n - 1
      p = first + i*stride
      if (i >= 2) rhs(p) = rhs(p) - ds(p)*rhs(first + (i - 2)*stride)
      if (i >= 1) rhs(p) = rhs(p) - dl(p)*rhs(first + (i - 1)*stride)
    end do

    do i = n - 1, 0, -1
      p = first + i*stride
      if (i + 1 < n) rhs(p) = rhs(p) - du(p)*rhs(first + (i + 1)*stride)
      if (i + 2 < n) rhs(p) = rhs(p) - dw(p)*rhs(first + (i + 2)*stride)
      rhs(p) = rhs(p)*d(p)
    end do
  end subroutine ys_solve_factored_penta_interleaved
end module y_line_solvers
