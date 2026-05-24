#include "header.h"

module convvelo

  use, intrinsic :: iso_c_binding
  use dnsdata, only: V, nPhi, nz, ny, der, nxd, izd, factor, iproc, D0mat, d240, d24m1, d24n, d24np1, &
                     COMPLEXderiv, LeftLU5div
  use pressure_output, only: compute_poisson, compute_dpdy
  use mpi_transpose, only: ny0, nyN, nx0, nxN, nxB, nzB, nx, has_average, ierr, sendbuf, recvbuf, &
                           pack_zTOx, unpack_zTOx, pack_xTOz, unpack_xTOz, alltoall, nzd
  use ffts, only: IFT, RFT, HFT, FFT, VVdx, VVdz
#ifdef HAVE_MPI
  use mpi_f08
#endif

  implicit none

  private

  integer(C_INT), parameter, public :: n_convvelo_velocity_fields = 33
  integer(C_INT), parameter, public :: n_convvelo_scalar_fields = 10
  integer(C_INT), parameter :: i_u_cross_u = 1
  integer(C_INT), parameter :: i_u_cross_dyu = 2
  integer(C_INT), parameter :: i_u_cross_v = 3
  integer(C_INT), parameter :: i_u_cross_dyv = 4
  integer(C_INT), parameter :: i_u_cross_w = 5
  integer(C_INT), parameter :: i_u_cross_dyw = 6
  integer(C_INT), parameter :: i_u_cross_dyyu = 7
  integer(C_INT), parameter :: i_v_cross_u = 8
  integer(C_INT), parameter :: i_v_cross_dyu = 9
  integer(C_INT), parameter :: i_v_cross_v = 10
  integer(C_INT), parameter :: i_v_cross_dyv = 11
  integer(C_INT), parameter :: i_v_cross_w = 12
  integer(C_INT), parameter :: i_v_cross_dyw = 13
  integer(C_INT), parameter :: i_v_cross_dyyv = 14
  integer(C_INT), parameter :: i_w_cross_u = 15
  integer(C_INT), parameter :: i_w_cross_dyu = 16
  integer(C_INT), parameter :: i_w_cross_v = 17
  integer(C_INT), parameter :: i_w_cross_dyv = 18
  integer(C_INT), parameter :: i_w_cross_w = 19
  integer(C_INT), parameter :: i_w_cross_dyw = 20
  integer(C_INT), parameter :: i_w_cross_dyyw = 21
  integer(C_INT), parameter :: i_t_theta_theta = 1
  integer(C_INT), parameter :: i_t_theta_u = 2
  integer(C_INT), parameter :: i_t_theta_v = 3
  integer(C_INT), parameter :: i_t_theta_w = 4
  integer(C_INT), parameter :: i_t_theta_dyytheta = 7
  integer(C_INT), parameter :: i_t_theta_dytheta = 9
  integer(C_INT), parameter :: i_t_theta_dyv = 10

  logical, save :: convvelo_initialized = .false.
  integer(C_INT), save, public :: n_convvelo_fields = 0
  integer(C_INT64_T), save :: n_mean_samples = 0_C_INT64_T

  complex(C_DOUBLE_COMPLEX), allocatable, save, public :: convvelo_stats(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save, public :: convvelo_work(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save, public :: component_means(:, :)
  real(C_DOUBLE), allocatable, save :: convvelo_real0(:, :, :)
  real(C_DOUBLE), allocatable, save :: convvelo_real1(:, :, :)
  real(C_DOUBLE), allocatable, save :: convvelo_real_prod(:, :, :)
  integer(C_INT64_T), allocatable, save :: n_field_samples(:)

  public :: init_convvelo, reset_convvelo_stats, update_convvelo_component_means, free_convvelo
  public :: start_convvelo_field, finish_convvelo_field
  public :: acc_convvelo_stats

contains

  subroutine init_convvelo()
    implicit none

    if (convvelo_initialized) return

    n_convvelo_fields = n_convvelo_velocity_fields + nPhi*n_convvelo_scalar_fields

    allocate (convvelo_stats(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, n_convvelo_fields))
    allocate (convvelo_work(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    allocate (component_means(ny0 - 2:nyN + 2, 1:3 + nPhi))
    allocate (convvelo_real0(2*(nxd + 1), nzB, ny + 3))
    allocate (convvelo_real1(2*(nxd + 1), nzB, ny + 3))
    allocate (convvelo_real_prod(2*(nxd + 1), nzB, ny + 3))
    allocate (n_field_samples(n_convvelo_fields))

    convvelo_stats = (0.0d0, 0.0d0)
    convvelo_work = (0.0d0, 0.0d0)
    component_means = (0.0d0, 0.0d0)
    convvelo_real0 = 0.0d0
    convvelo_real1 = 0.0d0
    convvelo_real_prod = 0.0d0
    n_field_samples = 0_C_INT64_T

    !$omp target enter data map(to: convvelo_stats, convvelo_work, component_means, convvelo_real0, convvelo_real1, convvelo_real_prod)

    n_mean_samples = 0_C_INT64_T
    convvelo_initialized = .true.
  end subroutine init_convvelo

  subroutine reset_convvelo_stats()
    implicit none

    if (.not. convvelo_initialized) return

    convvelo_stats = (0.0d0, 0.0d0)
    convvelo_work = (0.0d0, 0.0d0)
    component_means = (0.0d0, 0.0d0)
    convvelo_real0 = 0.0d0
    convvelo_real1 = 0.0d0
    convvelo_real_prod = 0.0d0
    n_field_samples = 0_C_INT64_T

    !$omp target update to(convvelo_stats, convvelo_work, component_means, convvelo_real0, convvelo_real1, convvelo_real_prod)

    n_mean_samples = 0_C_INT64_T
  end subroutine reset_convvelo_stats

  subroutine update_convvelo_component_means()
    implicit none

    complex(C_DOUBLE_COMPLEX) :: snapshot(ny0 - 2:nyN + 2, 1:3 + nPhi)
    real(C_DOUBLE) :: old_weight, new_weight

    if (.not. convvelo_initialized) call init_convvelo()

    snapshot = (0.0d0, 0.0d0)
    if (has_average) then
      !$omp target update from(V(ny0 - 2:nyN + 2, 0, 0, 1:3 + nPhi))
      snapshot(:, :) = V(:, 0, 0, :)
    end if

#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, snapshot, size(snapshot), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    n_mean_samples = n_mean_samples + 1_C_INT64_T
    old_weight = dble(n_mean_samples - 1_C_INT64_T)/dble(n_mean_samples)
    new_weight = 1.0d0/dble(n_mean_samples)
    component_means = old_weight*component_means + new_weight*snapshot

    !$omp target update to(component_means)
  end subroutine update_convvelo_component_means

  subroutine acc_convvelo_stats()
    implicit none
    integer(C_INT) :: iPhi, scalar_offset

    if (.not. convvelo_initialized) call init_convvelo()
    !$omp target update to(V)

    call accumulate_cross_components(i_u_cross_u, 1, 1)
    call accumulate_cross_derivative(i_u_cross_dyu, 1, 1, 1)
    call accumulate_cross_components(i_u_cross_v, 1, 2)
    call accumulate_cross_derivative(i_u_cross_dyv, 1, 2, 1)
    call accumulate_cross_components(i_u_cross_w, 1, 3)
    call accumulate_cross_derivative(i_u_cross_dyw, 1, 3, 1)
    call accumulate_cross_derivative(i_u_cross_dyyu, 1, 1, 2)

    call accumulate_cross_components(i_v_cross_u, 2, 1)
    call accumulate_cross_derivative(i_v_cross_dyu, 2, 1, 1)
    call accumulate_cross_components(i_v_cross_v, 2, 2)
    call accumulate_cross_derivative(i_v_cross_dyv, 2, 2, 1)
    call accumulate_cross_components(i_v_cross_w, 2, 3)
    call accumulate_cross_derivative(i_v_cross_dyw, 2, 3, 1)
    call accumulate_cross_derivative(i_v_cross_dyyv, 2, 2, 2)

    call accumulate_cross_components(i_w_cross_u, 3, 1)
    call accumulate_cross_derivative(i_w_cross_dyu, 3, 1, 1)
    call accumulate_cross_components(i_w_cross_v, 3, 2)
    call accumulate_cross_derivative(i_w_cross_dyv, 3, 2, 1)
    call accumulate_cross_components(i_w_cross_w, 3, 3)
    call accumulate_cross_derivative(i_w_cross_dyw, 3, 3, 1)
    call accumulate_cross_derivative(i_w_cross_dyyw, 3, 3, 2)

    do iPhi = 1, nPhi
      scalar_offset = n_convvelo_velocity_fields + (iPhi - 1)*n_convvelo_scalar_fields
      call accumulate_cross_components(scalar_offset + i_t_theta_theta, 3 + iPhi, 3 + iPhi)
      call accumulate_cross_components(scalar_offset + i_t_theta_u, 3 + iPhi, 1)
      call accumulate_cross_components(scalar_offset + i_t_theta_v, 3 + iPhi, 2)
      call accumulate_cross_components(scalar_offset + i_t_theta_w, 3 + iPhi, 3)
      call accumulate_cross_product_field(scalar_offset + 5, 3 + iPhi, 3*iPhi + 4)
      call accumulate_cross_product_field(scalar_offset + 6, 3 + iPhi, 3*iPhi + 6)
      call accumulate_cross_derivative(scalar_offset + i_t_theta_dyytheta, 3 + iPhi, 3 + iPhi, 2)
      call accumulate_cross_product_derivative_field(scalar_offset + 8, 3 + iPhi, 3*iPhi + 5)
      call accumulate_cross_derivative(scalar_offset + i_t_theta_dytheta, 3 + iPhi, 3 + iPhi, 1)
      call accumulate_cross_derivative(scalar_offset + i_t_theta_dyv, 3 + iPhi, 2, 1)
    end do
    call accumulate_cross_product_field(25, 1, 1)
    call accumulate_cross_product_field(26, 1, 6)
    call accumulate_cross_product_field(27, 2, 4)
    call accumulate_cross_product_field(28, 2, 5)
    call accumulate_cross_product_field(29, 3, 6)
    call accumulate_cross_product_field(30, 3, 3)
    call accumulate_cross_product_derivative_field(31, 1, 4)
    call accumulate_cross_product_derivative_field(32, 2, 2)
    call accumulate_cross_product_derivative_field(33, 3, 5)

    call accumulate_cross_pressure(22, 1, .false.)
    call accumulate_cross_pressure(23, 2, .true.)
    call accumulate_cross_pressure(24, 3, .false.)
  end subroutine acc_convvelo_stats

  subroutine start_convvelo_field()
    implicit none

    if (.not. convvelo_initialized) call init_convvelo()

    convvelo_work = (0.0d0, 0.0d0)
    !$omp target update to(convvelo_work)
  end subroutine start_convvelo_field

  subroutine finish_convvelo_field(field_index)
    implicit none

    integer(C_INT), intent(in) :: field_index
    real(C_DOUBLE) :: old_weight, new_weight

    if (.not. convvelo_initialized) call init_convvelo()
    if (field_index < 1 .or. field_index > n_convvelo_fields) then
      error stop "finish_convvelo_field: field_index out of range"
    end if

    n_field_samples(field_index) = n_field_samples(field_index) + 1_C_INT64_T
    old_weight = dble(n_field_samples(field_index) - 1_C_INT64_T)/dble(n_field_samples(field_index))
    new_weight = 1.0d0/dble(n_field_samples(field_index))
    convvelo_stats(:, :, :, field_index) = old_weight*convvelo_stats(:, :, :, field_index) + new_weight*convvelo_work

    !$omp target update to(convvelo_stats(:, :, :, field_index))
  end subroutine finish_convvelo_field

  subroutine load_component_to_work(component_index)
    implicit none
    integer(C_INT), intent(in) :: component_index

    convvelo_work(:, :, :) = V(:, :, :, component_index)
  end subroutine load_component_to_work

  subroutine apply_dy_to_work(component_index)
    implicit none
    integer(C_INT), intent(in) :: component_index
    integer(C_INT) :: iz, ix

    convvelo_work = (0.0d0, 0.0d0)
    do ix = nx0, nxN
      do iz = -nz, nz
        call COMPLEXderiv(V(:, iz, ix, component_index), convvelo_work(:, iz, ix), der, D0mat)
      end do
    end do
  end subroutine apply_dy_to_work

  subroutine apply_dyy_to_work(component_index)
    implicit none
    integer(C_INT), intent(in) :: component_index
    integer(C_INT) :: iy, iz, ix

    convvelo_work = (0.0d0, 0.0d0)
    do ix = nx0, nxN
      do iz = -nz, nz
        convvelo_work(0, iz, ix) = sum(d240(-2:2)*V(-1:3, iz, ix, component_index))
        convvelo_work(-1, iz, ix) = sum(d24m1(-2:2)*V(-1:3, iz, ix, component_index))
        convvelo_work(ny, iz, ix) = sum(d24n(-2:2)*V(ny - 3:ny + 1, iz, ix, component_index))
        convvelo_work(ny + 1, iz, ix) = sum(d24np1(-2:2)*V(ny - 3:ny + 1, iz, ix, component_index))
        do iy = ny0, nyN
          convvelo_work(iy, iz, ix) = sum(der(iy, 2, -2:2)*V(iy - 2:iy + 2, iz, ix, component_index))
        end do
        convvelo_work(1, iz, ix) = convvelo_work(1, iz, ix) - ( &
                                   der(1, 0, -1)*convvelo_work(0, iz, ix) + &
                                   der(1, 0, -2)*convvelo_work(-1, iz, ix))
        convvelo_work(2, iz, ix) = convvelo_work(2, iz, ix) - der(2, 0, -2)*convvelo_work(0, iz, ix)
        convvelo_work(ny - 1, iz, ix) = convvelo_work(ny - 1, iz, ix) - ( &
                                        der(ny - 1, 0, 1)*convvelo_work(ny, iz, ix) + &
                                        der(ny - 1, 0, 2)*convvelo_work(ny + 1, iz, ix))
        convvelo_work(ny - 2, iz, ix) = convvelo_work(ny - 2, iz, ix) - &
                                        der(ny - 2, 0, 2)*convvelo_work(ny, iz, ix)
        call LeftLU5div(convvelo_work(:, iz, ix), D0mat, convvelo_work(:, iz, ix))
      end do
    end do
  end subroutine apply_dyy_to_work

  subroutine multiply_work_by_conjugate(component_index)
    implicit none
    integer(C_INT), intent(in) :: component_index
    integer(C_INT) :: iy, iz, ix

    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0 - 2, nyN + 2
          convvelo_work(iy, iz, ix) = conjg(V(iy, iz, ix, component_index))*convvelo_work(iy, iz, ix)
        end do
      end do
    end do
  end subroutine multiply_work_by_conjugate

  subroutine load_convvelo_field_to_zbuf(component_index)
    implicit none
    integer(C_INT), intent(in) :: component_index
    integer(C_INT) :: iy, iz, ix, jx, izd_idx

    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(VVdz, nzd, nxB, ny) private(iy, jx, izd_idx)
    do iy = 1, ny + 3
      do jx = 1, nxB
        do izd_idx = 1, nzd
          VVdz(izd_idx, jx, iy, 1) = (0.0d0, 0.0d0)
        end do
      end do
    end do
    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(VVdz, V, izd, ny, nz, nx0, nxN, component_index) private(ix, iz, iy, jx)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = -1, ny + 1
          jx = ix - nx0 + 1
          VVdz(izd(iz) + 1, jx, iy + 2, 1) = V(iy, iz, ix, component_index)
        end do
      end do
    end do
  end subroutine load_convvelo_field_to_zbuf

  subroutine spectral_field_to_real_x(rx)
    implicit none
    real(C_DOUBLE), intent(out) :: rx(2*(nxd + 1), nzB, ny + 3)
#ifdef HAVE_MPI
    type(MPI_Request) :: request
    type(MPI_Status) :: status
#endif
    integer(C_INT) :: ix, iz, iy

    call IFT(VVdz(:, :, :, 1), ny)
    call pack_zTOx(VVdz(:, :, :, 1), sendbuf(:, 1), ny)
    call alltoall(sendbuf(:, 1), recvbuf(:, 1), request)
#ifdef HAVE_MPI
    call MPI_Wait(request, status, ierr)
#endif
    call unpack_zTOx(recvbuf(:, 1), VVdx(:, :, :, 1), ny)
    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(VVdx, nx, nxd, nzB, ny) private(ix, iz, iy)
    do iy = 1, ny + 3
      do iz = 1, nzB
        do ix = nx + 2, nxd + 1
          VVdx(ix, iz, iy, 1) = (0.0d0, 0.0d0)
        end do
      end do
    end do
    call RFT(VVdx(:, :, :, 1), rx, ny)
  end subroutine spectral_field_to_real_x

  subroutine accumulate_scaled_product(dst, lhs, rhs, scale)
    implicit none
    real(C_DOUBLE), intent(inout) :: dst(2*(nxd + 1), nzB, ny + 3)
    real(C_DOUBLE), intent(in) :: lhs(2*(nxd + 1), nzB, ny + 3)
    real(C_DOUBLE), intent(in) :: rhs(2*(nxd + 1), nzB, ny + 3)
    real(C_DOUBLE), intent(in) :: scale
    integer(C_INT) :: iy, iz, ix

    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(dst, lhs, rhs, scale, nxd, nzB, ny) private(iy, iz, ix)
    do iy = 1, ny + 3
      do iz = 1, nzB
        do ix = 1, 2*nxd
          dst(ix, iz, iy) = dst(ix, iz, iy) + scale*lhs(ix, iz, iy)*rhs(ix, iz, iy)
        end do
      end do
    end do
  end subroutine accumulate_scaled_product

  subroutine real_x_to_spectral_field(rx, field)
    implicit none
    real(C_DOUBLE), intent(in) :: rx(2*(nxd + 1), nzB, ny + 3)
    complex(C_DOUBLE_COMPLEX), intent(out) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
#ifdef HAVE_MPI
    type(MPI_Request) :: request
    type(MPI_Status) :: status
#endif
    integer(C_INT) :: ix, iz, iy

    call HFT(rx, VVdx(:, :, :, 1), ny)
    call pack_xTOz(VVdx(:, :, :, 1), sendbuf(:, 1), ny)
    call alltoall(sendbuf(:, 1), recvbuf(:, 1), request)
#ifdef HAVE_MPI
    call MPI_Wait(request, status, ierr)
#endif
    call unpack_xTOz(recvbuf(:, 1), VVdz(:, :, :, 1), ny)
    call FFT(VVdz(:, :, :, 1), ny)

    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(VVdz, nx0, nxN, ny, nz, field) private(ix, iz, iy)
    do ix = nx0, nxN
      do iy = -1, ny + 1
        do iz = 0, nz
          field(iy, iz, ix) = VVdz(iz + 1, ix - nx0 + 1, iy + 2, 1)
        end do
      end do
    end do
    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(VVdz, nx0, nxN, ny, nz, field, izd) private(ix, iz, iy)
    do ix = nx0, nxN
      do iy = -1, ny + 1
        do iz = -nz, -1
          field(iy, iz, ix) = VVdz(izd(iz) + 1, ix - nx0 + 1, iy + 2, 1)
        end do
      end do
    end do
  end subroutine real_x_to_spectral_field

  subroutine accumulate_cross_components(field_index, lhs_component, rhs_component)
    implicit none
    integer(C_INT), intent(in) :: field_index, lhs_component, rhs_component

    call load_component_to_work(rhs_component)
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_components

  subroutine accumulate_cross_derivative(field_index, lhs_component, rhs_component, derivative_order)
    implicit none
    integer(C_INT), intent(in) :: field_index, lhs_component, rhs_component, derivative_order

    select case (derivative_order)
    case (1)
      call apply_dy_to_work(rhs_component)
    case (2)
      call apply_dyy_to_work(rhs_component)
    case default
      error stop "accumulate_cross_derivative: unsupported derivative order"
    end select

    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_derivative

  subroutine accumulate_cross_pressure(field_index, lhs_component, use_dpdy)
    implicit none
    integer(C_INT), intent(in) :: field_index, lhs_component
    logical, intent(in) :: use_dpdy

    if (use_dpdy) then
      call compute_dpdy(convvelo_work)
    else
      call compute_poisson(convvelo_work)
    end if
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_pressure

  subroutine accumulate_cross_product_field(field_index, lhs_component, product_case)
    implicit none
    integer(C_INT), intent(in) :: field_index, lhs_component, product_case
    integer(C_INT) :: rhs0, rhs1
    integer(C_INT) :: iy, iz, ix

    select case (product_case)
    case (1)
      rhs0 = 1; rhs1 = 1
    case (2)
      rhs0 = 2; rhs1 = 2
    case (3)
      rhs0 = 3; rhs1 = 3
    case (4)
      rhs0 = 1; rhs1 = 2
    case (5)
      rhs0 = 2; rhs1 = 3
    case (6)
      rhs0 = 1; rhs1 = 3
    case default
      if (product_case < 7) error stop "accumulate_cross_product_field: unsupported product_case"
      rhs0 = mod(product_case - 4, 3) + 1
      rhs1 = 3 + (product_case - 4)/3
    end select

    call load_convvelo_field_to_zbuf(rhs0)
    call spectral_field_to_real_x(convvelo_real0)
    call load_convvelo_field_to_zbuf(rhs1)
    call spectral_field_to_real_x(convvelo_real1)

    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(convvelo_real_prod, nxd, nzB, ny) private(iy, iz, ix)
    do iy = 1, ny + 3
      do iz = 1, nzB
        do ix = 1, 2*(nxd + 1)
          convvelo_real_prod(ix, iz, iy) = 0.0d0
        end do
      end do
    end do

    call accumulate_scaled_product(convvelo_real_prod, convvelo_real0, convvelo_real1, factor)
    call real_x_to_spectral_field(convvelo_real_prod, convvelo_work)
    !$omp target update from(convvelo_work)
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_product_field

  subroutine accumulate_cross_product_derivative_field(field_index, lhs_component, product_case)
    implicit none
    integer(C_INT), intent(in) :: field_index, lhs_component, product_case

    call build_cross_product_work(product_case)
    call apply_dy_to_existing_work()
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_product_derivative_field

  subroutine build_cross_product_work(product_case)
    implicit none
    integer(C_INT), intent(in) :: product_case
    integer(C_INT) :: rhs0, rhs1
    integer(C_INT) :: iy, iz, ix

    select case (product_case)
    case (1)
      rhs0 = 1; rhs1 = 1
    case (2)
      rhs0 = 2; rhs1 = 2
    case (3)
      rhs0 = 3; rhs1 = 3
    case (4)
      rhs0 = 1; rhs1 = 2
    case (5)
      rhs0 = 2; rhs1 = 3
    case (6)
      rhs0 = 1; rhs1 = 3
    case default
      if (product_case < 7) error stop "build_cross_product_work: unsupported product_case"
      rhs0 = mod(product_case - 4, 3) + 1
      rhs1 = 3 + (product_case - 4)/3
    end select

    call load_convvelo_field_to_zbuf(rhs0)
    call spectral_field_to_real_x(convvelo_real0)
    call load_convvelo_field_to_zbuf(rhs1)
    call spectral_field_to_real_x(convvelo_real1)

    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(convvelo_real_prod, nxd, nzB, ny) private(iy, iz, ix)
    do iy = 1, ny + 3
      do iz = 1, nzB
        do ix = 1, 2*(nxd + 1)
          convvelo_real_prod(ix, iz, iy) = 0.0d0
        end do
      end do
    end do

    call accumulate_scaled_product(convvelo_real_prod, convvelo_real0, convvelo_real1, factor)
    call real_x_to_spectral_field(convvelo_real_prod, convvelo_work)
    !$omp target update from(convvelo_work)
  end subroutine build_cross_product_work

  subroutine apply_dy_to_existing_work()
    implicit none
    complex(C_DOUBLE_COMPLEX) :: tmp(ny0 - 2:nyN + 2)
    integer(C_INT) :: iz, ix

    do ix = nx0, nxN
      do iz = -nz, nz
        tmp = convvelo_work(:, iz, ix)
        call COMPLEXderiv(tmp, convvelo_work(:, iz, ix), der, D0mat)
      end do
    end do
  end subroutine apply_dy_to_existing_work

  subroutine free_convvelo()
    implicit none

    if (.not. convvelo_initialized) return

    !$omp target exit data map(delete: convvelo_stats, convvelo_work, component_means, convvelo_real0, convvelo_real1, convvelo_real_prod)
    deallocate (convvelo_stats, convvelo_work, component_means, convvelo_real0, convvelo_real1, convvelo_real_prod, n_field_samples)

    n_convvelo_fields = 0
    n_mean_samples = 0_C_INT64_T
    convvelo_initialized = .false.
  end subroutine free_convvelo

end module convvelo
