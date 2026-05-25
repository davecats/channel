#include "header.h"

module convvelo

  use, intrinsic :: iso_c_binding
  use dnsdata, only: V, nPhi, nz, ny, der, nxd, izd, factor, iproc, D0mat, d240, d24m1, d24n, d24np1, &
                     COMPLEXderiv, LeftLU5div
  use pressure_output, only: compute_poisson, compute_dpdy
  use mpi_transpose, only: ny0, nyN, nx0, nxN, nxB, nzB, nx, has_average, ierr, sendbuf, recvbuf, &
                           pack_zTOx, unpack_zTOx, pack_xTOz, unpack_xTOz, alltoall, nzd
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  use ffts, only: IFT, RFT, HFT, FFT, VVdx, VVdz
#else
  use dnsdata, only: VVdx, VVdz
  use ffts, only: IFT, RFT, HFT, FFT
#endif
#ifdef HAVE_MPI
  use mpi_f08
#endif

  implicit none

  private

  integer(C_INT), parameter, public :: n_convvelo_velocity_fields = 33
  integer(C_INT), parameter, public :: n_convvelo_scalar_fields = 10
  integer(C_INT), parameter, public :: n_convvelo_velocity_fields_minimal = 20
  integer(C_INT), parameter, public :: n_convvelo_scalar_fields_minimal = 9
  integer(C_INT), parameter :: i_u = 1
  integer(C_INT), parameter :: i_v = 2
  integer(C_INT), parameter :: i_w = 3
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
  integer(C_INT), parameter :: i_u_cross_p = 22
  integer(C_INT), parameter :: i_v_cross_dpdy = 23
  integer(C_INT), parameter :: i_w_cross_p = 24
  integer(C_INT), parameter :: i_u_cross_uu = 25
  integer(C_INT), parameter :: i_u_cross_uw = 26
  integer(C_INT), parameter :: i_v_cross_uv = 27
  integer(C_INT), parameter :: i_v_cross_vw = 28
  integer(C_INT), parameter :: i_w_cross_uw = 29
  integer(C_INT), parameter :: i_w_cross_ww = 30
  integer(C_INT), parameter :: i_u_cross_dyuv = 31
  integer(C_INT), parameter :: i_v_cross_dyvv = 32
  integer(C_INT), parameter :: i_w_cross_dyvw = 33
  integer(C_INT), parameter :: i_t_theta_theta = 1
  integer(C_INT), parameter :: i_t_theta_u = 2
  integer(C_INT), parameter :: i_t_theta_v = 3
  integer(C_INT), parameter :: i_t_theta_w = 4
  integer(C_INT), parameter :: i_t_theta_thetau = 5
  integer(C_INT), parameter :: i_t_theta_thetaw = 6
  integer(C_INT), parameter :: i_t_theta_dyytheta = 7
  integer(C_INT), parameter :: i_t_theta_dythetav = 8
  integer(C_INT), parameter :: i_t_theta_dytheta = 9
  integer(C_INT), parameter :: i_t_theta_dyv = 10
  integer(C_INT), parameter :: i_prod_uu = 1
  integer(C_INT), parameter :: i_prod_vv = 2
  integer(C_INT), parameter :: i_prod_ww = 3
  integer(C_INT), parameter :: i_prod_uv = 4
  integer(C_INT), parameter :: i_prod_vw = 5
  integer(C_INT), parameter :: i_prod_uw = 6
  integer(C_INT), parameter :: minimal_velocity_fields(n_convvelo_velocity_fields_minimal) = [ &
                               i_u_cross_u, i_u_cross_v, i_v_cross_v, i_w_cross_w, &
                               i_u_cross_p, i_v_cross_dpdy, i_w_cross_p, &
                               i_u_cross_uu, i_u_cross_uw, i_v_cross_uv, i_v_cross_vw, i_w_cross_uw, i_w_cross_ww, &
                               i_u_cross_dyyu, i_v_cross_dyyv, i_w_cross_dyyw, &
                               i_u_cross_dyuv, i_v_cross_dyvv, i_w_cross_dyvw, i_u_cross_dyv &
                               ]
  integer(C_INT), parameter :: minimal_scalar_fields(n_convvelo_scalar_fields_minimal) = [ &
                               i_t_theta_theta, i_t_theta_u, i_t_theta_v, i_t_theta_w, &
                               i_t_theta_thetau, i_t_theta_thetaw, i_t_theta_dyytheta, i_t_theta_dythetav, i_t_theta_dyv &
                               ]

  logical, save :: convvelo_initialized = .false.
  logical, save :: convvelo_dirty = .false.
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
  public :: acc_convvelo_stats, write_convvelo_output, convvelo_has_pending_output, sync_convvelo_output_to_host

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
    convvelo_dirty = .false.
    convvelo_initialized = .true.
  end subroutine init_convvelo

  subroutine reset_convvelo_stats()
    implicit none

    if (.not. convvelo_initialized) return

    call zero_convvelo_stats()
    call zero_convvelo_work()
    call zero_component_means()
    call zero_real_buffer(convvelo_real0)
    call zero_real_buffer(convvelo_real1)
    call zero_real_buffer(convvelo_real_prod)
    n_field_samples = 0_C_INT64_T

    n_mean_samples = 0_C_INT64_T
    convvelo_dirty = .false.
  end subroutine reset_convvelo_stats

  subroutine update_convvelo_component_means()
    implicit none

    complex(C_DOUBLE_COMPLEX) :: snapshot(ny0 - 2:nyN + 2, 1:3 + nPhi)
    real(C_DOUBLE) :: old_weight, new_weight
#ifndef HAVE_MPI
    integer(C_INT) :: iy, ic
#endif

    if (.not. convvelo_initialized) call init_convvelo()

    snapshot = (0.0d0, 0.0d0)
    if (has_average) then
#ifdef HAVE_MPI
      !$omp target update from(V(ny0 - 2:nyN + 2, 0, 0, 1:3 + nPhi))
      snapshot(:, :) = V(:, 0, 0, :)
#else
      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp shared(snapshot, V) private(ic, iy)
      do ic = 1, 3 + nPhi
        do iy = ny0 - 2, nyN + 2
          snapshot(iy, ic) = V(iy, 0, 0, ic)
        end do
      end do
#endif
    end if

#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, snapshot, size(snapshot), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    n_mean_samples = n_mean_samples + 1_C_INT64_T
    old_weight = dble(n_mean_samples - 1_C_INT64_T)/dble(n_mean_samples)
    new_weight = 1.0d0/dble(n_mean_samples)
#ifdef HAVE_MPI
    component_means = old_weight*component_means + new_weight*snapshot
    !$omp target update to(component_means)
#else
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(component_means, snapshot, old_weight, new_weight) private(ic, iy)
    do ic = 1, 3 + nPhi
      do iy = ny0 - 2, nyN + 2
        component_means(iy, ic) = old_weight*component_means(iy, ic) + new_weight*snapshot(iy, ic)
      end do
    end do
#endif
    convvelo_dirty = .true.
  end subroutine update_convvelo_component_means

  subroutine acc_convvelo_stats()
    implicit none
    integer(C_INT) :: iPhi, scalar_component, scalar_offset

    if (.not. convvelo_initialized) call init_convvelo()
    call accumulate_cross_components(i_u_cross_u, i_u, i_u)
    call accumulate_cross_derivative(i_u_cross_dyu, i_u, i_u, 1)
    call accumulate_cross_components(i_u_cross_v, i_u, i_v)
    call accumulate_cross_derivative(i_u_cross_dyv, i_u, i_v, 1)
    call accumulate_cross_components(i_u_cross_w, i_u, i_w)
    call accumulate_cross_derivative(i_u_cross_dyw, i_u, i_w, 1)
    call accumulate_cross_derivative(i_u_cross_dyyu, i_u, i_u, 2)

    call accumulate_cross_components(i_v_cross_u, i_v, i_u)
    call accumulate_cross_derivative(i_v_cross_dyu, i_v, i_u, 1)
    call accumulate_cross_components(i_v_cross_v, i_v, i_v)
    call accumulate_cross_derivative(i_v_cross_dyv, i_v, i_v, 1)
    call accumulate_cross_components(i_v_cross_w, i_v, i_w)
    call accumulate_cross_derivative(i_v_cross_dyw, i_v, i_w, 1)
    call accumulate_cross_derivative(i_v_cross_dyyv, i_v, i_v, 2)

    call accumulate_cross_components(i_w_cross_u, i_w, i_u)
    call accumulate_cross_derivative(i_w_cross_dyu, i_w, i_u, 1)
    call accumulate_cross_components(i_w_cross_v, i_w, i_v)
    call accumulate_cross_derivative(i_w_cross_dyv, i_w, i_v, 1)
    call accumulate_cross_components(i_w_cross_w, i_w, i_w)
    call accumulate_cross_derivative(i_w_cross_dyw, i_w, i_w, 1)
    call accumulate_cross_derivative(i_w_cross_dyyw, i_w, i_w, 2)

    do iPhi = 1, nPhi
      scalar_component = scalar_component_index(iPhi)
      scalar_offset = n_convvelo_velocity_fields + (iPhi - 1)*n_convvelo_scalar_fields
      call accumulate_cross_components(scalar_offset + i_t_theta_theta, scalar_component, scalar_component)
      call accumulate_cross_components(scalar_offset + i_t_theta_u, scalar_component, i_u)
      call accumulate_cross_components(scalar_offset + i_t_theta_v, scalar_component, i_v)
      call accumulate_cross_components(scalar_offset + i_t_theta_w, scalar_component, i_w)
      call accumulate_cross_product_field(scalar_offset + i_t_theta_thetau, scalar_component, scalar_product_case(i_u, iPhi))
      call accumulate_cross_product_field(scalar_offset + i_t_theta_thetaw, scalar_component, scalar_product_case(i_w, iPhi))
      call accumulate_cross_derivative(scalar_offset + i_t_theta_dyytheta, scalar_component, scalar_component, 2)
call accumulate_cross_product_derivative_field(scalar_offset + i_t_theta_dythetav, scalar_component, scalar_product_case(i_v, iPhi))
      call accumulate_cross_derivative(scalar_offset + i_t_theta_dytheta, scalar_component, scalar_component, 1)
      call accumulate_cross_derivative(scalar_offset + i_t_theta_dyv, scalar_component, i_v, 1)
    end do
    call accumulate_cross_product_field(i_u_cross_uu, i_u, i_prod_uu)
    call accumulate_cross_product_field(i_u_cross_uw, i_u, i_prod_uw)
    call accumulate_cross_product_field(i_v_cross_uv, i_v, i_prod_uv)
    call accumulate_cross_product_field(i_v_cross_vw, i_v, i_prod_vw)
    call accumulate_cross_product_field(i_w_cross_uw, i_w, i_prod_uw)
    call accumulate_cross_product_field(i_w_cross_ww, i_w, i_prod_ww)
    call accumulate_cross_product_derivative_field(i_u_cross_dyuv, i_u, i_prod_uv)
    call accumulate_cross_product_derivative_field(i_v_cross_dyvv, i_v, i_prod_vv)
    call accumulate_cross_product_derivative_field(i_w_cross_dyvw, i_w, i_prod_vw)

    call accumulate_cross_pressure(i_u_cross_p, i_u, .false.)
    call accumulate_cross_pressure(i_v_cross_dpdy, i_v, .true.)
    call accumulate_cross_pressure(i_w_cross_p, i_w, .false.)
  end subroutine acc_convvelo_stats

  subroutine start_convvelo_field()
    implicit none

    if (.not. convvelo_initialized) call init_convvelo()

    call zero_convvelo_work()
  end subroutine start_convvelo_field

  subroutine finish_convvelo_field(field_index)
    implicit none

    integer(C_INT), intent(in) :: field_index
    real(C_DOUBLE) :: old_weight, new_weight
    integer(C_INT) :: ix, iy, iz

    if (.not. convvelo_initialized) call init_convvelo()
    if (field_index < 1 .or. field_index > n_convvelo_fields) then
      error stop "finish_convvelo_field: field_index out of range"
    end if

    n_field_samples(field_index) = n_field_samples(field_index) + 1_C_INT64_T
    old_weight = dble(n_field_samples(field_index) - 1_C_INT64_T)/dble(n_field_samples(field_index))
    new_weight = 1.0d0/dble(n_field_samples(field_index))
    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(convvelo_stats, convvelo_work, field_index, old_weight, new_weight, ny0, nyN, nz, nx0, nxN) private(ix, iz, iy)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0 - 2, nyN + 2
          convvelo_stats(iy, iz, ix, field_index) = old_weight*convvelo_stats(iy, iz, ix, field_index) + &
                                                    new_weight*convvelo_work(iy, iz, ix)
        end do
      end do
    end do
  end subroutine finish_convvelo_field

  subroutine load_component_to_work(component_index)
    implicit none
    integer(C_INT), intent(in) :: component_index
    integer(C_INT) :: iy, iz, ix

    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(convvelo_work, V, component_index, ny0, nyN, nz, nx0, nxN) private(ix, iz, iy)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0 - 2, nyN + 2
          convvelo_work(iy, iz, ix) = V(iy, iz, ix, component_index)
        end do
      end do
    end do
  end subroutine load_component_to_work

  subroutine apply_dy_to_work(component_index)
    implicit none
    integer(C_INT), intent(in) :: component_index
    integer(C_INT) :: iz, ix

    call zero_convvelo_work()
    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(V, convvelo_work, component_index, der, D0mat, nx0, nxN, nz) private(ix, iz)
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

    call zero_convvelo_work()
    !$omp target teams distribute parallel do collapse(2) &
    !$omp shared(convvelo_work, V, component_index, d240, d24m1, d24n, d24np1, der, D0mat, ny0, nyN, ny, nx0, nxN, nz) private(ix, iz, iy)
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

    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(convvelo_work, V, component_index, ny0, nyN, nx0, nxN, nz) private(ix, iz, iy)
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

    call decode_product_case(product_case, rhs0, rhs1, "accumulate_cross_product_field")
    call build_cross_product_work(rhs0, rhs1)
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_product_field

  subroutine accumulate_cross_product_derivative_field(field_index, lhs_component, product_case)
    implicit none
    integer(C_INT), intent(in) :: field_index, lhs_component, product_case

    call build_cross_product_work_from_case(product_case)
    call apply_dy_to_existing_work()
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_product_derivative_field

  subroutine build_cross_product_work_from_case(product_case)
    implicit none
    integer(C_INT), intent(in) :: product_case
    integer(C_INT) :: rhs0, rhs1

    call decode_product_case(product_case, rhs0, rhs1, "build_cross_product_work")
    call build_cross_product_work(rhs0, rhs1)
  end subroutine build_cross_product_work_from_case

  subroutine build_cross_product_work(rhs0, rhs1)
    implicit none
    integer(C_INT), intent(in) :: rhs0, rhs1

    call load_convvelo_field_to_zbuf(rhs0)
    call spectral_field_to_real_x(convvelo_real0)
    call load_convvelo_field_to_zbuf(rhs1)
    call spectral_field_to_real_x(convvelo_real1)
    call zero_real_buffer(convvelo_real_prod)
    call accumulate_scaled_product(convvelo_real_prod, convvelo_real0, convvelo_real1, factor)
    call real_x_to_spectral_field(convvelo_real_prod, convvelo_work)
  end subroutine build_cross_product_work

  subroutine apply_dy_to_existing_work()
    implicit none
    complex(C_DOUBLE_COMPLEX) :: tmp(ny0 - 2:nyN + 2)
    integer(C_INT) :: iz, ix

    !$omp target teams distribute parallel do collapse(2) default(none) &
    !$omp shared(convvelo_work, der, D0mat, nx0, nxN, nz) private(ix, iz, tmp)
    do ix = nx0, nxN
      do iz = -nz, nz
        tmp = convvelo_work(:, iz, ix)
        call COMPLEXderiv(tmp, convvelo_work(:, iz, ix), der, D0mat)
      end do
    end do
  end subroutine apply_dy_to_existing_work

  subroutine zero_convvelo_work()
    implicit none
    integer(C_INT) :: iy, iz, ix

    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(convvelo_work, ny0, nyN, nx0, nxN, nz) private(ix, iz, iy)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0 - 2, nyN + 2
          convvelo_work(iy, iz, ix) = (0.0d0, 0.0d0)
        end do
      end do
    end do
  end subroutine zero_convvelo_work

  subroutine zero_convvelo_stats()
    implicit none
    integer(C_INT) :: field_index, iy, iz, ix

    !$omp target teams distribute parallel do collapse(4) &
    !$omp shared(convvelo_stats, n_convvelo_fields, ny0, nyN, nx0, nxN, nz) private(field_index, ix, iz, iy)
    do field_index = 1, n_convvelo_fields
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = ny0 - 2, nyN + 2
            convvelo_stats(iy, iz, ix, field_index) = (0.0d0, 0.0d0)
          end do
        end do
      end do
    end do
  end subroutine zero_convvelo_stats

  subroutine zero_component_means()
    implicit none
    integer(C_INT) :: iy, ic

    !$omp target teams distribute parallel do collapse(2) &
    !$omp shared(component_means, ny0, nyN, nPhi) private(ic, iy)
    do ic = 1, 3 + nPhi
      do iy = ny0 - 2, nyN + 2
        component_means(iy, ic) = (0.0d0, 0.0d0)
      end do
    end do
  end subroutine zero_component_means

  subroutine zero_real_buffer(buffer)
    implicit none
    real(C_DOUBLE), intent(inout) :: buffer(2*(nxd + 1), nzB, ny + 3)
    integer(C_INT) :: iy, iz, ix

    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(buffer, nxd, nzB, ny) private(ix, iz, iy)
    do iy = 1, ny + 3
      do iz = 1, nzB
        do ix = 1, 2*(nxd + 1)
          buffer(ix, iz, iy) = 0.0d0
        end do
      end do
    end do
  end subroutine zero_real_buffer

  subroutine free_convvelo()
    implicit none

    if (.not. convvelo_initialized) return

    !$omp target exit data map(delete: convvelo_stats, convvelo_work, component_means, convvelo_real0, convvelo_real1, convvelo_real_prod)
    deallocate (convvelo_stats, convvelo_work, component_means, convvelo_real0, convvelo_real1, convvelo_real_prod, n_field_samples)

    n_convvelo_fields = 0
    n_mean_samples = 0_C_INT64_T
    convvelo_dirty = .false.
    convvelo_initialized = .false.
  end subroutine free_convvelo

  logical function convvelo_has_pending_output()
    implicit none

    convvelo_has_pending_output = convvelo_initialized .and. convvelo_dirty .and. &
                                  (n_mean_samples > 0_C_INT64_T .or. any(n_field_samples > 0_C_INT64_T))
  end function convvelo_has_pending_output

  subroutine write_convvelo_output(filename, write_full_fields)
    implicit none

    character(len=*), intent(in) :: filename
    logical, intent(in) :: write_full_fields

    if (.not. convvelo_has_pending_output()) return

    if (write_full_fields) then
      call write_convvelo_raw_stats(filename)
    else
      call write_convvelo_component_means(filename)
    end if

    convvelo_dirty = .false.
  end subroutine write_convvelo_output

  subroutine sync_convvelo_output_to_host()
    implicit none

    if (.not. convvelo_initialized) return

    !$omp target update from(component_means, convvelo_stats)
  end subroutine sync_convvelo_output_to_host

  subroutine write_convvelo_raw_stats(filename)
    implicit none

    character(len=*), intent(in) :: filename
    integer(C_INT) :: field_index, iPhi

#ifdef HAVE_MPI
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    type(MPI_Datatype) :: file_type, mem_type, profile_file_type, profile_mem_type
    integer :: ierror
    integer, parameter :: ndims = 3
    integer, parameter :: ndims_profile = 1
    integer :: sizes(ndims), subsizes(ndims), starts(ndims)
    integer :: profile_sizes(ndims_profile), profile_subsizes(ndims_profile), profile_starts(ndims_profile)
    integer(MPI_OFFSET_KIND) :: disp, field_bytes, profile_bytes, total_bytes
#else
    integer :: io
#endif

    if (.not. convvelo_initialized) return

    !$omp target update from(component_means, convvelo_stats)

#ifdef HAVE_MPI
    sizes = [ny + 3, 2*nz + 1, nx + 1]
    subsizes = [ny + 3, 2*nz + 1, nxN - nx0 + 1]
    starts = [0, 0, nx0]
    call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, file_type, ierror)
    call MPI_Type_commit(file_type, ierror)

    sizes = [ny + 3, 2*nz + 1, nxN - nx0 + 1]
    subsizes = sizes
    starts = [0, 0, 0]
    call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, mem_type, ierror)
    call MPI_Type_commit(mem_type, ierror)

    profile_sizes = [ny + 3]
    profile_subsizes = [nyN - ny0 + 5]
    profile_starts = [ny0 - 1]
    call MPI_Type_create_subarray(ndims_profile, profile_sizes, profile_subsizes, profile_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, profile_file_type, ierror)
    call MPI_Type_commit(profile_file_type, ierror)

    profile_sizes = [nyN - ny0 + 5]
    profile_subsizes = profile_sizes
    profile_starts = [0]
    call MPI_Type_create_subarray(ndims_profile, profile_sizes, profile_subsizes, profile_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, profile_mem_type, ierror)
    call MPI_Type_commit(profile_mem_type, ierror)

    profile_bytes = int(16, MPI_OFFSET_KIND)*int(ny + 3, MPI_OFFSET_KIND)
    field_bytes = int(16, MPI_OFFSET_KIND)*int(ny + 3, MPI_OFFSET_KIND)* &
                  int(2*nz + 1, MPI_OFFSET_KIND)*int(nx + 1, MPI_OFFSET_KIND)
    total_bytes = int(1 + nPhi, MPI_OFFSET_KIND)*profile_bytes + int(n_convvelo_fields, MPI_OFFSET_KIND)*field_bytes

    call MPI_File_open(MPI_COMM_WORLD, trim(filename), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
    call MPI_File_set_size(fh, total_bytes)

    disp = 0_MPI_OFFSET_KIND
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, profile_file_type, 'native', MPI_INFO_NULL)
    call MPI_File_write_all(fh, component_means(:, 1), 1, profile_mem_type, status)
    do iPhi = 1, nPhi
      disp = int(iPhi, MPI_OFFSET_KIND)*profile_bytes
      call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, profile_file_type, 'native', MPI_INFO_NULL)
      call MPI_File_write_all(fh, component_means(:, 3 + iPhi), 1, profile_mem_type, status)
    end do

    do field_index = 1, n_convvelo_fields
      disp = int(1 + nPhi, MPI_OFFSET_KIND)*profile_bytes + int(field_index - 1, MPI_OFFSET_KIND)*field_bytes
      call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, file_type, 'native', MPI_INFO_NULL)
      call MPI_File_write_all(fh, convvelo_stats(:, :, :, field_index), 1, mem_type, status)
    end do

    call MPI_File_close(fh)
    call MPI_Type_free(file_type, ierror)
    call MPI_Type_free(mem_type, ierror)
    call MPI_Type_free(profile_file_type, ierror)
    call MPI_Type_free(profile_mem_type, ierror)
#else
    open (unit=99, file=trim(filename), form='unformatted', access='stream', status='replace', action='write', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open convvelo output file: ', trim(filename)
      stop 1
    end if

    write (99) component_means(:, 1)
    do iPhi = 1, nPhi
      write (99) component_means(:, 3 + iPhi)
    end do
    do field_index = 1, n_convvelo_fields
      write (99) convvelo_stats(:, :, :, field_index)
    end do
    close (99)
#endif
  end subroutine write_convvelo_raw_stats

  subroutine write_convvelo_component_means(filename)
    implicit none

    character(len=*), intent(in) :: filename
    integer :: io, iPhi
    integer(C_INT) :: scalar_offset

    if (.not. convvelo_initialized) return
    if (iproc /= 0) return

    !$omp target update from(component_means, convvelo_stats)

    open (unit=98, file=trim(filename), form='unformatted', access='stream', status='replace', action='write', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open convvelo means output file: ', trim(filename)
      stop 1
    end if

    write (98) component_means(:, 1)
    do iPhi = 1, nPhi
      write (98) component_means(:, 3 + iPhi)
    end do

    call write_selected_fields(98, 0_C_INT, minimal_velocity_fields)
    do iPhi = 1, nPhi
      scalar_offset = n_convvelo_velocity_fields + (iPhi - 1)*n_convvelo_scalar_fields
      call write_selected_fields(98, scalar_offset, minimal_scalar_fields)
    end do
    close (98)
  end subroutine write_convvelo_component_means

  subroutine decode_product_case(product_case, rhs0, rhs1, routine_name)
    implicit none
    integer(C_INT), intent(in) :: product_case
    integer(C_INT), intent(out) :: rhs0, rhs1
    character(len=*), intent(in) :: routine_name

    select case (product_case)
    case (i_prod_uu)
      rhs0 = i_u
      rhs1 = i_u
    case (i_prod_vv)
      rhs0 = i_v
      rhs1 = i_v
    case (i_prod_ww)
      rhs0 = i_w
      rhs1 = i_w
    case (i_prod_uv)
      rhs0 = i_u
      rhs1 = i_v
    case (i_prod_vw)
      rhs0 = i_v
      rhs1 = i_w
    case (i_prod_uw)
      rhs0 = i_u
      rhs1 = i_w
    case default
      if (product_case < 7) error stop trim(routine_name)//": unsupported product_case"
      rhs0 = mod(product_case - 4, 3) + 1
      rhs1 = 3 + (product_case - 4)/3
    end select
  end subroutine decode_product_case

  integer(C_INT) function scalar_component_index(iPhi)
    implicit none
    integer(C_INT), intent(in) :: iPhi

    scalar_component_index = i_w + iPhi
  end function scalar_component_index

  integer(C_INT) function scalar_product_case(velocity_component, iPhi)
    implicit none
    integer(C_INT), intent(in) :: velocity_component, iPhi

    scalar_product_case = i_prod_uw + 3*(iPhi - 1) + velocity_component
  end function scalar_product_case

  subroutine write_selected_fields(io_unit, field_offset, field_indices)
    implicit none
    integer, intent(in) :: io_unit
    integer(C_INT), intent(in) :: field_offset
    integer(C_INT), intent(in) :: field_indices(:)
    integer(C_INT) :: i

    do i = 1, size(field_indices)
      write (io_unit) convvelo_stats(:, :, :, field_offset + field_indices(i))
    end do
  end subroutine write_selected_fields

end module convvelo
