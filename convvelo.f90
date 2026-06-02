#include "header.h"

module convvelo

  use, intrinsic :: iso_c_binding
  use config, only: ini_config, has_section, get_string, get_real, lower
  use dnsdata, only: V, nPhi, nz, ny, der, nxd, izd, factor, iproc, D0mat, d240, d24m1, d24n, d24np1, &
                     apply_complex_derivative_with_y_pencil, has_terminal, &
                     time, deltat
  use pressure_output, only: compute_poisson, compute_dpdy
  use mpi_transpose, only: ny0, nyN, nx0, nxN, nxB, nzB, nx, has_average, ierr, sendbuf, recvbuf, &
                           pack_zTOx, unpack_zTOx, pack_xTOz, unpack_xTOz, alltoall, nzd, fft_transpose_is_local, &
                           repack_zTOx_local, repack_xTOz_local
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

  integer(C_INT), parameter :: n_convvelo_velocity_fields_total = 33
  integer(C_INT), parameter :: n_convvelo_scalar_fields_total = 10
  integer(C_INT), parameter :: n_convvelo_velocity_fields_minimal_total = 21
  integer(C_INT), parameter :: n_convvelo_scalar_fields_minimal_total = 9
  integer(C_INT), parameter :: n_convvelo_profile_fields = 4
  integer(C_INT), parameter :: i_u = 1
  integer(C_INT), parameter :: i_v = 2
  integer(C_INT), parameter :: i_w = 3
  integer(C_INT64_T), parameter :: convvelo_file_header_bytes = 2_C_INT64_T*8_C_INT64_T + 8_C_INT64_T
  character(len=*), parameter :: convvelo_runtime_filename = "convvelo.bin"
  character(len=32), parameter :: velocity_field_names(n_convvelo_velocity_fields_total) = [character(len=32) :: &
                               "u_cross_u", "u_cross_dyu", "u_cross_v", "u_cross_dyv", "u_cross_w", "u_cross_dyw", "u_cross_dyyu", &
                               "v_cross_u", "v_cross_dyu", "v_cross_v", "v_cross_dyv", "v_cross_w", "v_cross_dyw", "v_cross_dyyv", &
                               "w_cross_u", "w_cross_dyu", "w_cross_v", "w_cross_dyv", "w_cross_w", "w_cross_dyw", "w_cross_dyyw", &
                                 "u_cross_p", "v_cross_dpdy", "w_cross_p", "u_cross_uu", "u_cross_uw", "v_cross_uv", "v_cross_vw", &
                                                         "w_cross_uw", "w_cross_ww", "u_cross_dyuv", "v_cross_dyvv", "w_cross_dyvw"]
  character(len=32), parameter :: scalar_field_names(n_convvelo_scalar_fields_total) = [character(len=32) :: &
                                                                 "t_cross_t", "t_cross_u", "t_cross_v", "t_cross_w", "t_cross_tu", &
                                                         "t_cross_tw", "t_cross_dyyt", "t_cross_dytv", "t_cross_dyt", "t_cross_dyv"]
  character(len=32), parameter :: minimal_velocity_field_names(n_convvelo_velocity_fields_minimal_total) = [character(len=32) :: &
                                                                  "u_cross_u", "u_cross_v", "u_cross_w", "v_cross_v", "w_cross_w", &
                                                                                         "u_cross_p", "v_cross_dpdy", "w_cross_p", &
                                               "u_cross_uu", "u_cross_uw", "v_cross_uv", "v_cross_vw", "w_cross_uw", "w_cross_ww", &
                                                                                   "u_cross_dyyu", "v_cross_dyyv", "w_cross_dyyw", &
                                                                      "u_cross_dyuv", "v_cross_dyvv", "w_cross_dyvw", "u_cross_dyv"]
  character(len=32), parameter :: minimal_scalar_field_names(n_convvelo_scalar_fields_minimal_total) = [character(len=32) :: &
                                                                               "t_cross_t", "t_cross_u", "t_cross_v", "t_cross_w", &
                                                          "t_cross_tu", "t_cross_tw", "t_cross_dyyt", "t_cross_dytv", "t_cross_dyv"]
  character(len=16), parameter :: convvelo_profile_field_names(n_convvelo_profile_fields) = [character(len=16) :: &
                                                                                             "mean_u", "mean_v", "mean_w", "mean_t"]

  logical, save :: convvelo_initialized = .false.
  logical, save :: convvelo_dirty = .false.
  logical, save, public :: convvelo_enabled = .false.
  logical, save :: convvelo_write_full_fields = .true.
  integer(C_INT), save, public :: n_convvelo_velocity_fields = 0
  integer(C_INT), save, public :: n_convvelo_scalar_fields = 0
  integer(C_INT), save, public :: n_convvelo_fields = 0
  integer(C_INT64_T), save :: n_convvelo_profile_slots = 0_C_INT64_T
  integer(C_INT64_T), save :: convvelo_last_write_index = -1_C_INT64_T
  real(C_DOUBLE), save :: convvelo_t_start = 0.0d0
  real(C_DOUBLE), save :: convvelo_dt_compute = -1.0d0
  real(C_DOUBLE), save :: convvelo_dt_write = -1.0d0
  real(C_DOUBLE), save :: convvelo_average_start_time = 0.0d0
  real(C_DOUBLE), save :: convvelo_average_end_time = 0.0d0
  character(len=16), save :: convvelo_output_mode = "full"
  integer(C_INT), allocatable, save :: convvelo_velocity_field_ids(:)
  integer(C_INT), allocatable, save :: convvelo_scalar_field_ids(:)
  integer(C_INT), allocatable, save :: convvelo_field_map(:)
  integer(C_INT64_T), save :: n_mean_samples = 0_C_INT64_T

  complex(C_DOUBLE_COMPLEX), allocatable, save, public :: convvelo_stats(:, :, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save, public :: convvelo_work(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save, public :: component_means(:, :)
  real(C_DOUBLE), allocatable, save :: convvelo_real0(:, :, :)
  real(C_DOUBLE), allocatable, save :: convvelo_real1(:, :, :)
  real(C_DOUBLE), allocatable, save :: convvelo_real_prod(:, :, :)
  integer(C_INT64_T), allocatable, save :: n_field_samples(:)

  public :: init_convvelo, reset_convvelo_stats, update_convvelo_component_means, free_convvelo
  public :: finish_convvelo_field
  public :: acc_convvelo_stats, convvelo_has_pending_output
  public :: init_convvelo_runtime, advance_convvelo_runtime, finalize_convvelo_runtime
  public :: get_convvelo_memory_estimate, write_convvelo_raw_stats, write_convvelo_runtime_snapshot
  public :: configure_convvelo

contains

  subroutine configure_convvelo(cfg)
    implicit none
    type(ini_config), intent(in) :: cfg
    logical :: found
    character(len=16) :: convvelo_mode

    convvelo_enabled = .false.
    convvelo_write_full_fields = .true.
    convvelo_t_start = time
    convvelo_dt_compute = -1.0d0
    convvelo_dt_write = -1.0d0
    convvelo_output_mode = "full"
    convvelo_average_start_time = 0.0d0
    convvelo_average_end_time = 0.0d0

    convvelo_enabled = has_section(cfg, "convvelo")

    if (.not. convvelo_enabled) return

    convvelo_mode = convvelo_output_mode
    call get_string(cfg, "convvelo", "output_mode", convvelo_mode, found)
    if (found) convvelo_output_mode = trim(convvelo_mode)
    call get_real(cfg, "convvelo", "t_start", convvelo_t_start, found)
    call get_real(cfg, "convvelo", "dt_compute", convvelo_dt_compute, found)
    call get_real(cfg, "convvelo", "dt_write", convvelo_dt_write, found)

    select case (lower(adjustl(trim(convvelo_output_mode))))
    case ("full")
      convvelo_write_full_fields = .true.
      convvelo_output_mode = "full"
    case ("minimal")
      convvelo_write_full_fields = .false.
      convvelo_output_mode = "minimal"
    case default
      error stop "configure_convvelo: convvelo_output_mode must be 'full' or 'minimal'"
    end select

    if (convvelo_dt_compute <= 0.0d0) then
      error stop "configure_convvelo: convvelo_dt_compute must be > 0 when convection velocity output is enabled"
    end if
  end subroutine configure_convvelo

  subroutine get_convvelo_memory_estimate(n_floats)
    implicit none
    integer(C_INT64_T), intent(out) :: n_floats
    integer(C_INT64_T) :: local_y, spectral_planes, real_planes, n_fields

    if (.not. convvelo_enabled) then
      n_floats = 0_C_INT64_T
      return
    end if

    local_y = int(nyN - ny0 + 5, C_INT64_T)
    spectral_planes = local_y*int(2*nz + 1, C_INT64_T)*int(nxN - nx0 + 1, C_INT64_T)
    real_planes = int(2*(nxd + 1), C_INT64_T)*int(nzB, C_INT64_T)*int(ny + 3, C_INT64_T)
    if (convvelo_write_full_fields) then
      n_fields = int(n_convvelo_velocity_fields_total + nPhi*n_convvelo_scalar_fields_total, C_INT64_T)
    else
      n_fields = int(n_convvelo_velocity_fields_minimal_total + nPhi*n_convvelo_scalar_fields_minimal_total, C_INT64_T)
    end if

    n_floats = 0_C_INT64_T
    n_floats = n_floats + 2_C_INT64_T*spectral_planes*n_fields
    n_floats = n_floats + 2_C_INT64_T*spectral_planes
    n_floats = n_floats + 3_C_INT64_T*real_planes
  end subroutine get_convvelo_memory_estimate

  subroutine init_convvelo()
    implicit none

    if (convvelo_initialized) return

    call init_convvelo_field_layout()

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
    convvelo_last_write_index = -1_C_INT64_T
    convvelo_average_start_time = 0.0d0
    convvelo_average_end_time = 0.0d0
    convvelo_dirty = .false.
    convvelo_initialized = .true.
  end subroutine init_convvelo

  subroutine init_convvelo_field_layout()
    implicit none
    integer(C_INT) :: i, raw_field_index

    if (allocated(convvelo_velocity_field_ids)) deallocate (convvelo_velocity_field_ids)
    if (allocated(convvelo_scalar_field_ids)) deallocate (convvelo_scalar_field_ids)
    if (allocated(convvelo_field_map)) deallocate (convvelo_field_map)

    if (convvelo_write_full_fields) then
      allocate (convvelo_velocity_field_ids(n_convvelo_velocity_fields_total))
      allocate (convvelo_scalar_field_ids(n_convvelo_scalar_fields_total))
      convvelo_velocity_field_ids = [(i, i=1, n_convvelo_velocity_fields_total)]
      convvelo_scalar_field_ids = [(i, i=1, n_convvelo_scalar_fields_total)]
    else
      allocate (convvelo_velocity_field_ids(n_convvelo_velocity_fields_minimal_total))
      allocate (convvelo_scalar_field_ids(n_convvelo_scalar_fields_minimal_total))
      do i = 1, size(convvelo_velocity_field_ids)
        convvelo_velocity_field_ids(i) = field_name_index(minimal_velocity_field_names(i), velocity_field_names, "velocity")
      end do
      do i = 1, size(convvelo_scalar_field_ids)
        convvelo_scalar_field_ids(i) = field_name_index(minimal_scalar_field_names(i), scalar_field_names, "scalar")
      end do
    end if

    n_convvelo_velocity_fields = size(convvelo_velocity_field_ids)
    n_convvelo_scalar_fields = size(convvelo_scalar_field_ids)
    n_convvelo_fields = n_convvelo_velocity_fields + nPhi*n_convvelo_scalar_fields
    allocate (convvelo_field_map(n_convvelo_velocity_fields_total + nPhi*n_convvelo_scalar_fields_total))
    convvelo_field_map = 0_C_INT

    do i = 1, size(convvelo_velocity_field_ids)
      convvelo_field_map(convvelo_velocity_field_ids(i)) = i
    end do
    do i = 1, nPhi
      do raw_field_index = 1, size(convvelo_scalar_field_ids)
        convvelo_field_map(n_convvelo_velocity_fields_total + (i - 1)*n_convvelo_scalar_fields_total + convvelo_scalar_field_ids(raw_field_index)) = &
          size(convvelo_velocity_field_ids) + (i - 1)*size(convvelo_scalar_field_ids) + raw_field_index
      end do
    end do

    n_convvelo_profile_slots = int(3 + nPhi, C_INT64_T)
  end subroutine init_convvelo_field_layout

  subroutine reset_convvelo_stats()
    implicit none
    integer(C_INT) :: field_index, iy, iz, ix, ic

    if (.not. convvelo_initialized) return

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

    !$omp target teams distribute parallel do collapse(2) &
    !$omp shared(component_means, ny0, nyN, nPhi) private(ic, iy)
    do ic = 1, 3 + nPhi
      do iy = ny0 - 2, nyN + 2
        component_means(iy, ic) = (0.0d0, 0.0d0)
      end do
    end do
    n_field_samples = 0_C_INT64_T

    n_mean_samples = 0_C_INT64_T
    convvelo_average_start_time = 0.0d0
    convvelo_average_end_time = 0.0d0
    convvelo_dirty = .false.
  end subroutine reset_convvelo_stats

  subroutine update_convvelo_component_means()
    implicit none

    complex(C_DOUBLE_COMPLEX) :: snapshot(ny0 - 2:nyN + 2, 1:3 + nPhi)
    real(C_DOUBLE) :: old_weight, new_weight
    integer(C_INT) :: iy, ic

    if (.not. convvelo_initialized) call init_convvelo()

    snapshot = (0.0d0, 0.0d0)

    if (has_average) then
      !$omp target teams distribute parallel do collapse(2)  &
      !$omp shared(snapshot, V, ny0, nyN, nPhi) private(ic, iy) map(tofrom: snapshot)
      do ic = 1, 3 + nPhi
        do iy = ny0 - 2, nyN + 2
          snapshot(iy, ic) = V(iy, 0, 0, ic)
        end do
      end do
    end if

#ifdef HAVE_MPI
    call MPI_Allreduce(MPI_IN_PLACE, snapshot, size(snapshot), MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    if (n_mean_samples == 0_C_INT64_T) convvelo_average_start_time = time
    n_mean_samples = n_mean_samples + 1_C_INT64_T
    convvelo_average_end_time = time
    old_weight = dble(n_mean_samples - 1_C_INT64_T)/dble(n_mean_samples)
    new_weight = 1.0d0/dble(n_mean_samples)
    !$omp target teams distribute parallel do collapse(2)  &
    !$omp shared(component_means, snapshot, old_weight, new_weight, ny0, nyN, nPhi) private(ic, iy) map(to: snapshot)
    do ic = 1, 3 + nPhi
      do iy = ny0 - 2, nyN + 2
        component_means(iy, ic) = old_weight*component_means(iy, ic) + new_weight*snapshot(iy, ic)
      end do
    end do
    convvelo_dirty = .true.
  end subroutine update_convvelo_component_means

  subroutine acc_convvelo_stats()
    implicit none
    integer(C_INT) :: iPhi, scalar_component

    if (.not. convvelo_initialized) call init_convvelo()
    call accumulate_cross_components(field_name_index("u_cross_u", velocity_field_names, "velocity"), i_u, i_u)
    call accumulate_cross_derivative(field_name_index("u_cross_dyu", velocity_field_names, "velocity"), i_u, i_u, 1)
    call accumulate_cross_components(field_name_index("u_cross_v", velocity_field_names, "velocity"), i_u, i_v)
    call accumulate_cross_derivative(field_name_index("u_cross_dyv", velocity_field_names, "velocity"), i_u, i_v, 1)
    call accumulate_cross_components(field_name_index("u_cross_w", velocity_field_names, "velocity"), i_u, i_w)
    call accumulate_cross_derivative(field_name_index("u_cross_dyw", velocity_field_names, "velocity"), i_u, i_w, 1)
    call accumulate_cross_derivative(field_name_index("u_cross_dyyu", velocity_field_names, "velocity"), i_u, i_u, 2)

    call accumulate_cross_components(field_name_index("v_cross_u", velocity_field_names, "velocity"), i_v, i_u)
    call accumulate_cross_derivative(field_name_index("v_cross_dyu", velocity_field_names, "velocity"), i_v, i_u, 1)
    call accumulate_cross_components(field_name_index("v_cross_v", velocity_field_names, "velocity"), i_v, i_v)
    call accumulate_cross_derivative(field_name_index("v_cross_dyv", velocity_field_names, "velocity"), i_v, i_v, 1)
    call accumulate_cross_components(field_name_index("v_cross_w", velocity_field_names, "velocity"), i_v, i_w)
    call accumulate_cross_derivative(field_name_index("v_cross_dyw", velocity_field_names, "velocity"), i_v, i_w, 1)
    call accumulate_cross_derivative(field_name_index("v_cross_dyyv", velocity_field_names, "velocity"), i_v, i_v, 2)

    call accumulate_cross_components(field_name_index("w_cross_u", velocity_field_names, "velocity"), i_w, i_u)
    call accumulate_cross_derivative(field_name_index("w_cross_dyu", velocity_field_names, "velocity"), i_w, i_u, 1)
    call accumulate_cross_components(field_name_index("w_cross_v", velocity_field_names, "velocity"), i_w, i_v)
    call accumulate_cross_derivative(field_name_index("w_cross_dyv", velocity_field_names, "velocity"), i_w, i_v, 1)
    call accumulate_cross_components(field_name_index("w_cross_w", velocity_field_names, "velocity"), i_w, i_w)
    call accumulate_cross_derivative(field_name_index("w_cross_dyw", velocity_field_names, "velocity"), i_w, i_w, 1)
    call accumulate_cross_derivative(field_name_index("w_cross_dyyw", velocity_field_names, "velocity"), i_w, i_w, 2)

    do iPhi = 1, nPhi
      scalar_component = i_w + iPhi
      call accumulate_cross_components(raw_scalar_field_index(iPhi, "t_cross_t"), scalar_component, scalar_component)
      call accumulate_cross_components(raw_scalar_field_index(iPhi, "t_cross_u"), scalar_component, i_u)
      call accumulate_cross_components(raw_scalar_field_index(iPhi, "t_cross_v"), scalar_component, i_v)
      call accumulate_cross_components(raw_scalar_field_index(iPhi, "t_cross_w"), scalar_component, i_w)
      call accumulate_cross_product_field(raw_scalar_field_index(iPhi, "t_cross_tu"), scalar_component, scalar_component, i_u)
      call accumulate_cross_product_field(raw_scalar_field_index(iPhi, "t_cross_tw"), scalar_component, scalar_component, i_w)
      call accumulate_cross_derivative(raw_scalar_field_index(iPhi, "t_cross_dyyt"), scalar_component, scalar_component, 2)
      call accumulate_cross_product_derivative_field(raw_scalar_field_index(iPhi, "t_cross_dytv"), scalar_component, scalar_component, i_v)
      call accumulate_cross_derivative(raw_scalar_field_index(iPhi, "t_cross_dyt"), scalar_component, scalar_component, 1)
      call accumulate_cross_derivative(raw_scalar_field_index(iPhi, "t_cross_dyv"), scalar_component, i_v, 1)
    end do
    call accumulate_cross_product_field(field_name_index("u_cross_uu", velocity_field_names, "velocity"), i_u, i_u, i_u)
    call accumulate_cross_product_field(field_name_index("u_cross_uw", velocity_field_names, "velocity"), i_u, i_u, i_w)
    call accumulate_cross_product_field(field_name_index("v_cross_uv", velocity_field_names, "velocity"), i_v, i_u, i_v)
    call accumulate_cross_product_field(field_name_index("v_cross_vw", velocity_field_names, "velocity"), i_v, i_v, i_w)
    call accumulate_cross_product_field(field_name_index("w_cross_uw", velocity_field_names, "velocity"), i_w, i_u, i_w)
    call accumulate_cross_product_field(field_name_index("w_cross_ww", velocity_field_names, "velocity"), i_w, i_w, i_w)
   call accumulate_cross_product_derivative_field(field_name_index("u_cross_dyuv", velocity_field_names, "velocity"), i_u, i_u, i_v)
   call accumulate_cross_product_derivative_field(field_name_index("v_cross_dyvv", velocity_field_names, "velocity"), i_v, i_v, i_v)
   call accumulate_cross_product_derivative_field(field_name_index("w_cross_dyvw", velocity_field_names, "velocity"), i_w, i_v, i_w)

    call accumulate_cross_pressure(field_name_index("u_cross_p", velocity_field_names, "velocity"), i_u, .false.)
    call accumulate_cross_pressure(field_name_index("v_cross_dpdy", velocity_field_names, "velocity"), i_v, .true.)
    call accumulate_cross_pressure(field_name_index("w_cross_p", velocity_field_names, "velocity"), i_w, .false.)
  end subroutine acc_convvelo_stats

  subroutine finish_convvelo_field(field_index)
    implicit none

    integer(C_INT), intent(in) :: field_index
    integer(C_INT) :: storage_index
    real(C_DOUBLE) :: old_weight, new_weight
    integer(C_INT) :: ix, iy, iz

    if (.not. convvelo_initialized) call init_convvelo()
    if (field_index < 1 .or. field_index > size(convvelo_field_map)) then
      error stop "finish_convvelo_field: field_index out of range"
    end if
    storage_index = convvelo_field_map(field_index)
    if (storage_index == 0) return

    n_field_samples(storage_index) = n_field_samples(storage_index) + 1_C_INT64_T
    old_weight = dble(n_field_samples(storage_index) - 1_C_INT64_T)/dble(n_field_samples(storage_index))
    new_weight = 1.0d0/dble(n_field_samples(storage_index))
    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(convvelo_stats, convvelo_work, storage_index, old_weight, new_weight, ny0, nyN, nz, nx0, nxN) private(ix, iz, iy)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0 - 2, nyN + 2
          convvelo_stats(iy, iz, ix, storage_index) = old_weight*convvelo_stats(iy, iz, ix, storage_index) + &
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

    call apply_complex_derivative_with_y_pencil(V(:, :, :, component_index), convvelo_work)
  end subroutine apply_dy_to_work

  subroutine apply_dyy_to_work(component_index)
    implicit none
    integer(C_INT), intent(in) :: component_index
    integer(C_INT) :: iy, iz, ix

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
        do iy = nyN, ny0, -1
          convvelo_work(iy, iz, ix) = convvelo_work(iy, iz, ix) - ( &
                                      D0mat(iy, 1)*convvelo_work(iy + 1, iz, ix) + &
                                      D0mat(iy, 2)*convvelo_work(iy + 2, iz, ix))
          convvelo_work(iy, iz, ix) = convvelo_work(iy, iz, ix)*D0mat(iy, 0)
        end do
        do iy = ny0, nyN + 2
          convvelo_work(iy, iz, ix) = convvelo_work(iy, iz, ix) - ( &
                                      D0mat(iy, -2)*convvelo_work(iy - 2, iz, ix) + &
                                      D0mat(iy, -1)*convvelo_work(iy - 1, iz, ix))
        end do
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
    if (fft_transpose_is_local) then
      call repack_zTOx_local(VVdz(:, :, :, 1), VVdx(:, :, :, 1), ny)
    else
      call pack_zTOx(VVdz(:, :, :, 1), sendbuf(:, 1), ny)
      call alltoall(sendbuf(:, 1), recvbuf(:, 1), request)
    end if
#ifdef HAVE_MPI
    if (.not. fft_transpose_is_local) call MPI_Wait(request, status, ierr)
#endif
    if (.not. fft_transpose_is_local) call unpack_zTOx(recvbuf(:, 1), VVdx(:, :, :, 1), ny)
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
    if (fft_transpose_is_local) then
      call repack_xTOz_local(VVdx(:, :, :, 1), VVdz(:, :, :, 1), ny)
    else
      call pack_xTOz(VVdx(:, :, :, 1), sendbuf(:, 1), ny)
      call alltoall(sendbuf(:, 1), recvbuf(:, 1), request)
    end if
#ifdef HAVE_MPI
    if (.not. fft_transpose_is_local) call MPI_Wait(request, status, ierr)
#endif
    if (.not. fft_transpose_is_local) call unpack_xTOz(recvbuf(:, 1), VVdz(:, :, :, 1), ny)
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

    if (.not. convvelo_field_requested(field_index)) return
    call load_component_to_work(rhs_component)
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_components

  subroutine accumulate_cross_derivative(field_index, lhs_component, rhs_component, derivative_order)
    implicit none
    integer(C_INT), intent(in) :: field_index, lhs_component, rhs_component, derivative_order

    if (.not. convvelo_field_requested(field_index)) return
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

    if (.not. convvelo_field_requested(field_index)) return
    if (use_dpdy) then
      call compute_dpdy(convvelo_work)
    else
      call compute_poisson(convvelo_work)
    end if
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_pressure

  subroutine accumulate_cross_product_field(field_index, lhs_component, rhs0, rhs1)
    implicit none
    integer(C_INT), intent(in) :: field_index, lhs_component, rhs0, rhs1

    if (.not. convvelo_field_requested(field_index)) return
    call build_cross_product_work(rhs0, rhs1)
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_product_field

  subroutine accumulate_cross_product_derivative_field(field_index, lhs_component, rhs0, rhs1)
    implicit none
    integer(C_INT), intent(in) :: field_index, lhs_component, rhs0, rhs1

    if (.not. convvelo_field_requested(field_index)) return
    call build_cross_product_work(rhs0, rhs1)
    call apply_dy_to_existing_work()
    call multiply_work_by_conjugate(lhs_component)
    call finish_convvelo_field(field_index)
  end subroutine accumulate_cross_product_derivative_field

  subroutine build_cross_product_work(rhs0, rhs1)
    implicit none
    integer(C_INT), intent(in) :: rhs0, rhs1
    integer(C_INT) :: iy, iz, ix

    call load_convvelo_field_to_zbuf(rhs0)
    call spectral_field_to_real_x(convvelo_real0)
    call load_convvelo_field_to_zbuf(rhs1)
    call spectral_field_to_real_x(convvelo_real1)
    !$omp target teams distribute parallel do collapse(3) &
    !$omp shared(convvelo_real_prod, convvelo_real0, convvelo_real1, factor, nxd, nzB, ny) private(ix, iz, iy)
    do iy = 1, ny + 3
      do iz = 1, nzB
        do ix = 1, 2*(nxd + 1)
          convvelo_real_prod(ix, iz, iy) = factor*convvelo_real0(ix, iz, iy)*convvelo_real1(ix, iz, iy)
        end do
      end do
    end do
    call real_x_to_spectral_field(convvelo_real_prod, convvelo_work)
  end subroutine build_cross_product_work

  subroutine apply_dy_to_existing_work()
    implicit none
    complex(C_DOUBLE_COMPLEX), allocatable :: deriv(:, :, :)

    allocate (deriv(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    call apply_complex_derivative_with_y_pencil(convvelo_work, deriv)
    convvelo_work = deriv
    deallocate (deriv)
  end subroutine apply_dy_to_existing_work

  subroutine free_convvelo()
    implicit none

    if (.not. convvelo_initialized) return

    !$omp target exit data map(delete: convvelo_stats, convvelo_work, component_means, convvelo_real0, convvelo_real1, convvelo_real_prod)
    deallocate (convvelo_stats, convvelo_work, component_means, convvelo_real0, convvelo_real1, convvelo_real_prod, n_field_samples)
    if (allocated(convvelo_velocity_field_ids)) deallocate (convvelo_velocity_field_ids)
    if (allocated(convvelo_scalar_field_ids)) deallocate (convvelo_scalar_field_ids)
    if (allocated(convvelo_field_map)) deallocate (convvelo_field_map)

    n_convvelo_fields = 0
    n_convvelo_velocity_fields = 0
    n_convvelo_scalar_fields = 0
    n_convvelo_profile_slots = 0_C_INT64_T
    convvelo_last_write_index = -1_C_INT64_T
    n_mean_samples = 0_C_INT64_T
    convvelo_average_start_time = 0.0d0
    convvelo_average_end_time = 0.0d0
    convvelo_dirty = .false.
    convvelo_initialized = .false.
  end subroutine free_convvelo

  logical function convvelo_has_pending_output()
    implicit none

    convvelo_has_pending_output = convvelo_initialized .and. convvelo_dirty .and. &
                                  (n_mean_samples > 0_C_INT64_T .or. any(n_field_samples > 0_C_INT64_T))
  end function convvelo_has_pending_output

  subroutine init_convvelo_runtime()
    implicit none

    if (.not. convvelo_enabled) return

    call init_convvelo()
    call reset_convvelo_stats()
    if (convvelo_dt_write > 0.0d0) then
      convvelo_last_write_index = int(floor((time + 0.5d0*deltat)/convvelo_dt_write), C_INT64_T)
    else
      convvelo_last_write_index = -1_C_INT64_T
    end if
    call write_convvelo_field_layout(convvelo_runtime_filename)
  end subroutine init_convvelo_runtime

  subroutine advance_convvelo_runtime()
    implicit none

    if (.not. convvelo_enabled) return

    if (crossed_convvelo_interval(convvelo_dt_compute, convvelo_t_start)) then
      if (has_terminal) write (*, *) "Computing convvelo stats at time ", time
      call update_convvelo_component_means()
      call acc_convvelo_stats()
    end if

    if (convvelo_dt_write > 0.0d0) then
      if (convvelo_has_pending_output() .and. crossed_convvelo_interval(convvelo_dt_write, convvelo_t_start)) then
        if (has_terminal) write (*, *) "Writing convvelo snapshot at time ", time
        call write_convvelo_runtime_snapshot()
      end if
    end if
  end subroutine advance_convvelo_runtime

  subroutine finalize_convvelo_runtime()
    implicit none

    if (.not. convvelo_enabled) return

    if (convvelo_has_pending_output()) then
      call write_convvelo_runtime_snapshot()
    end if
    call free_convvelo()
  end subroutine finalize_convvelo_runtime

  subroutine write_convvelo_runtime_snapshot()
    implicit none
    character(len=256) :: snapshot_filename
    character(len=32) :: index_string
    logical :: exists

    convvelo_last_write_index = 0_C_INT64_T
    do
      write (index_string, '(I0)') convvelo_last_write_index
      snapshot_filename = "convvelo."//trim(index_string)//".bin"
      inquire (file=trim(snapshot_filename), exist=exists)
      if (.not. exists) exit
      convvelo_last_write_index = convvelo_last_write_index + 1_C_INT64_T
    end do

    if (has_terminal) write (*, *) "Writing "//trim(snapshot_filename)//" at time ", time
    call write_convvelo_raw_stats(snapshot_filename)
    convvelo_dirty = .false.
  end subroutine write_convvelo_runtime_snapshot

  subroutine write_convvelo_raw_stats(filename)
    implicit none

    character(len=*), intent(in) :: filename
    integer(C_INT) :: field_index, iPhi
    real(C_DOUBLE) :: header_times(2)
    integer(C_INT64_T) :: header_sample_count

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
    header_times = [convvelo_average_start_time, convvelo_average_end_time]
    header_sample_count = n_mean_samples

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
    total_bytes = convvelo_file_header_bytes + int(n_convvelo_profile_slots, MPI_OFFSET_KIND)*profile_bytes + &
                  int(n_convvelo_fields, MPI_OFFSET_KIND)*field_bytes

    call MPI_File_open(MPI_COMM_WORLD, trim(filename), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
    call MPI_File_set_size(fh, total_bytes)

    if (iproc == 0) then
      call MPI_File_write_at(fh, 0_MPI_OFFSET_KIND, header_times, 2, MPI_DOUBLE_PRECISION, status)
      call MPI_File_write_at(fh, 16_MPI_OFFSET_KIND, header_sample_count, 1, MPI_INTEGER8, status)
    end if

    disp = convvelo_file_header_bytes
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, profile_file_type, 'native', MPI_INFO_NULL)
    call MPI_File_write_all(fh, component_means(:, 1), 1, profile_mem_type, status)
    disp = convvelo_file_header_bytes + profile_bytes
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, profile_file_type, 'native', MPI_INFO_NULL)
    call MPI_File_write_all(fh, component_means(:, 2), 1, profile_mem_type, status)
    disp = convvelo_file_header_bytes + 2_MPI_OFFSET_KIND*profile_bytes
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, profile_file_type, 'native', MPI_INFO_NULL)
    call MPI_File_write_all(fh, component_means(:, 3), 1, profile_mem_type, status)
    do iPhi = 1, nPhi
      disp = convvelo_file_header_bytes + int(2 + iPhi, MPI_OFFSET_KIND)*profile_bytes
      call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, profile_file_type, 'native', MPI_INFO_NULL)
      call MPI_File_write_all(fh, component_means(:, 3 + iPhi), 1, profile_mem_type, status)
    end do

    do field_index = 1, n_convvelo_fields
      disp = convvelo_file_header_bytes + int(n_convvelo_profile_slots, MPI_OFFSET_KIND)*profile_bytes + &
             int(field_index - 1, MPI_OFFSET_KIND)*field_bytes
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

    write (99) header_times
    write (99) header_sample_count
    write (99) component_means(:, 1)
    write (99) component_means(:, 2)
    write (99) component_means(:, 3)
    do iPhi = 1, nPhi
      write (99) component_means(:, 3 + iPhi)
    end do
    do field_index = 1, n_convvelo_fields
      write (99) convvelo_stats(:, :, :, field_index)
    end do
    close (99)
#endif
  end subroutine write_convvelo_raw_stats

  subroutine write_convvelo_field_layout(filename)
    implicit none

    character(len=*), intent(in) :: filename
    character(len=512) :: layout_filename
    integer :: io, i

    if (iproc /= 0) return

    layout_filename = trim(filename)//".fields"
    open (unit=98, file=trim(layout_filename), status='replace', action='write', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open convvelo field layout file: ', trim(layout_filename)
      stop 1
    end if

    write (98, '(A)', advance='no') 'profile_fields:'
    do i = 1, size(convvelo_profile_field_names)
      write (98, '(1X,A)', advance='no') trim(convvelo_profile_field_names(i))
    end do
    write (98, *)

    write (98, '(A)', advance='no') 'velocity_fields:'
    do i = 1, size(convvelo_velocity_field_ids)
      write (98, '(1X,A)', advance='no') trim(velocity_field_names(convvelo_velocity_field_ids(i)))
    end do
    write (98, *)

    write (98, '(A)', advance='no') 'scalar_fields:'
    do i = 1, size(convvelo_scalar_field_ids)
      write (98, '(1X,A)', advance='no') trim(scalar_field_names(convvelo_scalar_field_ids(i)))
    end do
    write (98, *)
    close (98)
  end subroutine write_convvelo_field_layout

  integer(C_INT) function raw_scalar_field_index(i_phi, name)
    implicit none
    integer(C_INT), intent(in) :: i_phi
    character(len=*), intent(in) :: name

    raw_scalar_field_index = n_convvelo_velocity_fields_total + (i_phi - 1)*n_convvelo_scalar_fields_total + &
                             field_name_index(name, scalar_field_names, "scalar")
  end function raw_scalar_field_index

  integer(C_INT) function field_name_index(name, names, category)
    implicit none
    character(len=*), intent(in) :: name, category
    character(len=*), intent(in) :: names(:)
    integer(C_INT) :: i

    do i = 1, size(names)
      if (trim(names(i)) == trim(name)) then
        field_name_index = i
        return
      end if
    end do
    error stop "field_name_index: unknown "//trim(category)//" field name"
  end function field_name_index

  logical function convvelo_field_requested(field_index)
    implicit none
    integer(C_INT), intent(in) :: field_index

    if (.not. allocated(convvelo_field_map)) then
      convvelo_field_requested = .false.
      return
    end if
 convvelo_field_requested = field_index >= 1 .and. field_index <= size(convvelo_field_map) .and. convvelo_field_map(field_index) > 0
  end function convvelo_field_requested

  logical function crossed_convvelo_interval(period, t_start)
    implicit none
    real(C_DOUBLE), intent(in) :: period, t_start

    crossed_convvelo_interval = .false.
    if (period <= 0.0d0) return
    if (time + 0.5d0*deltat < t_start) return

    crossed_convvelo_interval = floor((time + 0.5d0*deltat - t_start)/period) > &
                                floor((time - 0.5d0*deltat - t_start)/period)
  end function crossed_convvelo_interval

end module convvelo
