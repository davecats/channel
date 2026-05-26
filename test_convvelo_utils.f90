module test_convvelo_utils
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use dnsdata
#ifdef HAVE_MPI
  use mpi_f08
#endif
  implicit none

  private

  character(len=*), parameter, public :: convvelo_reference_file = "tests/convvelo/raw_statistics.bin"
  character(len=*), parameter, public :: convvelo_reference_snapshots(3) = [character(len=32) :: &
                           "tests/convvelo/Dati.cart.39.out", "tests/convvelo/Dati.cart.40.out", "tests/convvelo/Dati.cart.41.out" &
                                                                            ]
  integer(C_INT), parameter, public :: n_convvelo_velocity_fields = 33
  integer(C_INT), parameter, public :: n_convvelo_scalar_fields = 10
  integer(C_INT), parameter, public :: n_convvelo_velocity_fields_minimal = 20
  integer(C_INT), parameter, public :: n_convvelo_scalar_fields_minimal = 9

  character(len=32), parameter, public :: velocity_field_names(n_convvelo_velocity_fields) = [character(len=32) :: &
                               "u_cross_u", "u_cross_dyu", "u_cross_v", "u_cross_dyv", "u_cross_w", "u_cross_dyw", "u_cross_dyyu", &
                               "v_cross_u", "v_cross_dyu", "v_cross_v", "v_cross_dyv", "v_cross_w", "v_cross_dyw", "v_cross_dyyv", &
                               "w_cross_u", "w_cross_dyu", "w_cross_v", "w_cross_dyv", "w_cross_w", "w_cross_dyw", "w_cross_dyyw", &
                                 "u_cross_p", "v_cross_dpdy", "w_cross_p", "u_cross_uu", "u_cross_uw", "v_cross_uv", "v_cross_vw", &
                                                        "w_cross_uw", "w_cross_ww", "u_cross_dyuv", "v_cross_dyvv", "w_cross_dyvw" &
                                                                                              ]
  character(len=32), parameter, public :: scalar_field_names(n_convvelo_scalar_fields) = [character(len=32) :: &
                                                                 "t_cross_t", "t_cross_u", "t_cross_v", "t_cross_w", "t_cross_tu", &
                                                        "t_cross_tw", "t_cross_dyyt", "t_cross_dytv", "t_cross_dyt", "t_cross_dyv" &
                                                                                          ]

  ! Minimal mode stores only the reduced zero-crossflow subset.
  character(len=32), parameter, public :: minimal_velocity_field_names(n_convvelo_velocity_fields_minimal) = [character(len=32) :: &
                                                                               "u_cross_u", "u_cross_v", "v_cross_v", "w_cross_w", &
                                                                                         "u_cross_p", "v_cross_dpdy", "w_cross_p", &
                                               "u_cross_uu", "u_cross_uw", "v_cross_uv", "v_cross_vw", "w_cross_uw", "w_cross_ww", &
                                                                                   "u_cross_dyyu", "v_cross_dyyv", "w_cross_dyyw", &
                                                                     "u_cross_dyuv", "v_cross_dyvv", "w_cross_dyvw", "u_cross_dyv" &
                                                                                                              ]
  character(len=32), parameter, public :: minimal_scalar_field_names(n_convvelo_scalar_fields_minimal) = [character(len=32) :: &
                                                                 "t_cross_t", "t_cross_u", "t_cross_v", "t_cross_w", "t_cross_tu", &
                                                                       "t_cross_tw", "t_cross_dyyt", "t_cross_dytv", "t_cross_dyv" &
                                                                                                          ]
  integer(C_INT), parameter :: n_convvelo_profile_header_slots = 3

  public :: compare_field_to_reference
  public :: compare_profile_to_reference
  public :: read_full_field
  public :: read_full_profile
  public :: read_minimal_field
  public :: compute_reference_mean_profiles
  public :: finish_convvelo_test
  public :: read_raw_stat_field_mpi
  public :: velocity_field_index
  public :: scalar_field_index

contains

  subroutine finish_convvelo_test(success_message, failure_message, nfail)
    implicit none
    character(len=*), intent(in) :: success_message, failure_message
    integer(C_INT), intent(in) :: nfail
    integer(C_INT) :: global_nfail

    global_nfail = nfail
#ifdef HAVE_MPI
    call MPI_Allreduce(nfail, global_nfail, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    if (iproc == 0) then
      if (global_nfail == 0) then
        write (*, *) trim(success_message)
      else
        write (*, *) trim(failure_message), global_nfail
      end if
    end if

#ifdef HAVE_MPI
    if (global_nfail /= 0) call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  end subroutine finish_convvelo_test

  subroutine compute_reference_mean_profiles(mean_profiles)
    implicit none
    complex(C_DOUBLE_COMPLEX), intent(out) :: mean_profiles(ny0 - 2:nyN + 2, 3 + nPhi)
    integer(C_INT) :: i_snapshot, i_phi

    mean_profiles = (0.0d0, 0.0d0)
    do i_snapshot = 1, size(convvelo_reference_snapshots)
      call read_restart_file(convvelo_reference_snapshots(i_snapshot), V)
      mean_profiles(:, 1) = mean_profiles(:, 1) + V(:, 0, 0, 1)
      mean_profiles(:, 2) = mean_profiles(:, 2) + V(:, 0, 0, 2)
      mean_profiles(:, 3) = mean_profiles(:, 3) + V(:, 0, 0, 3)
      do i_phi = 1, nPhi
        mean_profiles(:, 3 + i_phi) = mean_profiles(:, 3 + i_phi) + V(:, 0, 0, 3 + i_phi)
      end do
    end do
    mean_profiles = mean_profiles/dble(size(convvelo_reference_snapshots))
  end subroutine compute_reference_mean_profiles

  subroutine compare_field_to_reference(name, field, raw_field_index, tol_local, nfail_local)
    implicit none
    character(len=*), intent(in) :: name
    complex(C_DOUBLE_COMPLEX), intent(in) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer(C_INT), intent(in) :: raw_field_index
    real(C_DOUBLE), intent(in) :: tol_local
    integer(C_INT), intent(inout) :: nfail_local

    complex(C_DOUBLE_COMPLEX), allocatable :: reference(:, :, :)
    allocate (reference(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    call read_raw_stat_field_mpi(convvelo_reference_file, raw_field_index, reference)
    call compare_complex_fields(name, field, reference, tol_local, nfail_local)
    deallocate (reference)
  end subroutine compare_field_to_reference

  subroutine compare_profile_to_reference(name, generated_file, reference_file, generated_index, reference_index, tol_local, nfail_local)
    implicit none
    character(len=*), intent(in) :: name, generated_file, reference_file
    integer(C_INT), intent(in) :: generated_index, reference_index
    real(C_DOUBLE), intent(in) :: tol_local
    integer(C_INT), intent(inout) :: nfail_local

    complex(C_DOUBLE_COMPLEX), allocatable :: generated(:, :, :), reference(:, :, :)
    allocate (generated(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    allocate (reference(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    generated = (0.0d0, 0.0d0)
    reference = (0.0d0, 0.0d0)
    call read_profile(trim(generated_file), generated_index, generated(:, 0, 0))
    call read_profile(trim(reference_file), reference_index, reference(:, 0, 0))
    call compare_complex_fields(name, generated, reference, tol_local, nfail_local)
    deallocate (generated, reference)
  end subroutine compare_profile_to_reference

  integer(C_INT) function velocity_field_index(name)
    implicit none
    character(len=*), intent(in) :: name

    velocity_field_index = find_name_index(name, velocity_field_names)
  end function velocity_field_index

  integer(C_INT) function scalar_field_index(name)
    implicit none
    character(len=*), intent(in) :: name

    scalar_field_index = find_name_index(name, scalar_field_names)
  end function scalar_field_index

  integer(C_INT) function find_name_index(name, names)
    implicit none
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: names(:)
    integer(C_INT) :: i

    do i = 1, size(names)
      if (trim(name) == trim(names(i))) then
        find_name_index = i
        return
      end if
    end do

    write (*, *) "Unknown convvelo field name: ", trim(name)
    stop 2
  end function find_name_index

  subroutine compare_complex_fields(name, computed, reference, tol_local, nfail_local)
    implicit none
    character(len=*), intent(in) :: name
    complex(C_DOUBLE_COMPLEX), intent(in) :: computed(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    complex(C_DOUBLE_COMPLEX), intent(in) :: reference(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    real(C_DOUBLE), intent(in) :: tol_local
    integer(C_INT), intent(inout) :: nfail_local

    real(C_DOUBLE) :: local_max_err, global_max_err, diff
    logical :: local_has_nan, global_has_nan
    integer(C_INT) :: iy, iz, ix
    integer(C_INT) :: worst_iy, worst_iz, worst_ix
    complex(C_DOUBLE_COMPLEX) :: worst_computed, worst_reference

    local_max_err = 0.0d0
    local_has_nan = .false.
    worst_iy = ny0 - 2
    worst_iz = -nz
    worst_ix = nx0
    worst_computed = (0.0d0, 0.0d0)
    worst_reference = (0.0d0, 0.0d0)

    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0 - 2, nyN + 2
          if (ieee_is_nan(real(computed(iy, iz, ix))) .or. ieee_is_nan(aimag(computed(iy, iz, ix)))) local_has_nan = .true.
          diff = abs(computed(iy, iz, ix) - reference(iy, iz, ix))
          if (diff > local_max_err) then
            local_max_err = diff
            worst_iy = iy
            worst_iz = iz
            worst_ix = ix
            worst_computed = computed(iy, iz, ix)
            worst_reference = reference(iy, iz, ix)
          end if
        end do
      end do
    end do

#ifdef HAVE_MPI
    call MPI_Allreduce(local_max_err, global_max_err, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(local_has_nan, global_has_nan, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#else
    global_max_err = local_max_err
    global_has_nan = local_has_nan
#endif

    if (iproc == 0) write (*, '(A,": max_err=",ES12.4,", has_nan=",L1)') trim(name), global_max_err, global_has_nan

    if (global_has_nan .or. global_max_err >= tol_local) then
      nfail_local = nfail_local + 1
#ifdef HAVE_MPI
      if (abs(local_max_err - global_max_err) < 1.0d-15) then
        write (*, '(A,": worst rank/iy/iz/ix = ",4(I0,1X))') trim(name), iproc, worst_iy, worst_iz, worst_ix
        write (*, '(A,": computed / reference = ",2("(",ES20.12,",",ES20.12,")",1X))') &
          trim(name), worst_computed, worst_reference
      end if
#else
      write (*, '(A,": worst iy/iz/ix = ",3(I0,1X))') trim(name), worst_iy, worst_iz, worst_ix
      write (*, '(A,": computed / reference = ",2("(",ES20.12,",",ES20.12,")",1X))') &
        trim(name), worst_computed, worst_reference
#endif
    end if
  end subroutine compare_complex_fields

  subroutine read_profile(filename, profile_index, profile)
    implicit none
    character(len=*), intent(in) :: filename
    integer(C_INT), intent(in) :: profile_index
    complex(C_DOUBLE_COMPLEX), intent(out) :: profile(ny0 - 2:nyN + 2)
    integer :: io
    integer(C_INT64_T) :: profile_bytes, pos_bytes

    profile_bytes = int(16, C_INT64_T)*int(ny + 3, C_INT64_T)
    pos_bytes = int(profile_index, C_INT64_T)*profile_bytes + 1_C_INT64_T

    open (unit=97, file=filename, form='unformatted', access='stream', status='old', action='read', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open convvelo profile file: ', trim(filename)
      stop 1
    end if
    read (97, pos=pos_bytes) profile
    close (97)
  end subroutine read_profile

  subroutine read_full_profile(filename, profile_index, profile)
    implicit none
    character(len=*), intent(in) :: filename
    integer(C_INT), intent(in) :: profile_index
    complex(C_DOUBLE_COMPLEX), intent(out) :: profile(ny0 - 2:nyN + 2)

    call read_profile(filename, profile_index, profile)
  end subroutine read_full_profile

  subroutine read_full_field(filename, field_index, ref)
    implicit none
    character(len=*), intent(in) :: filename
    integer(C_INT), intent(in) :: field_index
    complex(C_DOUBLE_COMPLEX), intent(out) :: ref(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer :: io
    integer(C_INT64_T) :: profile_bytes, field_bytes, pos_bytes

    profile_bytes = int(16, C_INT64_T)*int(ny + 3, C_INT64_T)
    field_bytes = int(16, C_INT64_T)*int(ny + 3, C_INT64_T)*int(2*nz + 1, C_INT64_T)*int(nx + 1, C_INT64_T)
    pos_bytes = int(n_convvelo_profile_header_slots + nPhi, C_INT64_T)*profile_bytes + int(field_index, C_INT64_T)*field_bytes + 1_C_INT64_T

    open (unit=94, file=filename, form='unformatted', access='stream', status='old', action='read', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open convvelo full file: ', trim(filename)
      stop 1
    end if
    read (94, pos=pos_bytes) ref
    close (94)
  end subroutine read_full_field

  subroutine read_minimal_field(filename, field_index, ref)
    implicit none
    character(len=*), intent(in) :: filename
    integer(C_INT), intent(in) :: field_index
    complex(C_DOUBLE_COMPLEX), intent(out) :: ref(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer :: io
    integer(C_INT64_T) :: profile_bytes, field_bytes, pos_bytes

    profile_bytes = int(16, C_INT64_T)*int(ny + 3, C_INT64_T)
    field_bytes = int(16, C_INT64_T)*int(ny + 3, C_INT64_T)*int(2*nz + 1, C_INT64_T)*int(nx + 1, C_INT64_T)
    pos_bytes = int(n_convvelo_profile_header_slots + nPhi, C_INT64_T)*profile_bytes + int(field_index, C_INT64_T)*field_bytes + 1_C_INT64_T

    open (unit=96, file=filename, form='unformatted', access='stream', status='old', action='read', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open convvelo minimal file: ', trim(filename)
      stop 1
    end if
    read (96, pos=pos_bytes) ref
    close (96)
  end subroutine read_minimal_field

  subroutine read_raw_stat_field_mpi(filename, field_index, ref)
    use, intrinsic :: iso_c_binding
    implicit none
    character(len=*), intent(in) :: filename
    integer(C_INT), intent(in) :: field_index
    complex(C_DOUBLE_COMPLEX), intent(out) :: ref(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
#ifdef HAVE_MPI
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    type(MPI_Datatype) :: file_type, mem_type
    integer :: ierror
    integer, parameter :: ndims = 3
    integer :: sizes(ndims), subsizes(ndims), starts(ndims)
    integer(MPI_OFFSET_KIND) :: disp
    integer(MPI_OFFSET_KIND) :: field_bytes, profile_bytes

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

    profile_bytes = int(16, MPI_OFFSET_KIND)*int(ny + 3, MPI_OFFSET_KIND)
    field_bytes = int(16, MPI_OFFSET_KIND)*int(ny + 3, MPI_OFFSET_KIND)* &
                  int(2*nz + 1, MPI_OFFSET_KIND)*int(nx + 1, MPI_OFFSET_KIND)
    disp = int(n_convvelo_profile_header_slots + nPhi, MPI_OFFSET_KIND)*profile_bytes + &
           int(field_index, MPI_OFFSET_KIND)*field_bytes

    call MPI_File_open(MPI_COMM_WORLD, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh)
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, file_type, 'native', MPI_INFO_NULL)
    call MPI_File_read_all(fh, ref, 1, mem_type, status)
    call MPI_File_close(fh)
    call MPI_Type_free(file_type, ierror)
    call MPI_Type_free(mem_type, ierror)
#else
    integer :: io
    integer(C_INT64_T) :: profile_bytes, field_bytes, pos_bytes

    profile_bytes = int(16, C_INT64_T)*int(ny + 3, C_INT64_T)
    field_bytes = int(16, C_INT64_T)*int(ny + 3, C_INT64_T)*int(2*nz + 1, C_INT64_T)*int(nx + 1, C_INT64_T)
    pos_bytes = int(n_convvelo_profile_header_slots + nPhi, C_INT64_T)*profile_bytes + int(field_index, C_INT64_T)*field_bytes + 1_C_INT64_T

    open (unit=95, file=filename, form='unformatted', access='stream', status='old', action='read', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open convvelo reference file: ', trim(filename)
      stop 1
    end if
    read (95, pos=pos_bytes) ref
    close (95)
#endif
  end subroutine read_raw_stat_field_mpi

end module test_convvelo_utils
