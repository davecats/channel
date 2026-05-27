program test_post_convvelo
  use, intrinsic :: iso_c_binding
  use dnsdata
  use convvelo, only: free_convvelo
  use pressure_output
  use driver
  use test_convvelo_utils
  implicit none

  character(len=256) :: config_file, restart_file, generated_file
  character(len=128) :: full_name
  integer(C_INT64_T) :: average_count
  integer(C_INT) :: i_field, i_phi, nfail, raw_field_index
  real(C_DOUBLE), parameter :: expected_start_time = 38.86540781764443d0
  real(C_DOUBLE), parameter :: expected_end_time = 40.99369209898264d0
  real(C_DOUBLE) :: average_end_time, average_start_time, field_tol, header_tol, mean_tol
  complex(C_DOUBLE_COMPLEX), allocatable :: generated(:, :, :)

  call get_command_argument(1, config_file)
  call get_command_argument(2, restart_file)
  call get_command_argument(3, generated_file)
  if (len_trim(config_file) == 0 .or. len_trim(restart_file) == 0 .or. len_trim(generated_file) == 0) then
    write (*, *) "Usage: test_post_convvelo <config_file> <restart_file> <generated_file>"
    stop 2
  end if

  field_tol = 2.0d-12
  header_tol = 1.0d-12
  mean_tol = 1.0d-12
  nfail = 0

  call initialize(trim(config_file), trim(restart_file), .false.)
  allocate (generated(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))

  call read_convvelo_header(trim(generated_file), average_start_time, average_end_time, average_count)
  if (abs(average_start_time - expected_start_time) >= header_tol) then
    write (*, '(A,ES20.12,A,ES20.12)') "convvelo average_start_time mismatch: got ", average_start_time, ", expected ", expected_start_time
    nfail = nfail + 1
  end if
  if (abs(average_end_time - expected_end_time) >= header_tol) then
   write (*, '(A,ES20.12,A,ES20.12)') "convvelo average_end_time mismatch: got ", average_end_time, ", expected ", expected_end_time
    nfail = nfail + 1
  end if
  if (average_count /= 3_C_INT64_T) then
    write (*, '(A,I0,A,I0)') "convvelo average_count mismatch: got ", average_count, ", expected ", 3_C_INT64_T
    nfail = nfail + 1
  end if

  call compare_profile_to_reference("mean_u", trim(generated_file), convvelo_reference_file, 0_C_INT, 0_C_INT, mean_tol, nfail)
  call compare_profile_to_reference("mean_v", trim(generated_file), convvelo_reference_file, 1_C_INT, 1_C_INT, mean_tol, nfail)
  call compare_profile_to_reference("mean_w", trim(generated_file), convvelo_reference_file, 2_C_INT, 2_C_INT, mean_tol, nfail)
  do i_phi = 1, nPhi
    write (full_name, '(A,"[phi=",I0,"]")') "mean_theta", i_phi
    call compare_profile_to_reference(trim(full_name), trim(generated_file), convvelo_reference_file, 2_C_INT + i_phi, 2_C_INT + i_phi, mean_tol, nfail)
  end do

  do i_field = 1, n_convvelo_velocity_fields
    call read_full_field(trim(generated_file), i_field - 1, generated)
    call compare_field_to_reference(trim(velocity_field_names(i_field)), generated, i_field - 1, field_tol, nfail)
  end do
  do i_phi = 1, nPhi
    do i_field = 1, n_convvelo_scalar_fields
      write (full_name, '(A,"[phi=",I0,"]")') trim(scalar_field_names(i_field)), i_phi
      raw_field_index = n_convvelo_velocity_fields + (i_phi - 1)*n_convvelo_scalar_fields + (i_field - 1)
      call read_full_field(trim(generated_file), raw_field_index, generated)
      call compare_field_to_reference(trim(full_name), generated, raw_field_index, field_tol, nfail)
    end do
  end do

  call finish_convvelo_test("post_convvelo regression PASSED", "post_convvelo regression FAILED, failing fields = ", nfail)

  deallocate (generated)
  call free_convvelo()
  call free_pressure_output()
  call free_memory(.false.)
#ifdef HAVE_MPI
  call MPI_Finalize()
#endif

end program test_post_convvelo
