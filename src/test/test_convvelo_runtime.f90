program test_convvelo_runtime
  use, intrinsic :: iso_c_binding
  use case_setup, only: nPhi, ny0, nyN, nz, nx0, nxN, time, free_memory
  use convvelo, only: free_convvelo, convvelo_has_pending_output, write_convvelo_runtime_snapshot, convvelo_enabled, &
                      update_convvelo_component_means, acc_convvelo_stats
  use pressure_output, only: free_pressure_output
  use driver, only: initialize, timeloop
  use test_convvelo_utils
  use mpi_f08
  implicit none

  character(len=32) :: mode
  character(len=256) :: config_file, generated_file, reference_file
  character(len=128) :: full_name
  integer(C_INT64_T) :: average_count
  integer :: command_status
  integer(C_INT) :: i_field, i_phi, nfail, raw_field_index, minimal_scalar_index
  real(C_DOUBLE), parameter :: expected_start_time = 38.86540781764443d0
  real(C_DOUBLE), parameter :: expected_end_time = 40.99369209898264d0
  real(C_DOUBLE) :: average_end_time, average_start_time, field_tol, header_tol, mean_tol
  real(C_DOUBLE) :: reset_end_time, reset_start_time
  complex(C_DOUBLE_COMPLEX), allocatable :: generated(:, :, :)

  call get_command_argument(1, mode)
  if (len_trim(mode) == 0) mode = "full"

  field_tol = 2.0d-12
  header_tol = 1.0d-12
  mean_tol = 1.0d-12
  nfail = 0

  select case (trim(mode))
  case ("full")
    config_file = "tests/convvelo/dns_runtime.in"
    generated_file = "convvelo.0.bin"
    reference_file = "tests/convvelo/raw_statistics.bin"
  case ("minimal")
    config_file = "tests/convvelo/dns_runtime_minimal.in"
    generated_file = "convvelo.0.bin"
    reference_file = "tests/convvelo/convvelo_runtime_minimal.bin"
  case default
    write (*, *) "Unknown convvelo runtime mode: ", trim(mode)
    stop 2
  end select

  call execute_command_line("rm -f convvelo.bin convvelo.*.bin", exitstat=command_status)
  if (command_status /= 0) then
    write (*, *) "Could not remove old convvelo runtime files"
    stop 3
  end if

  call initialize(trim(config_file), "tests/data/start_field_scalar.out")
  allocate (generated(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
  call timeloop()
  if (convvelo_enabled .and. convvelo_has_pending_output()) then
    call write_convvelo_runtime_snapshot()
  end if

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
  if (convvelo_enabled) then
    call update_convvelo_component_means()
    call acc_convvelo_stats()
    call write_convvelo_runtime_snapshot()
    call read_convvelo_header("convvelo.1.bin", reset_start_time, reset_end_time, average_count)
    if (abs(reset_start_time - time) >= header_tol) then
      write (*, '(A,ES20.12,A,ES20.12)') "convvelo reset start_time mismatch: got ", reset_start_time, ", expected ", time
      nfail = nfail + 1
    end if
    if (abs(reset_end_time - time) >= header_tol) then
      write (*, '(A,ES20.12,A,ES20.12)') "convvelo reset end_time mismatch: got ", reset_end_time, ", expected ", time
      nfail = nfail + 1
    end if
    if (average_count /= 1_C_INT64_T) then
      write (*, '(A,I0,A,I0)') "convvelo reset average_count mismatch: got ", average_count, ", expected ", 1_C_INT64_T
      nfail = nfail + 1
    end if
  end if

  if (trim(mode) == "full") then
    call compare_profile_to_reference("mean_u", trim(generated_file), trim(reference_file), 0_C_INT, 0_C_INT, mean_tol, nfail)
    call compare_profile_to_reference("mean_v", trim(generated_file), trim(reference_file), 1_C_INT, 1_C_INT, mean_tol, nfail)
    call compare_profile_to_reference("mean_w", trim(generated_file), trim(reference_file), 2_C_INT, 2_C_INT, mean_tol, nfail)
    do i_phi = 1, nPhi
      write (full_name, '(A,"[phi=",I0,"]")') "mean_theta", i_phi
      call compare_profile_to_reference(trim(full_name), trim(generated_file), trim(reference_file), 2_C_INT + i_phi, 2_C_INT + i_phi, mean_tol, nfail)
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
    call finish_convvelo_test("Integrated convvelo regression PASSED", "Integrated convvelo regression FAILED, failing fields = ", nfail)
  else
    call compare_profile_to_reference("mean_v", trim(generated_file), trim(reference_file), 1_C_INT, 1_C_INT, mean_tol, nfail)
    call compare_profile_to_reference("mean_w", trim(generated_file), trim(reference_file), 2_C_INT, 2_C_INT, mean_tol, nfail)
    call compare_profile_to_reference("mean_u", trim(generated_file), trim(reference_file), 0_C_INT, 0_C_INT, mean_tol, nfail)
    do i_phi = 1, nPhi
      write (full_name, '(A,"[phi=",I0,"]")') "mean_theta", i_phi
      call compare_profile_to_reference(trim(full_name), trim(generated_file), trim(reference_file), 2_C_INT + i_phi, 2_C_INT + i_phi, mean_tol, nfail)
    end do
    do i_field = 1, n_convvelo_velocity_fields_minimal
      call read_minimal_field(trim(generated_file), i_field - 1, generated)
      call compare_field_to_reference( &
        trim(minimal_velocity_field_names(i_field)), generated, velocity_field_index(trim(minimal_velocity_field_names(i_field))) - 1, field_tol, nfail &
        )
    end do
    do i_phi = 1, nPhi
      do i_field = 1, n_convvelo_scalar_fields_minimal
        write (full_name, '(A,"[phi=",I0,"]")') trim(minimal_scalar_field_names(i_field)), i_phi
        minimal_scalar_index = n_convvelo_velocity_fields_minimal + (i_phi - 1)*n_convvelo_scalar_fields_minimal + (i_field - 1)
        raw_field_index = n_convvelo_velocity_fields + (i_phi - 1)*n_convvelo_scalar_fields + &
                          scalar_field_index(trim(minimal_scalar_field_names(i_field))) - 1
        call read_minimal_field(trim(generated_file), minimal_scalar_index, generated)
        call compare_field_to_reference(trim(full_name), generated, raw_field_index, field_tol, nfail)
      end do
    end do
    call finish_convvelo_test("Integrated convvelo minimal regression PASSED", "Integrated convvelo minimal regression FAILED, failing fields = ", nfail)
  end if

  deallocate (generated)
  if (convvelo_enabled) call free_convvelo()
  call free_pressure_output()
  call free_memory(.true.)
  call MPI_Finalize()

end program test_convvelo_runtime
