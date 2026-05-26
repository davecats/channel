program test_convvelo_runtime
  use, intrinsic :: iso_c_binding
  use dnsdata
  use convvelo, only: free_convvelo, convvelo_has_pending_output, write_convvelo_raw_stats
  use pressure_output
  use driver
  use test_convvelo_utils
  implicit none

  character(len=32) :: mode
  character(len=256) :: config_file, generated_file
  character(len=128) :: full_name
  integer(C_INT) :: i_field, i_phi, nfail, raw_field_index, minimal_scalar_index
  real(C_DOUBLE) :: field_tol, mean_tol
  complex(C_DOUBLE_COMPLEX), allocatable :: mean_profiles(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: generated(:, :, :)

  call get_command_argument(1, mode)
  if (len_trim(mode) == 0) mode = "full"

  field_tol = 2.0d-12
  mean_tol = 1.0d-12
  nfail = 0

  select case (trim(mode))
  case ("full")
    config_file = "tests/convvelo/dns_runtime.in"
    generated_file = "tests/convvelo/raw_statistics.bin"
  case ("minimal")
    config_file = "tests/convvelo/dns_runtime_minimal.in"
    generated_file = "tests/convvelo/convvelo_runtime_minimal.bin"
  case default
    write (*, *) "Unknown convvelo runtime mode: ", trim(mode)
    stop 2
  end select

  call initialize(trim(config_file), "tests/data/start_field_scalar.out")
  allocate (mean_profiles(ny0 - 2:nyN + 2, 3 + nPhi))
  allocate (generated(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
  call timeloop()
  if (convvelo_enabled .and. convvelo_has_pending_output()) then
    call write_convvelo_raw_stats(trim(convvelo_output_file))
  end if

  call compute_reference_mean_profiles(mean_profiles)
  if (trim(mode) == "full") then
    call compare_profile_to_expected("mean_u", trim(generated_file), 0_C_INT, mean_profiles(:, 1), mean_tol, nfail)
    call compare_profile_to_expected("mean_v", trim(generated_file), 1_C_INT, mean_profiles(:, 2), mean_tol, nfail)
    call compare_profile_to_expected("mean_w", trim(generated_file), 2_C_INT, mean_profiles(:, 3), mean_tol, nfail)
    do i_phi = 1, nPhi
      write (full_name, '(A,"[phi=",I0,"]")') "mean_theta", i_phi
      call compare_profile_to_expected(trim(full_name), trim(generated_file), 2_C_INT + i_phi, mean_profiles(:, 3 + i_phi), mean_tol, nfail)
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
    call compare_profile_to_expected("mean_v", trim(generated_file), 1_C_INT, mean_profiles(:, 2), mean_tol, nfail)
    call compare_profile_to_expected("mean_w", trim(generated_file), 2_C_INT, mean_profiles(:, 3), mean_tol, nfail)
    call compare_profile_to_expected("mean_u", trim(generated_file), 0_C_INT, mean_profiles(:, 1), mean_tol, nfail)
    do i_phi = 1, nPhi
      write (full_name, '(A,"[phi=",I0,"]")') "mean_theta", i_phi
      call compare_profile_to_expected(trim(full_name), trim(generated_file), 2_C_INT + i_phi, mean_profiles(:, 3 + i_phi), mean_tol, nfail)
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

  deallocate (mean_profiles, generated)
  if (convvelo_enabled) call free_convvelo()
  call free_pressure_output()
  call free_memory(.true.)
#ifdef HAVE_MPI
  call MPI_Finalize()
#endif

end program test_convvelo_runtime
