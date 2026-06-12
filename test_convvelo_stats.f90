program test_convvelo_stats
  use, intrinsic :: iso_c_binding
  use dnsdata, only: nPhi, ny0, nyN, nz, nx0, nxN, V, read_restart_file, sync_velocity_to_device, free_memory
  use convvelo, only: init_convvelo, reset_convvelo_stats, acc_convvelo_stats, convvelo_stats, free_convvelo
  use pressure_output, only: free_pressure_output
  use driver, only: initialize
  use test_convvelo_utils
#ifdef HAVE_MPI
  use mpi_f08
#endif
  implicit none

  character(len=*), parameter :: config_file = "tests/convvelo/dns.in"
  character(len=*), parameter :: restart_in = "tests/convvelo/Dati.cart.out"

  integer(C_INT) :: i_field, i_phi, nfail, raw_field_index
  real(C_DOUBLE) :: tol
  character(len=128) :: full_name

  tol = 1.0d-12
  nfail = 0

  call initialize(config_file, restart_in)
  call init_convvelo()
  call reset_convvelo_stats()
  call read_restart_file("tests/convvelo/Dati.cart.39.out", V)
  call sync_velocity_to_device()
  call acc_convvelo_stats()
  call read_restart_file("tests/convvelo/Dati.cart.40.out", V)
  call sync_velocity_to_device()
  call acc_convvelo_stats()
  call read_restart_file("tests/convvelo/Dati.cart.41.out", V)
  call sync_velocity_to_device()
  call acc_convvelo_stats()
  !$omp target update from(convvelo_stats)

  do i_field = 1, n_convvelo_velocity_fields
    call compare_field_to_reference(trim(velocity_field_names(i_field)), convvelo_stats(:, :, :, i_field), i_field - 1, tol, nfail)
  end do

  do i_phi = 1, nPhi
    do i_field = 1, n_convvelo_scalar_fields
      write (full_name, '(A,"[phi=",I0,"]")') trim(scalar_field_names(i_field)), i_phi
      raw_field_index = n_convvelo_velocity_fields + (i_phi - 1)*n_convvelo_scalar_fields + (i_field - 1)
      call compare_field_to_reference(trim(full_name), convvelo_stats(:, :, :, raw_field_index + 1), raw_field_index, tol, nfail)
    end do
  end do

  call finish_convvelo_test("Convvelo statistics regression PASSED", "Convvelo statistics regression FAILED, failing fields = ", nfail)

  call free_pressure_output()
  call free_convvelo()
  call free_memory(.true.)
#ifdef HAVE_MPI
  call MPI_Finalize()
#endif

end program test_convvelo_stats
