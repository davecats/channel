program test_convvelo_mpi_io
  use convvelo, only: init_convvelo, fill_convvelo_synthetic_state_for_test, write_convvelo_raw_stats, free_convvelo
  use dnsdata, only: free_memory
  use driver, only: initialize
  use pressure_output, only: free_pressure_output
  use mpi_f08
  implicit none
  character(len=256) :: config_file, restart_file

  call get_command_argument(1, config_file)
  call get_command_argument(2, restart_file)
  if (len_trim(config_file) == 0 .or. len_trim(restart_file) == 0) then
    write (*, *) "Usage: test_convvelo_mpi_io <config_file> <restart_file>"
    stop 2
  end if

  call initialize(trim(config_file), trim(restart_file), .false.)
  call init_convvelo()
  call fill_convvelo_synthetic_state_for_test()
  call write_convvelo_raw_stats("convvelo_synthetic.bin")

  call free_convvelo()
  call free_pressure_output()
  call free_memory(.true.)
  call MPI_Finalize()
end program test_convvelo_mpi_io
