#include "header.h"

program post_convvelo
  use, intrinsic :: iso_c_binding
  use dnsdata, only: iproc, read_restart_file, free_memory, V
  use convvelo, only: convvelo_enabled, free_convvelo, reset_convvelo_stats, update_convvelo_component_means, &
                      acc_convvelo_stats, convvelo_has_pending_output, write_convvelo_runtime_snapshot
  use pressure_output, only: free_pressure_output
  use driver, only: initialize
#ifdef HAVE_FFTW
  use ffts, only: free_fft, VVdz, VVdx, rVVdx
#endif
#ifdef HAVE_MPI
  use mpi_f08
#endif
  implicit none

  character(len=*), parameter :: config_file = "dns.in"
  character(len=256), allocatable :: restart_files(:)
  integer :: nfiles, ifile

  call discover_restart_files(restart_files, nfiles)
  if (nfiles == 0) then
    write (*, *) "post_convvelo: No restart files matching Dati.cart.*.out found in the current directory."
    stop 1
  end if

  call initialize(config_file, trim(restart_files(1)), .false.)
  if (.not. convvelo_enabled) then
    write (*, *) "post_convvelo: convvelo output is not enabled in dns.in."
    stop 1
  end if

  call reset_convvelo_stats()

  do ifile = 1, nfiles
    if (ifile > 1) then
      if (iproc == 0) write (*, '(A,1X,A)') "Reading", trim(restart_files(ifile))
      call read_restart_file(trim(restart_files(ifile)), V)
      !$omp target update to(V)
    end if

    if (iproc == 0) write (*, '(A,1X,A)') "Processing", trim(restart_files(ifile))
    call update_convvelo_component_means()
    call acc_convvelo_stats()
    if (iproc == 0) write (*, '(A,1X,A)') "Processing finished", trim(restart_files(ifile))
  end do

  if (.not. convvelo_has_pending_output()) then
    write (*, *) "post_convvelo: No convvelo output was accumulated."
    stop 1
  end if
  call write_convvelo_runtime_snapshot()

  deallocate (restart_files)

  call free_convvelo()
  call free_pressure_output()
#ifdef HAVE_FFTW
  call free_fft(VVdz, VVdx, rVVdx)
#endif
  call free_memory(.false.)
#ifdef HAVE_MPI
  call MPI_Finalize()
#endif

contains

  subroutine discover_restart_files(files, nfiles)
    character(len=256), allocatable, intent(out) :: files(:)
    integer, intent(out) :: nfiles
    character(len=512) :: command, listfile, line
    integer :: unit, io, count, exitstat, cmdstat
    integer(C_INT) :: pid

    pid = getpid()
    write (listfile, '(A,I0)') ".post_convvelo_files.", pid
    command = "find . -maxdepth 1 -name 'Dati.cart.*.out' -printf '%f\n' | sort -V > "//trim(listfile)
    call execute_command_line(trim(command), wait=.true., exitstat=exitstat, cmdstat=cmdstat)
    if (cmdstat /= 0 .or. exitstat /= 0) then
      write (*, *) "post_convvelo: Failed to scan restart files in the current directory."
      stop 1
    end if

    open (newunit=unit, file=trim(listfile), status="old", action="read", iostat=io)
    if (io /= 0) then
      write (*, *) "post_convvelo: Failed to open temporary restart-file list."
      stop 1
    end if

    count = 0
    do
      read (unit, '(A)', iostat=io) line
      if (io /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do
    rewind (unit)

    allocate (files(count))
    files = ""
    do io = 1, count
      read (unit, '(A)') files(io)
    end do
    close (unit, status="delete")
    nfiles = count
  end subroutine discover_restart_files

  integer(C_INT) function getpid()
    interface
      integer(C_INT) function c_getpid() bind(C, name="getpid")
        use, intrinsic :: iso_c_binding
      end function c_getpid
    end interface
    getpid = c_getpid()
  end function getpid

end program post_convvelo
