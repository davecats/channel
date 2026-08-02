#include "header.h"

PROGRAM post_pressure
  USE, intrinsic :: iso_c_binding
  USE dnsdata, ONLY: iproc, read_restart_file, free_memory, V, ny, nz, nx, ny0, nyN, nx0, nxN
  USE pressure_output, ONLY: free_pressure_output, compute_pressure_output
  USE driver, ONLY: initialize
#ifdef HAVE_FFTW
  USE ffts, ONLY: free_fft, VVdz, VVdx, rVVdx
#endif
#ifdef HAVE_MPI
  USE mpi_f08
#endif
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: config_file = 'dns.in'
  CHARACTER(len=256), ALLOCATABLE :: restart_files(:)
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: p(:, :, :), dpdy(:, :, :)
  CHARACTER(len=128) :: step_id
  CHARACTER(len=256) :: p_filename, dpdy_filename
  INTEGER :: nfiles, ifile

  call discover_restart_files(restart_files, nfiles)
  if (nfiles == 0) then
    write (*, *) 'post_pressure: No restart files matching Dati.cart.*.out found in the current directory.'
    stop 1
  end if

  call initialize(config_file, trim(restart_files(1)), .FALSE.)

  allocate (p(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
  allocate (dpdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))

  do ifile = 1, nfiles
    if (ifile > 1) then
      if (iproc == 0) write (*, '(A,1X,A)') 'Reading', trim(restart_files(ifile))
      call read_restart_file(trim(restart_files(ifile)), V)
      !$omp target update to(V)
    end if

    if (iproc == 0) write (*, '(A,1X,A)') 'Processing', trim(restart_files(ifile))

    call compute_pressure_output(p_out=p, dpdy_out=dpdy)
    !$omp target update from(p, dpdy)

    step_id = extract_timestep(trim(restart_files(ifile)))
    p_filename = 'pField'//trim(step_id)//'.dat'
    dpdy_filename = 'dpdyField'//trim(step_id)//'.dat'

    if (iproc == 0) write (*, '(A,1X,A)') 'Processing finished', trim(restart_files(ifile))
    call write_field_mpi(trim(p_filename), p)
    call write_field_mpi(trim(dpdy_filename), dpdy)
  end do

  deallocate (p, dpdy, restart_files)

  call free_pressure_output()
#ifdef HAVE_FFTW
  call free_fft()
#endif
  call free_memory(.FALSE.)
#ifdef HAVE_MPI
  call MPI_Finalize()
#endif

CONTAINS

  SUBROUTINE discover_restart_files(files, nfiles)
    CHARACTER(len=256), ALLOCATABLE, INTENT(OUT) :: files(:)
    INTEGER, INTENT(OUT) :: nfiles
    CHARACTER(len=512) :: command, listfile, line
    INTEGER :: unit, io, count, exitstat, cmdstat
    INTEGER(C_INT) :: pid

    pid = process_id()
    write (listfile, '(A,I0)') '.post_pressure_files.', pid
    command = "find . -maxdepth 1 -type f -name 'Dati.cart.*.out' -printf '%f\n' | sort -V > "//trim(listfile)
    call execute_command_line(trim(command), wait=.true., exitstat=exitstat, cmdstat=cmdstat)
    if (cmdstat /= 0 .or. exitstat /= 0) then
      write (*, *) 'post_pressure: Failed to scan restart files in the current directory.'
      stop 1
    end if

    open (newunit=unit, file=trim(listfile), status='old', action='read', iostat=io)
    if (io /= 0) then
      write (*, *) 'post_pressure: Failed to open temporary restart-file list.'
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
    files = ''
    do io = 1, count
      read (unit, '(A)') files(io)
    end do
    close (unit, status='delete')
    nfiles = count
  END SUBROUTINE discover_restart_files

  FUNCTION extract_timestep(filename) RESULT(step)
    CHARACTER(len=*), INTENT(IN) :: filename
    CHARACTER(len=128) :: step
    INTEGER :: start_idx, end_idx

    step = trim(filename)
    start_idx = len('Dati.cart.') + 1
    end_idx = len_trim(filename) - len('.out')
    if (len_trim(filename) >= len('Dati.cart.') + len('.out') .and. filename(1:start_idx - 1) == 'Dati.cart.' .and. end_idx >= start_idx) then
      step = filename(start_idx:end_idx)
    end if
  END FUNCTION extract_timestep

  SUBROUTINE write_field_mpi(filename, field)
    CHARACTER(len=*), INTENT(IN) :: filename
    COMPLEX(C_DOUBLE_COMPLEX), INTENT(IN) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
#ifdef HAVE_MPI
    TYPE(MPI_File) :: fh
    TYPE(MPI_Status) :: status
    TYPE(MPI_Datatype) :: file_type, mem_type
    INTEGER :: ierror
    INTEGER, PARAMETER :: ndims = 3
    INTEGER :: sizes(ndims), subsizes(ndims), starts(ndims)
    INTEGER(MPI_OFFSET_KIND) :: disp

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

    disp = 0_MPI_OFFSET_KIND
    call MPI_File_open(MPI_COMM_WORLD, trim(filename), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
    call MPI_File_set_size(fh, 0_MPI_OFFSET_KIND)
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, file_type, 'native', MPI_INFO_NULL)
    call MPI_File_write_all(fh, field, 1, mem_type, status)
    call MPI_File_close(fh)
    call MPI_Type_free(file_type, ierror)
    call MPI_Type_free(mem_type, ierror)
#else
    INTEGER :: io
    open (newunit=io, file=filename, form='unformatted', access='stream', status='replace', action='write')
    write (io) field
    close (io)
#endif
  END SUBROUTINE write_field_mpi

  INTEGER(C_INT) FUNCTION process_id()
    INTERFACE
      INTEGER(C_INT) FUNCTION c_getpid() BIND(C, NAME='getpid')
        USE, intrinsic :: iso_c_binding
      END FUNCTION c_getpid
    END INTERFACE
    process_id = c_getpid()
  END FUNCTION process_id

END PROGRAM post_pressure
