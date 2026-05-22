program test_pressure_dpdy
  use, intrinsic :: iso_c_binding
  use dnsdata
  use driver
  implicit none

  character(len=256) :: restart_in, config_file
  integer(C_INT) :: iy, iz, ix
  complex(C_DOUBLE_COMPLEX), allocatable :: p(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: dpdy(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: p_ref(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: dpdy_ref(:, :, :)
  real(C_DOUBLE) :: local_max_err_p, global_max_err_p
  real(C_DOUBLE) :: local_max_err_dpdy, global_max_err_dpdy
  real(C_DOUBLE) :: diff, tol

  config_file = "tests/data/dns_test.in"
  restart_in = "tests/poisson/Dati.cart.35.out"

  call initialize(config_file, restart_in)

  allocate (p(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
  allocate (dpdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
  allocate (p_ref(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
  allocate (dpdy_ref(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))

!  call compute_poisson(p)
!  call compute_dpdy(dpdy)

  call read_reference_field_mpi('tests/poisson/pField35.fld', p_ref)
  call read_reference_field_mpi('tests/poisson/dpdyField35.fld', dpdy_ref)

  local_max_err_p = 0.0d0
  local_max_err_dpdy = 0.0d0

  do ix = nx0, nxN
    do iz = -nz, nz
      do iy = ny0 - 2, nyN + 2
        diff = abs(p(iy, iz, ix) - p_ref(iy, iz, ix))
        if (diff > local_max_err_p) local_max_err_p = diff

        diff = abs(dpdy(iy, iz, ix) - dpdy_ref(iy, iz, ix))
        if (diff > local_max_err_dpdy) local_max_err_dpdy = diff
      end do
    end do
  end do

  call MPI_Allreduce(local_max_err_p, global_max_err_p, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(local_max_err_dpdy, global_max_err_dpdy, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  if (iproc == 0) then
    write (*, *) 'Max error p   = ', global_max_err_p
    write (*, *) 'Max error dpdy= ', global_max_err_dpdy
  end if

  tol = 1.0e-12_C_DOUBLE
  if (global_max_err_dpdy < tol .and. global_max_err_p < tol) then
    write (*, *) 'Regression test PASSED'
  else
    write (*, *) 'Regression test FAILED'
#ifdef HAVE_MPI
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  end if

  call free_memory(.TRUE.)
#ifdef HAVE_MPI
  call MPI_Finalize()
#endif

contains

  subroutine read_reference_field_mpi(filename, ref)
    use, intrinsic :: iso_c_binding
    implicit none
    character(len=*), intent(in) :: filename
    complex(C_DOUBLE_COMPLEX), intent(out) :: ref(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
#ifdef HAVE_MPI
    type(MPI_File) :: fh
    type(MPI_Status) :: status
    type(MPI_Datatype) :: file_type, mem_type
    integer :: ierror
    integer, parameter :: ndims = 3
    integer :: sizes(ndims), subsizes(ndims), starts(ndims)
    integer(MPI_OFFSET_KIND) :: disp

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
    call MPI_File_open(MPI_COMM_WORLD, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh)
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, file_type, 'native', MPI_INFO_NULL)
    call MPI_File_read_all(fh, ref, 1, mem_type, status)
    call MPI_File_close(fh)
    call MPI_Type_free(file_type, ierror)
    call MPI_Type_free(mem_type, ierror)
#else
    integer :: io
    open (unit=99, file=filename, form='unformatted', access='stream', status='old', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open reference field: ', trim(filename)
      stop 1
    end if
    read (99) ref
    close (99)
#endif
  end subroutine read_reference_field_mpi

end program test_pressure_dpdy
