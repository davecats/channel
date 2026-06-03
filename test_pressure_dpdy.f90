program test_pressure_dpdy
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use dnsdata
  use pressure_output
  use driver
  implicit none

  character(len=256) :: restart_in, config_file
  character(len=256) :: arg
  integer(C_INT) :: iy, iz, ix
  complex(C_DOUBLE_COMPLEX), allocatable :: p(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: dpdy(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: p_ref(:, :, :)
  complex(C_DOUBLE_COMPLEX), allocatable :: dpdy_ref(:, :, :)
  real(C_DOUBLE) :: local_max_err_p, global_max_err_p
  real(C_DOUBLE) :: local_max_err_dpdy, global_max_err_dpdy
  real(C_DOUBLE) :: diff, tol
  logical :: local_has_nan_p, global_has_nan_p
  logical :: local_has_nan_dpdy, global_has_nan_dpdy
  integer(C_INT) :: max_ix_p, max_iz_p, max_iy_p
  integer(C_INT) :: max_ix_dpdy, max_iz_dpdy, max_iy_dpdy
  complex(C_DOUBLE_COMPLEX) :: max_p_val, max_p_ref
  complex(C_DOUBLE_COMPLEX) :: max_dpdy_val, max_dpdy_ref

  config_file = "tests/data/dns_test.in"
  restart_in = "tests/poisson/Dati.cart.35.out"
  call get_command_argument(1, arg)
  if (len_trim(arg) > 0) config_file = trim(arg)

  call initialize(config_file, restart_in, .false.)

  allocate (p(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
  allocate (dpdy(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
  allocate (p_ref(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
  allocate (dpdy_ref(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))

  call compute_pressure_output(p_out=p, dpdy_out=dpdy)

  call read_reference_field_mpi('tests/poisson/pField35.fld', p_ref)
  call read_reference_field_mpi('tests/poisson/dpdyField35.fld', dpdy_ref)

  local_max_err_p = 0.0d0
  local_max_err_dpdy = 0.0d0
  local_has_nan_p = .false.
  local_has_nan_dpdy = .false.
  max_ix_p = nx0; max_iz_p = -nz; max_iy_p = ny0 - 2
  max_ix_dpdy = nx0; max_iz_dpdy = -nz; max_iy_dpdy = ny0 - 2
  max_p_val = (0.0d0, 0.0d0)
  max_p_ref = (0.0d0, 0.0d0)
  max_dpdy_val = (0.0d0, 0.0d0)
  max_dpdy_ref = (0.0d0, 0.0d0)

  do ix = nx0, nxN
    do iz = -nz, nz
      do iy = ny0 - 2, nyN + 2
        if (ieee_is_nan(real(p(iy, iz, ix))) .or. ieee_is_nan(aimag(p(iy, iz, ix)))) then
          local_has_nan_p = .true.
        end if
        if (ieee_is_nan(real(dpdy(iy, iz, ix))) .or. ieee_is_nan(aimag(dpdy(iy, iz, ix)))) then
          local_has_nan_dpdy = .true.
        end if

        diff = abs(p(iy, iz, ix) - p_ref(iy, iz, ix))
        if (diff > local_max_err_p) then
          local_max_err_p = diff
          max_ix_p = ix
          max_iz_p = iz
          max_iy_p = iy
          max_p_val = p(iy, iz, ix)
          max_p_ref = p_ref(iy, iz, ix)
        end if

        diff = abs(dpdy(iy, iz, ix) - dpdy_ref(iy, iz, ix))
        if (diff > local_max_err_dpdy) then
          local_max_err_dpdy = diff
          max_ix_dpdy = ix
          max_iz_dpdy = iz
          max_iy_dpdy = iy
          max_dpdy_val = dpdy(iy, iz, ix)
          max_dpdy_ref = dpdy_ref(iy, iz, ix)
        end if
      end do
    end do
  end do

  call MPI_Allreduce(local_max_err_p, global_max_err_p, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(local_max_err_dpdy, global_max_err_dpdy, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(local_has_nan_p, global_has_nan_p, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
  call MPI_Allreduce(local_has_nan_dpdy, global_has_nan_dpdy, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)

  if (iproc == 0) then
    write (*, *) 'Max error p   = ', global_max_err_p
    write (*, *) 'Max error dpdy= ', global_max_err_dpdy
    write (*, *) 'NaN in p      = ', global_has_nan_p
    write (*, *) 'NaN in dpdy   = ', global_has_nan_dpdy
  end if

#ifdef HAVE_MPI
  if (abs(local_max_err_p - global_max_err_p) < 1.0d-15) then
    write (*, *) 'Worst p error rank/iy/iz/ix = ', iproc, max_iy_p, max_iz_p, max_ix_p
    write (*, *) 'Computed p / ref p          = ', max_p_val, max_p_ref
  end if
  if (abs(local_max_err_dpdy - global_max_err_dpdy) < 1.0d-15) then
    write (*, *) 'Worst dpdy error rank/iy/iz/ix = ', iproc, max_iy_dpdy, max_iz_dpdy, max_ix_dpdy
    write (*, *) 'Computed dpdy / ref dpdy       = ', max_dpdy_val, max_dpdy_ref
  end if
#else
  write (*, *) 'Worst p error iy/iz/ix = ', max_iy_p, max_iz_p, max_ix_p
  write (*, *) 'Computed p / ref p     = ', max_p_val, max_p_ref
  write (*, *) 'Worst dpdy error iy/iz/ix = ', max_iy_dpdy, max_iz_dpdy, max_ix_dpdy
  write (*, *) 'Computed dpdy / ref dpdy   = ', max_dpdy_val, max_dpdy_ref
#endif

  if (iproc == 0) then
    write (*, *) 'Sample p(:,0,0) computed/refo at iy=0,1,8,15,16 = '
    if (lbound(p, 1) <= 0 .and. ubound(p, 1) >= 0) write (*, *) p(0, 0, 0), p_ref(0, 0, 0)
    if (lbound(p, 1) <= 1 .and. ubound(p, 1) >= 1) write (*, *) p(1, 0, 0), p_ref(1, 0, 0)
    if (lbound(p, 1) <= 8 .and. ubound(p, 1) >= 8) write (*, *) p(8, 0, 0), p_ref(8, 0, 0)
    if (lbound(p, 1) <= 15 .and. ubound(p, 1) >= 15) write (*, *) p(15, 0, 0), p_ref(15, 0, 0)
    if (lbound(p, 1) <= 16 .and. ubound(p, 1) >= 16) write (*, *) p(16, 0, 0), p_ref(16, 0, 0)
    write (*, *) 'Sample dpdy(:,0,0) computed/refo at iy=0,1,8,15,16 = '
    if (lbound(dpdy, 1) <= 0 .and. ubound(dpdy, 1) >= 0) write (*, *) dpdy(0, 0, 0), dpdy_ref(0, 0, 0)
    if (lbound(dpdy, 1) <= 1 .and. ubound(dpdy, 1) >= 1) write (*, *) dpdy(1, 0, 0), dpdy_ref(1, 0, 0)
    if (lbound(dpdy, 1) <= 8 .and. ubound(dpdy, 1) >= 8) write (*, *) dpdy(8, 0, 0), dpdy_ref(8, 0, 0)
    if (lbound(dpdy, 1) <= 15 .and. ubound(dpdy, 1) >= 15) write (*, *) dpdy(15, 0, 0), dpdy_ref(15, 0, 0)
    if (lbound(dpdy, 1) <= 16 .and. ubound(dpdy, 1) >= 16) write (*, *) dpdy(16, 0, 0), dpdy_ref(16, 0, 0)
  end if

  tol = 1.0e-12_C_DOUBLE
  if ((.not. global_has_nan_p) .and. (.not. global_has_nan_dpdy) .and. global_max_err_p < tol .and. global_max_err_dpdy < tol) then
    write (*, *) 'Regression test PASSED'
  else
    write (*, *) 'Regression test FAILED'
#ifdef HAVE_MPI
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  end if

  call free_pressure_output()
  call free_memory(.FALSE.)
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
    subsizes = [nyN - ny0 + 5, 2*nz + 1, nxN - nx0 + 1]
    starts = [ny0 - 1, 0, nx0]
    call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, file_type, ierror)
    call MPI_Type_commit(file_type, ierror)

    sizes = [nyN - ny0 + 5, 2*nz + 1, nxN - nx0 + 1]
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
