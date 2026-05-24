program test_convvelo_stats
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic
  use dnsdata
  use pressure_output
  use driver
  implicit none

  character(len=*), parameter :: config_file = "tests/convvelo/dns.in"
  character(len=*), parameter :: restart_in = "tests/convvelo/Dati.cart.out"
  character(len=*), parameter :: raw_stats_file = "tests/convvelo/raw_statistics.bin"
  integer(C_INT), parameter :: n_velocity_fields = 33
  integer(C_INT), parameter :: n_scalar_fields = 10

  complex(C_DOUBLE_COMPLEX), allocatable :: convvelo(:, :, :, :)
  integer(C_INT) :: i_field, i_phi, nfail, global_nfail
  real(C_DOUBLE) :: tol

  ! This is a placeholder regression harness for the online convvelo
  ! statistics. The future accumulation entry points are sketched below and
  ! the computed storage is kept local to this test until the implementation
  ! exists in the solver.

  tol = 1.0d-12
  nfail = 0

  call initialize(config_file, restart_in)
  allocate (convvelo(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, n_velocity_fields + nPhi*n_scalar_fields))
  convvelo = (0.0d0, 0.0d0)

  ! Future entry points:
  ! call reset_convvelo_stats()
  ! call read_restart_file("tests/convvelo/Dati.cart.39.out", V)
  ! call acc_convvelo_stats()
  ! call read_restart_file("tests/convvelo/Dati.cart.40.out", V)
  ! call acc_convvelo_stats()
  ! call read_restart_file("tests/convvelo/Dati.cart.41.out", V)
  ! call acc_convvelo_stats()

  call read_restart_file("tests/convvelo/Dati.cart.39.out", V)
  call read_restart_file("tests/convvelo/Dati.cart.40.out", V)
  call read_restart_file("tests/convvelo/Dati.cart.41.out", V)

  i_field = 0

  call compare_field("u_cross_u", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("u_cross_dyu", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("u_cross_v", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("u_cross_dyv", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("u_cross_w", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("u_cross_dyw", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("u_cross_dyyu", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1

  call compare_field("v_cross_u", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_dyu", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_v", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_dyv", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_w", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_dyw", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_dyyv", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1

  call compare_field("w_cross_u", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_dyu", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_v", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_dyv", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_w", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_dyw", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_dyyw", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1

  call compare_field("u_cross_p", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_dpdy", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_p", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1

  call compare_field("u_cross_uu", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("u_cross_uw", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_uv", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_vw", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_uw", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_ww", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1

  call compare_field("u_cross_dyuv", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("v_cross_dyvv", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1
  call compare_field("w_cross_dyvw", convvelo(:, :, :, i_field + 1), i_field, tol, nfail); i_field = i_field + 1

  do i_phi = 1, nPhi
    call compare_scalar_field("t_theta_theta", convvelo, i_phi, 0, tol, nfail)
    call compare_scalar_field("t_theta_u", convvelo, i_phi, 1, tol, nfail)
    call compare_scalar_field("t_theta_v", convvelo, i_phi, 2, tol, nfail)
    call compare_scalar_field("t_theta_w", convvelo, i_phi, 3, tol, nfail)
    call compare_scalar_field("t_theta_thetau", convvelo, i_phi, 4, tol, nfail)
    call compare_scalar_field("t_theta_thetaw", convvelo, i_phi, 5, tol, nfail)
    call compare_scalar_field("t_theta_dyytheta", convvelo, i_phi, 6, tol, nfail)
    call compare_scalar_field("t_theta_dythetav", convvelo, i_phi, 7, tol, nfail)
    call compare_scalar_field("t_theta_dytheta", convvelo, i_phi, 8, tol, nfail)
    call compare_scalar_field("t_theta_dyv", convvelo, i_phi, 9, tol, nfail)
  end do

  global_nfail = nfail
#ifdef HAVE_MPI
  call MPI_Allreduce(nfail, global_nfail, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

  if (iproc == 0) then
    if (global_nfail == 0) then
      write (*, *) "Convvelo statistics regression PASSED"
    else
      write (*, *) "Convvelo statistics regression FAILED, failing fields = ", global_nfail
    end if
  end if

  if (global_nfail /= 0) then
#ifdef HAVE_MPI
    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
#endif
  end if

  call free_pressure_output()
  deallocate (convvelo)
  call free_memory(.true.)
#ifdef HAVE_MPI
  call MPI_Finalize()
#endif

contains

  subroutine compare_scalar_field(base_name, computed_all, i_phi_local, scalar_field_index, tol_local, nfail_local)
    implicit none
    character(len=*), intent(in) :: base_name
   complex(C_DOUBLE_COMPLEX), intent(in) :: computed_all(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, n_velocity_fields + nPhi*n_scalar_fields)
    integer(C_INT), intent(in) :: i_phi_local, scalar_field_index
    real(C_DOUBLE), intent(in) :: tol_local
    integer(C_INT), intent(inout) :: nfail_local
    character(len=128) :: full_name
    integer(C_INT) :: raw_field_index

    write (full_name, '(A,"[phi=",I0,"]")') trim(base_name), i_phi_local
    raw_field_index = n_velocity_fields + (i_phi_local - 1)*n_scalar_fields + scalar_field_index
    call compare_field(trim(full_name), computed_all(:, :, :, raw_field_index + 1), raw_field_index, tol_local, nfail_local)
  end subroutine compare_scalar_field

  subroutine compare_field(name, computed, raw_field_index, tol_local, nfail_local)
    use, intrinsic :: ieee_arithmetic
    implicit none
    character(len=*), intent(in) :: name
    complex(C_DOUBLE_COMPLEX), intent(in) :: computed(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    integer(C_INT), intent(in) :: raw_field_index
    real(C_DOUBLE), intent(in) :: tol_local
    integer(C_INT), intent(inout) :: nfail_local

    complex(C_DOUBLE_COMPLEX), allocatable :: reference(:, :, :)
    real(C_DOUBLE) :: local_max_err, global_max_err, diff
    logical :: local_has_nan, global_has_nan
    integer(C_INT) :: iy, iz, ix
    integer(C_INT) :: worst_iy, worst_iz, worst_ix
    complex(C_DOUBLE_COMPLEX) :: worst_computed, worst_reference

    allocate (reference(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))
    call read_raw_stat_field_mpi(raw_stats_file, raw_field_index, reference)

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
          if (ieee_is_nan(real(computed(iy, iz, ix))) .or. ieee_is_nan(aimag(computed(iy, iz, ix)))) then
            local_has_nan = .true.
          end if

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

    if (iproc == 0) then
      write (*, '(A,": max_err=",ES12.4,", has_nan=",L1)') trim(name), global_max_err, global_has_nan
    end if

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

    deallocate (reference)
  end subroutine compare_field

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
    integer(MPI_OFFSET_KIND) :: field_bytes

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

    field_bytes = int(16, MPI_OFFSET_KIND)*int(ny + 3, MPI_OFFSET_KIND)* &
                  int(2*nz + 1, MPI_OFFSET_KIND)*int(nx + 1, MPI_OFFSET_KIND)
    disp = int(field_index, MPI_OFFSET_KIND)*field_bytes

    call MPI_File_open(MPI_COMM_WORLD, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh)
    call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, file_type, 'native', MPI_INFO_NULL)
    call MPI_File_read_all(fh, ref, 1, mem_type, status)
    call MPI_File_close(fh)
    call MPI_Type_free(file_type, ierror)
    call MPI_Type_free(mem_type, ierror)
#else
    integer :: io
    integer(C_INT64_T) :: field_bytes, pos_bytes

    field_bytes = int(16, C_INT64_T)*int(ny + 3, C_INT64_T)*int(2*nz + 1, C_INT64_T)*int(nx + 1, C_INT64_T)
    pos_bytes = int(field_index, C_INT64_T)*field_bytes + 1_C_INT64_T

    open (unit=99, file=filename, form='unformatted', access='stream', status='old', action='read', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open convvelo reference file: ', trim(filename)
      stop 1
    end if
    read (99, pos=pos_bytes) ref
    close (99)
#endif
  end subroutine read_raw_stat_field_mpi

end program test_convvelo_stats
