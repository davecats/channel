#include "build_options.h"

! Writing the convection-velocity statistics to disk.
!
! The raw statistics file is a header -- the start and end time of the
! averaging window, then the number of samples behind the mean profiles --
! followed by the 3 + nPhi wall-normal mean profiles, then one full spectral
! field per accumulated statistic.  Profiles and fields are written through
! MPI-IO subarray views, so the file looks the same whatever the rank layout
! is; a companion `<name>.fields` text file names the fields in the order
! they appear, which is how the post-processing tools find them.
!
! What is on disk is the average, but what the accumulator holds is the
! running sum, so every field is divided by its own sample count on the way
! out.  That is average_convvelo_field, which convvelo also uses to hand a
! single averaged field back to a caller.
!
! Everything this module writes is passed in rather than read from convvelo,
! so that it sits below it and the dependency runs one way.  The caller is
! responsible for having the statistics on the host.
module convvelo_io

  use, intrinsic :: iso_c_binding
  use channel_grid
  use roctx, only: roctxPush, roctxPop
  use mpi_f08

  implicit none
  private

  ! Two doubles for the averaging window, then an 8-byte sample count.
  integer(C_INT64_T), parameter, public :: convvelo_file_header_bytes = 2_C_INT64_T*8_C_INT64_T + 8_C_INT64_T

  public :: write_convvelo_raw_file, write_convvelo_layout_file, average_convvelo_field

contains

  ! Writes the header, the mean profiles and every accumulated field.
  subroutine write_convvelo_raw_file(filename, n_fields, stats, n_field_samples, means, &
                                     n_mean_samples, average_start_time, average_end_time)
    implicit none

    character(len=*), intent(in) :: filename
    integer(C_INT), intent(in) :: n_fields
    integer(C_INT64_T), intent(in) :: n_field_samples(n_fields)
    integer(C_INT64_T), intent(in) :: n_mean_samples
    real(C_DOUBLE), intent(in) :: average_start_time, average_end_time
    complex(C_DOUBLE_COMPLEX), intent(in) :: stats(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, n_fields)
    complex(C_DOUBLE_COMPLEX), intent(in) :: means(ny0 - 2:nyN + 2, 1:3 + nPhi)

    integer(C_INT) :: field_index, i_profile, n_profiles
    real(C_DOUBLE) :: header_times(2)
    integer(C_INT64_T) :: header_sample_count
    complex(C_DOUBLE_COMPLEX), allocatable :: field_out(:, :, :)

    type(MPI_File) :: fh
    type(MPI_Status) :: status
    type(MPI_Datatype) :: file_type, mem_type, profile_file_type, profile_mem_type
    integer :: ierror
    integer, parameter :: ndims = 3
    integer, parameter :: ndims_profile = 1
    integer :: sizes(ndims), subsizes(ndims), starts(ndims)
    integer :: profile_sizes(ndims_profile), profile_subsizes(ndims_profile), profile_starts(ndims_profile)
    integer(C_INT) :: write_y0, write_yN, write_y_count, local_y_count
    integer(MPI_OFFSET_KIND) :: disp, field_bytes, profile_bytes, total_bytes

    n_profiles = 3 + nPhi
    header_times = [average_start_time, average_end_time]
    header_sample_count = n_mean_samples
    allocate (field_out(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))

    ! The two ghost rows outside the channel belong to the ranks that own the
    ! walls, so only those widen their share of the file.
    write_y0 = ny0
    write_yN = nyN
    if (ny0 == 1) write_y0 = -1_C_INT
    if (nyN == ny - 1) write_yN = ny + 1_C_INT
    write_y_count = write_yN - write_y0 + 1_C_INT
    local_y_count = nyN - ny0 + 5_C_INT

    sizes = [ny + 3, 2*nz + 1, nx + 1]
    subsizes = [write_y_count, 2*nz + 1, nxN - nx0 + 1]
    starts = [write_y0 + 1, 0, nx0]
    call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, file_type, ierror)
    call MPI_Type_commit(file_type, ierror)

    sizes = [local_y_count, 2*nz + 1, nxN - nx0 + 1]
    subsizes = [write_y_count, 2*nz + 1, nxN - nx0 + 1]
    starts = [write_y0 - (ny0 - 2), 0, 0]
    call MPI_Type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, mem_type, ierror)
    call MPI_Type_commit(mem_type, ierror)

    profile_sizes = [ny + 3]
    profile_subsizes = [write_y_count]
    profile_starts = [write_y0 + 1]
    call MPI_Type_create_subarray(ndims_profile, profile_sizes, profile_subsizes, profile_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, profile_file_type, ierror)
    call MPI_Type_commit(profile_file_type, ierror)

    profile_sizes = [local_y_count]
    profile_subsizes = [write_y_count]
    profile_starts = [write_y0 - (ny0 - 2)]
    call MPI_Type_create_subarray(ndims_profile, profile_sizes, profile_subsizes, profile_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, profile_mem_type, ierror)
    call MPI_Type_commit(profile_mem_type, ierror)

    profile_bytes = int(16, MPI_OFFSET_KIND)*int(ny + 3, MPI_OFFSET_KIND)
    field_bytes = int(16, MPI_OFFSET_KIND)*int(ny + 3, MPI_OFFSET_KIND)* &
                  int(2*nz + 1, MPI_OFFSET_KIND)*int(nx + 1, MPI_OFFSET_KIND)
    total_bytes = convvelo_file_header_bytes + int(n_profiles, MPI_OFFSET_KIND)*profile_bytes + &
                  int(n_fields, MPI_OFFSET_KIND)*field_bytes

    call MPI_File_open(MPI_COMM_WORLD, trim(filename), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
    call MPI_File_set_size(fh, total_bytes)

    if (iproc == 0) then
      call roctxPush("MPI_File_write_at convvelo_header")
      call MPI_File_write_at(fh, 0_MPI_OFFSET_KIND, header_times, 2, MPI_DOUBLE_PRECISION, status)
      call MPI_File_write_at(fh, 16_MPI_OFFSET_KIND, header_sample_count, 1, MPI_INTEGER8, status)
      call roctxPop("MPI_File_write_at convvelo_header")
    end if

    call roctxPush("MPI_File_write_all convvelo_profiles")
    do i_profile = 1, n_profiles
      disp = convvelo_file_header_bytes + int(i_profile - 1, MPI_OFFSET_KIND)*profile_bytes
      call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, profile_file_type, 'native', MPI_INFO_NULL)
      call write_profile_collective(fh, means(:, i_profile), profile_mem_type, status)
    end do
    call roctxPop("MPI_File_write_all convvelo_profiles")

    call roctxPush("MPI_File_write_all convvelo_fields")
    do field_index = 1, n_fields
      disp = convvelo_file_header_bytes + int(n_profiles, MPI_OFFSET_KIND)*profile_bytes + &
             int(field_index - 1, MPI_OFFSET_KIND)*field_bytes
      call MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, file_type, 'native', MPI_INFO_NULL)
      call average_convvelo_field(n_fields, stats, n_field_samples, field_index, field_out)
      call MPI_File_write_all(fh, field_out, 1, mem_type, status)
    end do
    call roctxPop("MPI_File_write_all convvelo_fields")

    call MPI_File_close(fh)
    call MPI_Type_free(file_type, ierror)
    call MPI_Type_free(mem_type, ierror)
    call MPI_Type_free(profile_file_type, ierror)
    call MPI_Type_free(profile_mem_type, ierror)
    deallocate (field_out)
  end subroutine write_convvelo_raw_file

  ! Every rank takes part -- MPI_File_write_all is collective -- but only the
  ! one holding the mean profile contributes anything.
  subroutine write_profile_collective(fh, profile, profile_mem_type, status)
    implicit none
    type(MPI_File), intent(in) :: fh
    complex(C_DOUBLE_COMPLEX), intent(in) :: profile(ny0 - 2:nyN + 2)
    type(MPI_Datatype), intent(in) :: profile_mem_type
    type(MPI_Status), intent(out) :: status
    integer :: write_count

    write_count = 0
    if (has_average) write_count = 1
    call MPI_File_write_all(fh, profile, write_count, profile_mem_type, status)
  end subroutine write_profile_collective

  ! Divides one accumulated field by its own sample count.  A field nothing
  ! was ever accumulated into comes back zero rather than as a division by
  ! zero, which is what a run that never reached its start time produces.
  subroutine average_convvelo_field(n_fields, stats, n_field_samples, field_index, field)
    implicit none
    integer(C_INT), intent(in) :: n_fields, field_index
    integer(C_INT64_T), intent(in) :: n_field_samples(n_fields)
    complex(C_DOUBLE_COMPLEX), intent(in) :: stats(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, n_fields)
    complex(C_DOUBLE_COMPLEX), intent(out) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    real(C_DOUBLE) :: inv_samples
    integer(C_INT) :: ix, iy, iz

    if (n_field_samples(field_index) <= 0_C_INT64_T) then
      field = (0.0d0, 0.0d0)
      return
    end if

    inv_samples = 1.0d0/dble(n_field_samples(field_index))
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0 - 2, nyN + 2
          field(iy, iz, ix) = inv_samples*stats(iy, iz, ix, field_index)
        end do
      end do
    end do
  end subroutine average_convvelo_field

  ! Names the profiles and fields, in file order, next to the file itself.
  subroutine write_convvelo_layout_file(filename, profile_names, velocity_names, scalar_names)
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: profile_names(:)
    character(len=*), intent(in) :: velocity_names(:)
    character(len=*), intent(in) :: scalar_names(:)
    character(len=512) :: layout_filename
    integer :: io, i

    if (iproc /= 0) return

    layout_filename = trim(filename)//".fields"
    open (unit=98, file=trim(layout_filename), status='replace', action='write', iostat=io)
    if (io /= 0) then
      write (*, *) 'ERROR: could not open convvelo field layout file: ', trim(layout_filename)
      stop 1
    end if

    write (98, '(A)', advance='no') 'profile_fields:'
    do i = 1, size(profile_names)
      write (98, '(1X,A)', advance='no') trim(profile_names(i))
    end do
    write (98, *)

    write (98, '(A)', advance='no') 'velocity_fields:'
    do i = 1, size(velocity_names)
      write (98, '(1X,A)', advance='no') trim(velocity_names(i))
    end do
    write (98, *)

    write (98, '(A)', advance='no') 'scalar_fields:'
    do i = 1, size(scalar_names)
      write (98, '(1X,A)', advance='no') trim(scalar_names(i))
    end do
    write (98, *)
    close (98)
  end subroutine write_convvelo_layout_file

end module convvelo_io
