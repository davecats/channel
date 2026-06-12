!============================================!
!                                            !
!            Distributed Traspose            !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
!
! Author: Dr. Davide Gatti
! Date  : 15/Apr/2019
!

#include "header.h"

#define SPLIT_BLOCK(total, nparts, part, start, count) \
  count = (total)/(nparts); \
  if ((part) < mod((total), (nparts))) count = count + 1; \
  start = (part)*((total)/(nparts)) + min((part), mod((total), (nparts))) + 1

#define YSLAB_RANGE_FROM_CACHE(rank, ny, lower_pad, upper_pad, include_physical_ghosts, first_y, last_y) \
  first_y = yslab_active_first_y_by_rank(rank) - (lower_pad); \
  if ((include_physical_ghosts) .and. (rank) == 0) first_y = min(first_y, -1_C_INT); \
  if ((rank) == npy_grid - 1) then; \
    last_y = (ny) - 1 + (upper_pad); \
  else; \
    last_y = yslab_active_first_y_by_rank((rank) + 1) - 1 + (upper_pad); \
  end if; \
  if ((include_physical_ghosts) .and. (rank) == npy_grid - 1) last_y = max(last_y, (ny) + 1)

  MODULE mpi_transpose

    USE, intrinsic :: iso_c_binding
    USE, intrinsic :: iso_fortran_env
#ifdef HAVE_MPI
    USE mpi_f08
#endif
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    use omp_lib
#endif
    use roctx

    IMPLICIT NONE

#ifdef HAVE_MPI
    TYPE(MPI_Comm) :: MPI_CART_COMM, MPI_COMM_X, MPI_COMM_Y
#endif
#if defined(HAVE_HIP)
    complex(C_DOUBLE_COMPLEX), pointer:: sendbuf(:, :), recvbuf(:, :)
#else
    complex(C_DOUBLE_COMPLEX), allocatable :: sendbuf(:, :), recvbuf(:, :)
#endif
    complex(C_DOUBLE_COMPLEX), allocatable, target, save :: yslab_workspace(:, :)
    integer(C_INT), save :: yslab_scratch_rows = -1
    integer(C_INT), save :: yslab_scratch_lines = -1
    integer(C_INT), allocatable, save :: yslab_first_line_by_rank(:), yslab_line_count_by_rank(:)
    integer(C_INT), allocatable, save :: yslab_active_first_y_by_rank(:)
    integer(C_INT), save :: yslab_owned_first_line = 1
    integer(C_INT), save :: yslab_owned_line_count = 0
    integer(C_INT), save :: nproc, iproc, ierr, nzd, nx
    integer(C_INT), save :: npy_grid = 1, npxz = 1, ipy = 0, ipxz = 0
    integer(C_INT), save :: nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN, sendcount
    !$omp declare target(npy_grid, npxz, ipy, ipxz, nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN, sendcount)

    logical, save :: has_terminal, has_average, fft_transpose_is_local
#ifdef HAVE_MPI
    TYPE(MPI_Datatype), save :: writeview_type, owned2write_type, vel_read_type, vel_field_type
#endif

  CONTAINS

    subroutine setup_y_decomposition(ny, nlines_z, nx_lines)
      implicit none
      integer(C_INT), intent(in) :: ny, nlines_z, nx_lines
      integer(C_INT) :: rank, first_x, x_count, base_count, remainder

      ! For y slabs, the x subrange is further decomposed instead of a pencil decomposition
      ! That way we can always keep the full z range in the local workspace, which makes the code simpler
      ! Computationally it makes no difference since the line solve does not require any x/z ghost information,
      ! and the number of systems we solve per rank is identical

      if (allocated(yslab_first_line_by_rank)) then
        !$omp target exit data map(delete: yslab_first_line_by_rank, yslab_line_count_by_rank, &
        !$omp& yslab_active_first_y_by_rank)
        deallocate (yslab_first_line_by_rank, yslab_line_count_by_rank, &
                    yslab_active_first_y_by_rank)
      end if

      allocate (yslab_first_line_by_rank(0:npy_grid - 1), yslab_line_count_by_rank(0:npy_grid - 1), &
                yslab_active_first_y_by_rank(0:npy_grid - 1))

      base_count = nx_lines/npy_grid
      remainder = mod(nx_lines, npy_grid)

      do rank = 0, npy_grid - 1
        x_count = base_count
        if (rank < remainder) x_count = x_count + 1
        first_x = rank*base_count + min(rank, remainder)

        yslab_line_count_by_rank(rank) = x_count*nlines_z
        yslab_first_line_by_rank(rank) = first_x*nlines_z + 1

        yslab_active_first_y_by_rank(rank) = 1 + rank*(ny - 1)/npy_grid
      end do

      yslab_owned_first_line = yslab_first_line_by_rank(ipy)
      yslab_owned_line_count = yslab_line_count_by_rank(ipy)

      !$omp target enter data map(to: yslab_first_line_by_rank, yslab_line_count_by_rank, &
      !$omp& yslab_active_first_y_by_rank)
    end subroutine setup_y_decomposition

    subroutine prepare_yslab_scratch(nrows, nlines)
      implicit none
      integer(C_INT), intent(in) :: nrows, nlines

      if (nlines <= 0) return
      if (allocated(yslab_workspace)) then
        ! Treat nlines as a minimum capacity.  This lets pressure_output reserve
        ! extra y-slab columns for transposed V while the common line solve reuses
        ! the first block for the RHS/solution.
        if (yslab_scratch_rows /= nrows .or. yslab_scratch_lines < nlines) then
          !$omp target exit data map(delete: yslab_workspace)
          deallocate (yslab_workspace)
          yslab_scratch_rows = -1
          yslab_scratch_lines = -1
        end if
      end if
      if (.not. allocated(yslab_workspace)) then
        allocate (yslab_workspace(-1:nrows - 2, nlines))
        !$omp target enter data map(alloc: yslab_workspace)
        yslab_scratch_rows = nrows
        yslab_scratch_lines = nlines
      end if
    end subroutine prepare_yslab_scratch

    subroutine release_yslab_scratch()
      implicit none

      if (.not. allocated(yslab_workspace)) return
      !$omp target exit data map(delete: yslab_workspace)
      deallocate (yslab_workspace)
      yslab_scratch_rows = -1
      yslab_scratch_lines = -1
    end subroutine release_yslab_scratch

#ifdef HAVE_MPI
    subroutine yslab_transpose_to_full(field, slab, ny, nz, nlines, nlines_z, include_physical_ghosts)
      implicit none
      integer(C_INT), intent(in) :: ny, nz, nlines, nlines_z
      complex(C_DOUBLE_COMPLEX), intent(in) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: slab(-1:, :)
      logical, intent(in) :: include_physical_ghosts
      complex(C_DOUBLE_COMPLEX), allocatable :: send(:), recv(:)
      integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
      integer(C_INT) :: dest, src, first_line, line_count
      integer(C_INT) :: y_first, y_last, rows, local_rows, peer_rows, total_send, total_recv
      integer(C_INT) :: ilocal, iline, ix, iz, iy, p

      if (npy_grid == 1) then
        call roctxPush("yslab_to_full local_copy")
        call yslab_copy_to_full(field, slab, ny, nz, 1_C_INT, nlines, nlines_z)
        call roctxPop("yslab_to_full local_copy")
        return
      end if

      call roctxPush("yslab_to_full setup_counts")
      allocate (send_counts(npy_grid), recv_counts(npy_grid), send_displs(npy_grid), recv_displs(npy_grid))
      total_send = 0
      total_recv = 0
      YSLAB_RANGE_FROM_CACHE(ipy, ny, 0_C_INT, 0_C_INT, include_physical_ghosts, y_first, y_last)
      local_rows = y_last - y_first + 1
      do dest = 0, npy_grid - 1
        line_count = yslab_line_count_by_rank(dest)
        send_counts(dest + 1) = line_count*local_rows
        send_displs(dest + 1) = total_send
        total_send = total_send + send_counts(dest + 1)

        YSLAB_RANGE_FROM_CACHE(dest, ny, 0_C_INT, 0_C_INT, include_physical_ghosts, y_first, y_last)
        peer_rows = y_last - y_first + 1
        recv_counts(dest + 1) = yslab_owned_line_count*peer_rows
        recv_displs(dest + 1) = total_recv
        total_recv = total_recv + recv_counts(dest + 1)
      end do
      call roctxPop("yslab_to_full setup_counts")

      call roctxPush("yslab_to_full allocate_buffers")
      allocate (send(max(1, total_send)), recv(max(1, total_recv)))
      !$omp target enter data map(alloc: send, recv) map(to: send_displs, recv_displs)
      call roctxPop("yslab_to_full allocate_buffers")

      YSLAB_RANGE_FROM_CACHE(ipy, ny, 0_C_INT, 0_C_INT, include_physical_ghosts, y_first, y_last)
      rows = y_last - y_first + 1
      call roctxPush("yslab_to_full pack")
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(field, send, send_displs, nlines_z, rows, y_first, y_last, nlines, nz, nx0, npy_grid) &
      !$omp shared(yslab_first_line_by_rank, yslab_line_count_by_rank) &
      !$omp private(dest, ilocal, iy, first_line, line_count, iline, ix, iz, p)
      do dest = 0, npy_grid - 1
        do ilocal = 1, nlines
          do iy = y_first, y_last
            first_line = yslab_first_line_by_rank(dest)
            line_count = yslab_line_count_by_rank(dest)
            if (ilocal > line_count) cycle
            iline = first_line + ilocal - 1
            ix = (iline - 1)/nlines_z + nx0
            iz = mod(iline - 1, nlines_z) - nz
            p = send_displs(dest + 1) + (ilocal - 1)*rows + (iy - y_first + 1)
            send(p) = field(iy, iz, ix)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
      call roctxPop("yslab_to_full pack")

      call roctxPush("MPI_Alltoallv yslab_to_full")
#ifndef HAVE_HIP
      !$omp target data use_device_addr(send, recv)
#endif
      call MPI_Alltoallv(send, send_counts, send_displs, MPI_DOUBLE_COMPLEX, &
                         recv, recv_counts, recv_displs, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr)
#ifndef HAVE_HIP
      !$omp end target data
#endif
      call roctxPop("MPI_Alltoallv yslab_to_full")

      call roctxPush("yslab_to_full unpack")
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(slab, recv, recv_displs, yslab_owned_line_count, include_physical_ghosts, ny, npy_grid) &
      !$omp shared(yslab_active_first_y_by_rank) &
      !$omp private(src, ilocal, iy, y_first, y_last, rows, p)
      do src = 0, npy_grid - 1
        do ilocal = 1, yslab_owned_line_count
          do iy = -1, ny + 1
            YSLAB_RANGE_FROM_CACHE(src, ny, 0_C_INT, 0_C_INT, include_physical_ghosts, y_first, y_last)
            if (iy < y_first .or. iy > y_last) cycle
            rows = y_last - y_first + 1
            p = recv_displs(src + 1) + (ilocal - 1)*rows + (iy - y_first + 1)
            slab(iy, ilocal) = recv(p)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
      call roctxPop("yslab_to_full unpack")

      call roctxPush("yslab_to_full cleanup")
      !$omp target exit data map(delete: send, recv, send_displs, recv_displs)
      deallocate (send, recv, send_counts, recv_counts, send_displs, recv_displs)
      call roctxPop("yslab_to_full cleanup")
    end subroutine yslab_transpose_to_full

    subroutine yslab_transpose_from_full(slab, field, ny, nz, nlines, nlines_z)
      implicit none
      integer(C_INT), intent(in) :: ny, nz, nlines, nlines_z
      complex(C_DOUBLE_COMPLEX), intent(in) :: slab(-1:, :)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
      complex(C_DOUBLE_COMPLEX), allocatable :: send(:), recv(:)
      integer, allocatable :: send_counts(:), recv_counts(:), send_displs(:), recv_displs(:)
      integer(C_INT) :: dest, src, first_line, line_count
      integer(C_INT) :: y_first, y_last, rows, local_rows, peer_rows, total_send, total_recv
      integer(C_INT) :: ilocal, iline, ix, iz, iy, p

      if (npy_grid == 1) then
        call roctxPush("yslab_from_full local_copy")
        call yslab_copy_from_full(slab, field, ny, nz, 1_C_INT, nlines, nlines_z)
        call roctxPop("yslab_from_full local_copy")
        return
      end if

      call roctxPush("yslab_from_full setup_counts")
      allocate (send_counts(npy_grid), recv_counts(npy_grid), send_displs(npy_grid), recv_displs(npy_grid))
      total_send = 0
      total_recv = 0
      YSLAB_RANGE_FROM_CACHE(ipy, ny, 2_C_INT, 2_C_INT, .true., y_first, y_last)
      local_rows = y_last - y_first + 1
      do dest = 0, npy_grid - 1
        YSLAB_RANGE_FROM_CACHE(dest, ny, 2_C_INT, 2_C_INT, .true., y_first, y_last)
        peer_rows = y_last - y_first + 1
        send_counts(dest + 1) = yslab_owned_line_count*peer_rows
        send_displs(dest + 1) = total_send
        total_send = total_send + send_counts(dest + 1)

        line_count = yslab_line_count_by_rank(dest)
        recv_counts(dest + 1) = line_count*local_rows
        recv_displs(dest + 1) = total_recv
        total_recv = total_recv + recv_counts(dest + 1)
      end do
      call roctxPop("yslab_from_full setup_counts")

      call roctxPush("yslab_from_full allocate_buffers")
      allocate (send(max(1, total_send)), recv(max(1, total_recv)))
      !$omp target enter data map(alloc: send, recv) map(to: send_displs, recv_displs)
      call roctxPop("yslab_from_full allocate_buffers")

      call roctxPush("yslab_from_full pack")
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(slab, send, send_displs, yslab_owned_line_count, ny, npy_grid) &
      !$omp shared(yslab_active_first_y_by_rank) &
      !$omp private(dest, ilocal, iy, y_first, y_last, rows, p)
      do dest = 0, npy_grid - 1
        do ilocal = 1, yslab_owned_line_count
          do iy = -1, ny + 1
            YSLAB_RANGE_FROM_CACHE(dest, ny, 2_C_INT, 2_C_INT, .true., y_first, y_last)
            if (iy < y_first .or. iy > y_last) cycle
            rows = y_last - y_first + 1
            p = send_displs(dest + 1) + (ilocal - 1)*rows + (iy - y_first + 1)
            send(p) = slab(iy, ilocal)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
      call roctxPop("yslab_from_full pack")

      call roctxPush("MPI_Alltoallv yslab_from_full")
#ifndef HAVE_HIP
      !$omp target data use_device_addr(send, recv)
#endif
      call MPI_Alltoallv(send, send_counts, send_displs, MPI_DOUBLE_COMPLEX, &
                         recv, recv_counts, recv_displs, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr)
#ifndef HAVE_HIP
      !$omp end target data
#endif
      call roctxPop("MPI_Alltoallv yslab_from_full")

      YSLAB_RANGE_FROM_CACHE(ipy, ny, 2_C_INT, 2_C_INT, .true., y_first, y_last)
      rows = y_last - y_first + 1
      call roctxPush("yslab_from_full unpack")
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(field, recv, recv_displs, nlines_z, rows, y_first, y_last, nz, nx0, nlines, npy_grid) &
      !$omp shared(yslab_first_line_by_rank, yslab_line_count_by_rank) &
      !$omp private(src, ilocal, iy, first_line, line_count, iline, ix, iz, p)
      do src = 0, npy_grid - 1
        do ilocal = 1, nlines
          do iy = y_first, y_last
            first_line = yslab_first_line_by_rank(src)
            line_count = yslab_line_count_by_rank(src)
            if (ilocal > line_count) cycle
            iline = first_line + ilocal - 1
            ix = (iline - 1)/nlines_z + nx0
            iz = mod(iline - 1, nlines_z) - nz
            p = recv_displs(src + 1) + (ilocal - 1)*rows + (iy - y_first + 1)
            field(iy, iz, ix) = recv(p)
          end do
        end do
      end do
      !$omp end target teams distribute parallel do
      call roctxPop("yslab_from_full unpack")

      call roctxPush("yslab_from_full cleanup")
      !$omp target exit data map(delete: send, recv, send_displs, recv_displs)
      deallocate (send, recv, send_counts, recv_counts, send_displs, recv_displs)
      call roctxPop("yslab_from_full cleanup")
    end subroutine yslab_transpose_from_full
#endif

    subroutine yslab_copy_to_full(field, slab, ny, nz, first_line, line_count, nlines_z)
      implicit none
      integer(C_INT), intent(in) :: ny, nz, first_line, line_count, nlines_z
      complex(C_DOUBLE_COMPLEX), intent(in) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: slab(-1:, :)
      integer(C_INT) :: ilocal, iline, ix, iz, iy

      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp shared(field, slab, first_line, line_count, nlines_z, nz, nx0, ny) &
      !$omp private(ilocal, iline, ix, iz, iy)
      do ilocal = 1, line_count
        do iy = -1, ny + 1
          iline = first_line + ilocal - 1
          ix = (iline - 1)/nlines_z + nx0
          iz = mod(iline - 1, nlines_z) - nz
          slab(iy, ilocal) = field(iy, iz, ix)
        end do
      end do
      !$omp end target teams distribute parallel do
    end subroutine yslab_copy_to_full

    subroutine yslab_copy_from_full(slab, field, ny, nz, first_line, line_count, nlines_z)
      implicit none
      integer(C_INT), intent(in) :: ny, nz, first_line, line_count, nlines_z
      complex(C_DOUBLE_COMPLEX), intent(in) :: slab(-1:, :)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
      integer(C_INT) :: ilocal, iline, ix, iz, iy

      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp shared(slab, field, first_line, line_count, nlines_z, nz, nx0, ny) &
      !$omp private(ilocal, iline, ix, iz, iy)
      do ilocal = 1, line_count
        do iy = -1, ny + 1
          iline = first_line + ilocal - 1
          ix = (iline - 1)/nlines_z + nx0
          iz = mod(iline - 1, nlines_z) - nz
          field(iy, iz, ix) = slab(iy, ilocal)
        end do
      end do
      !$omp end target teams distribute parallel do
    end subroutine yslab_copy_from_full

    SUBROUTINE repack_zTOx_local(Vz, Vx, ny)
      use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
      implicit none
      complex(C_DOUBLE_COMPLEX), intent(in) :: Vz(1:, 1:, :)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(1:, 1:, :)
      integer(C_INT), intent(in) :: ny
      integer(C_SIZE_T) :: iy, ix, iz
      integer(C_INT) :: ny_batch

      ny_batch = size(Vz, 3)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(Vz, Vx) shared(ny_batch, nxB, nzd) private(iy, ix, iz)
      do iy = 1, ny_batch
        do ix = 1, nxB
          do iz = 1, nzd
            Vx(ix, iz, iy) = Vz(iz, ix, iy)
          end do
        end do
      end do
    END SUBROUTINE repack_zTOx_local

    SUBROUTINE repack_xTOz_local(Vx, Vz, ny)
      use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
      implicit none
      complex(C_DOUBLE_COMPLEX), intent(in) :: Vx(1:, 1:, :)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(1:, 1:, :)
      integer(C_INT), intent(in) :: ny
      integer(C_SIZE_T) :: iy, ix, iz
      integer(C_INT) :: ny_batch

      ny_batch = size(Vx, 3)
      !$omp target teams distribute parallel do collapse(3) default(none) &
      !$omp shared(Vx, Vz) shared(ny_batch, nxB, nzd) private(iy, ix, iz)
      do iy = 1, ny_batch
        do iz = 1, nzd
          do ix = 1, nxB
            Vz(iz, ix, iy) = Vx(ix, iz, iy)
          end do
        end do
      end do
    END SUBROUTINE repack_xTOz_local

    SUBROUTINE pack_zTOx(Vz, send, ny)
      use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
      implicit none
      complex(C_DOUBLE_COMPLEX), intent(in)  :: Vz(1:, 1:, :)
      complex(C_DOUBLE_COMPLEX), intent(out) :: send(:)
      integer(C_INT), intent(in)  :: ny
      integer(C_SIZE_T) :: iy, ix, iz, dest, p
      integer(C_INT) :: ny_batch

      ny_batch = size(Vz, 3)
      !$omp target teams distribute parallel do collapse(4) default(none) &
      !$omp shared(Vz, send) shared(ny_batch, nxB, nzB, npxz, sendcount) private(iy, ix, iz, dest, p)
      do dest = 0, npxz - 1
        do iy = 1, ny_batch
          do ix = 1, nxB
            do iz = 1, nzB
              p = dest*sendcount + iz + (nzB*(ix - 1)) + (nzB*nxB*(iy - 1))
              send(p) = Vz(dest*nzB + iz, ix, iy)
            end do
          end do
        end do
      end do

    END SUBROUTINE pack_zTOx

    SUBROUTINE unpack_zTOx(recv, Vx, ny)
      use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
      implicit none
      complex(C_DOUBLE_COMPLEX), intent(in)  :: recv(:)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(1:, 1:, :)
      integer(C_INT), intent(in)  :: ny
      integer(C_SIZE_T) :: iy, ix, iz, src, p
      integer(C_INT) :: ny_batch

      ny_batch = size(Vx, 3)
      !$omp target teams distribute parallel do collapse(4) default(none) &
      !$omp shared(Vx, recv) shared(ny_batch, nxB, nzB, npxz, sendcount) private(iy, ix, iz, src, p)
      do src = 0, npxz - 1
        do iy = 1, ny_batch
          do ix = 1, nxB
            do iz = 1, nzB
              p = src*sendcount + iz + (nzB*(ix - 1)) + (nzB*nxB*(iy - 1))
              Vx(ix + src*nxB, iz, iy) = recv(p)
            end do
          end do
        end do
      end do
    END SUBROUTINE unpack_zTOx

    SUBROUTINE pack_xTOz(Vx, send, ny)
      use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
      implicit none
      complex(C_DOUBLE_COMPLEX), intent(in)  :: Vx(1:, 1:, :)
      complex(C_DOUBLE_COMPLEX), intent(out) :: send(:)
      integer(C_INT), intent(in)  :: ny
      integer(C_SIZE_T) :: iy, ix, iz, dest, p
      integer(C_INT) :: ny_batch

      ny_batch = size(Vx, 3)
      !$omp target teams distribute parallel do collapse(4) default(none) &
      !$omp shared(Vx, send) shared(ny_batch, nxB, nzB, npxz, sendcount) private(iy, ix, iz, dest, p)
      do dest = 0, npxz - 1
        do iy = 1, ny_batch
          do iz = 1, nzB
            do ix = 1, nxB
              p = dest*sendcount + ix + (nxB*(iz - 1)) + (nxB*nzB*(iy - 1))
              send(p) = Vx(dest*nxB + ix, iz, iy)
            end do
          end do
        end do
      end do
    END SUBROUTINE pack_xTOz

    SUBROUTINE unpack_xTOz(recv, Vz, ny)
      use iso_c_binding, only: C_INT, C_SIZE_T, C_DOUBLE_COMPLEX
      implicit none
      complex(C_DOUBLE_COMPLEX), intent(in)  :: recv(:)
      complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(1:, 1:, :)
      integer(C_INT), intent(in)  :: ny
      integer(C_SIZE_T) :: iy, ix, iz, src, p
      integer(C_INT) :: ny_batch

      ny_batch = size(Vz, 3)
      !$omp target teams distribute parallel do collapse(4) default(none) &
      !$omp shared(Vz, recv) shared(ny_batch, nxB, nzB, npxz, sendcount) private(iy, ix, iz, src, p)
      do src = 0, npxz - 1
        do iy = 1, ny_batch
          do iz = 1, nzB
            do ix = 1, nxB
              p = src*sendcount + ix + (nxB*(iz - 1)) + (nxB*nzB*(iy - 1))
              Vz(iz + src*nzB, ix, iy) = recv(p)
            end do
          end do
        end do
      end do
    END SUBROUTINE unpack_xTOz

    SUBROUTINE alltoall(send, recv, request, label)
      IMPLICIT NONE
      complex(C_DOUBLE_COMPLEX), intent(out) :: recv(:)
      complex(C_DOUBLE_COMPLEX), intent(in)  :: send(:)
      type(MPI_Request), intent(inout) :: request
      character(len=*), intent(in), optional :: label
      character(len=96) :: range_name

      range_name = "MPI_Ialltoall fft_transpose"
      if (present(label)) range_name = "MPI_Ialltoall "//trim(label)
      call roctxPush(range_name)
#ifndef HAVE_HIP
      !$omp target data use_device_addr(send, recv)
#endif
      call MPI_IALLTOALL(send, sendcount, MPI_DOUBLE_COMPLEX, &
                         recv, sendcount, MPI_DOUBLE_COMPLEX, MPI_COMM_X, request, ierr)
#ifndef HAVE_HIP
      !$omp end target data
#endif
      call roctxPop(range_name)

    END SUBROUTINE alltoall

    SUBROUTINE gather_full_y_line(ny, local_line, full_line)
      IMPLICIT NONE
      integer(C_INT), intent(in) :: ny
      complex(C_DOUBLE_COMPLEX), intent(in) :: local_line(ny0 - 2:nyN + 2)
      complex(C_DOUBLE_COMPLEX), intent(out) :: full_line(-1:ny + 1)

#ifdef HAVE_MPI
      complex(C_DOUBLE_COMPLEX) :: send_line(-1:ny + 1)
      integer(C_INT) :: send_start, send_end
#else
      integer(C_INT) :: send_start, send_end
#endif

      send_start = ny0
      send_end = nyN

      if (ipy == 0) send_start = -1
      if (ipy == npy_grid - 1) send_end = ny + 1

#ifdef HAVE_MPI
      send_line = (0.0d0, 0.0d0)

      send_line(send_start:send_end) = local_line(send_start:send_end)

      call roctxPush("MPI_Allreduce gather_full_y_line")
      call MPI_Allreduce(send_line, full_line, ny + 3, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_Y, ierr)
      call roctxPop("MPI_Allreduce gather_full_y_line")
#else
      full_line(-1:ny + 1) = local_line(-1:ny + 1)
#endif
    END SUBROUTINE gather_full_y_line

    SUBROUTINE allgather_y_device_complex_rows(send_rows, recv_rows, nrows, ncols, marker)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: send_rows(:, :)
      complex(C_DOUBLE_COMPLEX), intent(inout) :: recv_rows(:, :, :)
      integer(C_INT), intent(in) :: nrows, ncols
      character(len=*), intent(in) :: marker
#ifndef HAVE_MPI
      integer(C_INT) :: irow, icol
#endif

#ifdef HAVE_MPI
      call roctxPush(marker)
#ifndef HAVE_HIP
      !$omp target data use_device_addr(send_rows, recv_rows)
#endif
      call MPI_Allgather(send_rows, nrows*ncols, MPI_DOUBLE_COMPLEX, recv_rows, nrows*ncols, MPI_DOUBLE_COMPLEX, MPI_COMM_Y, ierr)
#ifndef HAVE_HIP
      !$omp end target data
#endif
      call roctxPop(marker)
#else
      !$omp target teams distribute parallel do collapse(2) default(none) &
      !$omp shared(send_rows, recv_rows, nrows, ncols) private(irow, icol)
      do icol = 1, ncols
        do irow = 1, nrows
          recv_rows(irow, icol, 1) = send_rows(irow, icol)
        end do
      end do
      !$omp end target teams distribute parallel do
#endif
    END SUBROUTINE allgather_y_device_complex_rows

    !------- Divide the problem in pencils ---------!
    !-----------------------------------------------!
    SUBROUTINE init_MPI(nxpp, nz, ny, nzd, nPhi, overlapping, npy_requested)
      integer(C_INT), intent(in)  :: nxpp, nz, ny, nzd, nPhi, npy_requested
      logical, intent(in) :: overlapping
      integer, parameter :: ndims = 4
      integer :: i, color, key
      integer :: array_of_sizes(ndims), array_of_subsizes(ndims), array_of_starts(ndims), ierror
      integer(C_INT) :: miny_local, maxy_local
#if defined(HAVE_HIP)
      type(c_ptr) :: sendptr, recvptr
#endif
      integer(c_size_t) :: sendsize, recvsize
      ! Define which process write on screen
      has_terminal = (iproc == 0)
      npy_grid = npy_requested
#ifdef HAVE_MPI
      if (npy_grid < 1) then
        if (has_terminal) print *, "Error: npy must be >= 1."
        call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
      end if
      if (mod(nproc, npy_grid) /= 0) then
        if (has_terminal) then
          print *, "Error: nproc must be divisible by npy."
          print *, "       Received nproc=", nproc, " npy=", npy_grid
        end if
        call MPI_Abort(MPI_COMM_WORLD, 1, ierror)
      end if
#else
      if (npy_grid /= 1) error stop "init_MPI: npy > 1 requires MPI"
#endif
      npxz = nproc/npy_grid
      ipy = iproc/npxz
      ipxz = mod(iproc, npxz)

#ifdef HAVE_MPI
      color = ipy
      key = ipxz
      call MPI_Comm_split(MPI_COMM_WORLD, color, key, MPI_COMM_X, ierr)
      color = ipxz
      key = ipy
      call MPI_Comm_split(MPI_COMM_WORLD, color, key, MPI_COMM_Y, ierr)
#endif
      ! Calculate domain division in wall-normal direction.
      ! ny0:nyN are the only persistent y-partition variables.
      ny0 = 1 + ipy*(ny - 1)/npy_grid
      nyN = (ipy + 1)*(ny - 1)/npy_grid

      ! Local write extent, including physical ghost rows only on the end ranks.
      miny_local = ny0
      maxy_local = nyN
      if (ipy == 0) miny_local = ny0 - 2
      if (ipy == npy_grid - 1) maxy_local = nyN + 2
      if (npy_grid == 1) then
        miny_local = ny0 - 2
        maxy_local = nyN + 2
      end if

      ! Calculate domain division
      nx0 = ipxz*(nxpp)/npxz; nxN = (ipxz + 1)*(nxpp)/npxz - 1; nxB = nxN - nx0 + 1; 
      nz0 = ipxz*nzd/npxz; nzN = (ipxz + 1)*nzd/npxz - 1; nzB = nzN - nz0 + 1; 
      call setup_y_decomposition(ny, 2*nz + 1, nxB)
      has_average = (nx0 == 0)
      !$omp target update to(npy_grid, npxz, ipy, ipxz, nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN)
      fft_transpose_is_local = (nzB == nzd)
#ifdef HAVE_MPI
#ifdef mpiverbose
      DO i = 0, nproc - 1
       IF (iproc==i) WRITE(*,*) "iproc=",iproc," ipxz=",ipxz," ipy=",ipy," nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, "ny0=", ny0, "nyN=", nyN
        CALL MPI_Barrier(MPI_COMM_WORLD)
      END DO
      FLUSH (output_unit)
#endif
      ! The pairwise all-to-all transpose uses one shared sendcount for every rank,
      ! so both decomposed dimensions must divide evenly across MPI ranks.
      if (mod(nxpp, npxz) /= 0 .or. mod(nzd, npxz) /= 0) then
        if (has_terminal) then
          print *, "Error: FFT transpose requires npxz to divide both nx+1 and nzd."
          print *, "       Received nx+1=", nxpp, " nzd=", nzd, " npxz=", npxz, " npy=", npy_grid
          print *, "       This run would create uneven transpose counts and can fail in MPI_Ialltoall."
        end if
        CALL MPI_Abort(MPI_COMM_WORLD, 1, ierror)
      end if
      if (int(npxz, 8)*int(nxB, 8)*int(nzB, 8)*int(nyN - ny0 + 5, 8) > huge(0_C_INT)) then
        if (has_terminal) then
          print *, "Error: problem too large for MPI transpose (integer overflow). Try to increase the number of processes."
        end if
        CALL MPI_Abort(MPI_COMM_WORLD, 1, ierror)
      end if
      sendsize = npxz*nxB*nzB*(nyN - ny0 + 5)
      recvsize = npxz*nxB*nzB*(nyN - ny0 + 5)

      sendcount = nxB*nzB*(nyN - ny0 + 5)
      !$omp target update to(sendcount)

      ! Allocate buffers for transposes*int(16, c_size_t)
#if defined(HAVE_HIP)
      ! On HIP with HSA_XNACK=1, the use_device_ptr statements around the MPI calls are ignored.
      ! Hence, MPI does a CPU mpi copy! So we need to allocate it explicity on the device.
      sendptr = omp_target_alloc(sendsize*int(16*merge(2, 1, overlapping), c_size_t), omp_get_default_device())
      recvptr = omp_target_alloc(recvsize*int(16*merge(2, 1, overlapping), c_size_t), omp_get_default_device())
      call c_f_pointer(sendptr, sendbuf, [sendsize, int(merge(2, 1, overlapping), c_size_t)])
      call c_f_pointer(recvptr, recvbuf, [recvsize, int(merge(2, 1, overlapping), c_size_t)])
#else
      ALLOCATE (sendbuf(sendsize, merge(2, 1, overlapping))); sendbuf = 0
      ALLOCATE (recvbuf(recvsize, merge(2, 1, overlapping))); recvbuf = 0
      !$omp target enter data map(alloc: sendbuf, recvbuf)
#endif

      ! For READING VELOCITY, SETTING VIEW: datatype that maps velocity on disk to memory (it differs from writing: halo cells are read twice!)
    CALL MPI_Type_create_subarray(ndims, [ny+3, 2*nz+1, nxpp, 3+nPhi], [nyN-ny0+5, 2*nz+1, nxB, 3+nPhi], [ny0-1,0,nx0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_read_type, ierror)
      CALL MPI_Type_commit(vel_read_type, ierror)
      ! For READING VELOCITY, datatype in memory (avoids overflow) - NOTICE THAT THIS ARRAY IS FULL (NOT REALLY SUBARRAY)
    CALL MPI_Type_create_subarray(ndims, [nyN-ny0+5, 2*nz+1, nxB, 3+nPhi], [nyN-ny0+5, 2*nz+1, nxB, 3+nPhi], [0,0,0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_field_type, ierror)
      CALL MPI_Type_commit(vel_field_type, ierror)
      ! For WRITING VELOCITY, SETTING VIEW: datatype to map distributed velocity array to file
      array_of_sizes = [ny + 3, 2*nz + 1, nxpp, 3 + nPhi] ! size along each dimension of the WHOLE array ON DISK
      array_of_subsizes = [maxy_local - miny_local + 1, 2*nz + 1, nxB, 3 + nPhi] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
      array_of_starts = [miny_local + 1, 0, nx0, 0] ! starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, writeview_type, ierror)
      CALL MPI_Type_commit(writeview_type, ierror)
      ! For WRITING VELOCITY, SKIPPING HALO CELLS: datatype with holes to skip halo cells and select only data to be written
      array_of_sizes = [(nyN + 2) - (ny0 - 2) + 1, 2*nz + 1, nxB, 3 + nPhi] ! size along each dimension of the array IN MEMORY owned by each process
      array_of_subsizes = [maxy_local - miny_local + 1, 2*nz + 1, nxB, 3 + nPhi] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
      array_of_starts = [miny_local - (ny0 - 2), 0, 0, 0] ! starting position of each component; !!! IT'S ZERO BASED AND WRT TO ARRAY IN MEMORY !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, owned2write_type, ierror)
      CALL MPI_Type_commit(owned2write_type, ierror)
#endif
    END SUBROUTINE init_MPI

  END MODULE mpi_transpose
