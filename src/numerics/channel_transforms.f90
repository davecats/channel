#include "build_options.h"

! The spectral <-> physical round trip each Runge-Kutta substep runs, and the
! CFL number read off the physical-space velocities while they exist.
!
! transform_to_physical takes V forward: assemble the z-padded buffer, inverse
! transform in z, transpose z->x, zero the x padding, inverse transform in x.
! transform_back_and_build_rhs takes the nonlinear products the other way and
! accumulates the right-hand side.  Both walk one component per iteration and,
! when overlapping is on, run one iteration ahead so the alltoall of component
! m overlaps the transforms of component m-1 -- which is why the loop bound and
! the two buffer indices are MERGE expressions rather than constants.
!
! compute_cfl lives here because it consumes rVVdx inside the same fft
! workspace window as the forward transform; the cross-rank reduction of what
! it leaves in cfl, and the timestep adaptation, are in statistics.

module channel_transforms

  use, intrinsic :: iso_c_binding
  use channel_grid
  use channel_state
  use channel_equations, only: build_products, buildrhs
  use mpi_transpose
  use roctx, only: roctxPush, roctxPop
  use ffts

  IMPLICIT NONE

CONTAINS

  SUBROUTINE assemble_vvdz(m, to)
    IMPLICIT NONE
    integer(C_INT) :: i, j, k, y_first, y_last
    integer(C_INT), intent(in) :: m, to
    y_first = ny0
    y_last = nyN

    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, VVdz) shared(nxB, nzd, nx0, nxN, y_first, y_last, nz, m, to) private(i,j,k)
    DO i = y_first - 2, y_last + 2
      DO j = 1, nxB
        DO k = 1, nzd
          VVdz(k, j, i, to) = 0.0d0
        END DO
      END DO
    END DO
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(V, VVdz) shared(nx0, nxN, nz, nzd, y_first, y_last, m, to) private(i,j,k)
    DO i = y_first - 2, y_last + 2
      DO j = nx0, nxN
        DO k = 1, nzd
          IF (k <= nz + 1) THEN
            VVdz(k, j - nx0 + 1, i, to) = V(i, k - 1, j, m)
          ELSEIF (k >= nz + 2 .AND. k <= nzd - nz) THEN
            VVdz(k, j - nx0 + 1, i, to) = 0.0d0
          ELSE
            VVdz(k, j - nx0 + 1, i, to) = V(i, k - nzd - 1, j, m)
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE assemble_vvdz

  SUBROUTINE zero_vvdx_hft(to)
    IMPLICIT NONE
    integer(C_INT) :: i, j, k, y_first, y_last
    integer(C_INT), intent(in) :: to
    y_first = ny0
    y_last = nyN
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdx, nx, nxd, nzB, y_first, y_last, to) private(i,j,k)
    DO i = y_first - 2, y_last + 2
      DO j = 1, nzB
        DO k = nx + 2, nxd + 1
          VVdx(k, j, i, to) = 0.0
        END DO
      END DO
    END DO
  END SUBROUTINE zero_vvdx_hft

  subroutine compute_cfl()
    implicit none
    integer(C_INT) :: i, j, k, y_first, y_last
    real(C_DOUBLE) :: tmp
    y_first = ny0
    y_last = nyN
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp private(i,j,k,tmp) reduction(max:cfl) &
    !$omp shared(rVVdx, dx, dy, dz, ny, y_first, y_last, nxd, nzB)
    do j = 1, 2*nxd
      do k = 1, nzB
        do i = max(y_first - 2, 1_C_INT), min(y_last + 2, ny - 1)
          tmp = abs(rVVdx(j, k, i, 1))/dx + abs(rVVdx(j, k, i, 2))/dy(i) + abs(rVVdx(j, k, i, 3))/dz
          cfl = max(cfl, tmp)
        end do
      end do
    end do
  END SUBROUTINE compute_cfl

  SUBROUTINE transform_to_physical()
    IMPLICIT NONE
    integer(C_INT) ::  m, to, from, mm1
    type(MPI_Request), dimension(:) :: requests(3 + nPhi)
    type(MPI_Status)  :: status
    DO m = 1, MERGE(3 + nPhi + 1, 3 + nPhi, overlapping)

      ! Compute indices depending on overlap mode
      mm1 = MERGE(m - 1, m, overlapping)
      to = MERGE(mod(m - 1, 2) + 1, 1, overlapping)
      from = MERGE(mod(mm1 - 1, 2) + 1, 1, overlapping)

      ! Step 1: assemble, pack, post alltoall (only if in range)
      if (m <= 3 + nPhi) then
        call roctxPush("transform_to_physical assemble_vvdz")
        CALL assemble_vvdz(m, to)
        call roctxPop("transform_to_physical assemble_vvdz")
        call roctxPush("transform_to_physical IFT")
        CALL IFT(VVdz(:, :, :, to))
        call roctxPop("transform_to_physical IFT")
        if (fft_transpose_is_local) then
          call roctxPush("transform_to_physical repack_zTOx_local")
          call repack_zTOx_local(VVdz(:, :, :, to), VVdx(:, :, :, to))
          call roctxPop("transform_to_physical repack_zTOx_local")
        else
          call roctxPush("transform_to_physical pack_zTOx")
          CALL pack_zTOx(VVdz(:, :, :, to), sendbuf(:, to))
          call roctxPop("transform_to_physical pack_zTOx")
          CALL alltoall(sendbuf(:, to), recvbuf(:, to), requests(m), "zTOx transform_to_physical")
        end if
      end if

      ! Step 2: wait, unpack, FFT (depending on overlap)
      if (MERGE(m > 1, .true., overlapping)) then
        if (.not. fft_transpose_is_local) then
          call roctxPush("MPI_Wait zTOx transform_to_physical")
          CALL MPI_WAIT(requests(mm1), status, ierr)
          call roctxPop("MPI_Wait zTOx transform_to_physical")
          call roctxPush("transform_to_physical unpack_zTOx")
          CALL unpack_zTOx(recvbuf(:, from), VVdx(:, :, :, from))
          call roctxPop("transform_to_physical unpack_zTOx")
        end if
        call roctxPush("transform_to_physical zero_vvdx_hft")
        CALL zero_vvdx_hft(from)
        call roctxPop("transform_to_physical zero_vvdx_hft")
        call roctxPush("transform_to_physical RFT")
        CALL RFT(VVdx(:, :, :, from), rVVdx(:, :, :, mm1))
        call roctxPop("transform_to_physical RFT")
      end if
    END DO
  END SUBROUTINE transform_to_physical

  SUBROUTINE transform_back_and_build_rhs(ODE)
    IMPLICIT NONE
    integer(C_INT) ::  m, mm1, to, from
    real(C_DOUBLE), intent(in) :: ODE(1:3)
    type(MPI_Request), dimension(:) :: requests(6 + 3*nPhi)
    type(MPI_Status)  :: status
    ! Reverse pass to build rVVdx
    DO m = 1, MERGE(6 + 3*nPhi + 1, 6 + 3*nPhi, overlapping)

      ! Buffer indices depend on overlap mode

      mm1 = MERGE(m - 1, m, overlapping)
      to = MERGE(mod(m - 1, 2) + 1, 1, overlapping)
      from = MERGE(mod(m, 2) + 1, 1, overlapping)

      ! Step 1: Build, HFT, pack, and post alltoall
      if (m <= 6 + 3*nPhi) then
        call roctxPush("transform_back build_products")
        call build_products(m, to)
        call roctxPop("transform_back build_products")
        call roctxPush("transform_back HFT")
        call HFT(products(:, :, :, to), VVdx(:, :, :, to))
        call roctxPop("transform_back HFT")
        if (fft_transpose_is_local) then
          call roctxPush("transform_back repack_xTOz_local")
          call repack_xTOz_local(VVdx(:, :, :, to), VVdz(:, :, :, to))
          call roctxPop("transform_back repack_xTOz_local")
        else
          call roctxPush("transform_back pack_xTOz")
          call pack_xTOz(VVdx(:, :, :, to), sendbuf(:, to))
          call roctxPop("transform_back pack_xTOz")
          call alltoall(sendbuf(:, to), recvbuf(:, to), requests(m), "xTOz transform_back_and_build_rhs")
        end if
      end if

      ! Step 2: Wait, unpack, FFT, and build RHS
      if (MERGE(m > 1, .true., overlapping)) then
        if (.not. fft_transpose_is_local) then
          call roctxPush("MPI_Wait xTOz transform_back_and_build_rhs")
          call MPI_WAIT(requests(mm1), status, ierr)
          call roctxPop("MPI_Wait xTOz transform_back_and_build_rhs")
          call roctxPush("transform_back unpack_xTOz")
          call unpack_xTOz(recvbuf(:, from), VVdz(:, :, :, from))
          call roctxPop("transform_back unpack_xTOz")
        end if
        call roctxPush("transform_back FFT")
        call FFT(VVdz(:, :, :, from))
        call roctxPop("transform_back FFT")
        call roctxPush("transform_back buildrhs")
        call buildrhs(ODE, mm1, from)
        call roctxPop("transform_back buildrhs")
      end if
    END DO
  END SUBROUTINE transform_back_and_build_rhs

  ! The same round trip for a single field, without the pipelining above: one
  ! field is already sitting in VVdz(:,:,:,1), so there is nothing to overlap
  ! the alltoall with and the four steps stay in sequence.  convvelo and
  ! pressure_output both drive it -- they differ only in what they load into
  ! the buffer beforehand and what they do with the result -- so the label
  ! argument is what shows up in the profile.
  SUBROUTINE spectral_field_to_real_x(rx, label)
    IMPLICIT NONE
    real(C_DOUBLE), intent(out) :: rx(:, :, ny0 - 2:)
    character(len=*), intent(in) :: label

    call IFT(VVdz(:, :, :, 1))
    call transpose_zTOx(VVdz(:, :, :, 1), VVdx(:, :, :, 1), label)
    call zero_vvdx_hft(1)
    call RFT(VVdx(:, :, :, 1), rx)
  END SUBROUTINE spectral_field_to_real_x

  SUBROUTINE real_x_to_spectral_field(rx, field, label)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: rx(:, :, ny0 - 2:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN)
    character(len=*), intent(in) :: label
    integer(C_INT) :: ix, iz, iy, y_first, y_last
    y_first = ny0
    y_last = nyN

    call HFT(rx, VVdx(:, :, :, 1))
    call transpose_xTOz(VVdx(:, :, :, 1), VVdz(:, :, :, 1), label)
    call FFT(VVdz(:, :, :, 1))
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz, nx0, nxN, ny, y_first, y_last, nz, field) private(ix, iz, iy)
    do ix = nx0, nxN
      do iy = y_first - 2, y_last + 2
        do iz = 0, nz
          field(iy, iz, ix) = VVdz(iz + 1, ix - nx0 + 1, iy, 1)
        end do
      end do
    end do
    !$omp target teams distribute parallel do collapse(3) default(none) &
    !$omp shared(VVdz, nx0, nxN, ny, y_first, y_last, nz, field, izd) private(ix, iz, iy)
    do ix = nx0, nxN
      do iy = y_first - 2, y_last + 2
        do iz = -nz, -1
          field(iy, iz, ix) = VVdz(izd(iz) + 1, ix - nx0 + 1, iy, 1)
        end do
      end do
    end do
  END SUBROUTINE real_x_to_spectral_field
end module channel_transforms
