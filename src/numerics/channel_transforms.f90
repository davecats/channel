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
  use channel_equations, only: build_products, product_factors, buildrhs, buildrhs_finish
  use mpi_transpose
  use roctx, only: roctxPush, roctxPop
  use ffts

  IMPLICIT NONE

CONTAINS

  SUBROUTINE assemble_vvdz(m, to)
    IMPLICIT NONE
    integer(C_INT) :: i, j, k, y_first, y_last, ib, kb, ii, ki
    integer(C_INT), intent(in) :: m, to
    y_first = ny0
    y_last = nyN

    ! Every VVdz element is written below (the dealiasing band as zero), so no
    ! separate zeroing pass.  8x8 (k, i) tiles, k fastest: VVdz is z-fastest and V
    ! is y-fastest, so a wavefront writes full z-lines and reads full y-lines.
    !$omp target teams distribute parallel do collapse(5) default(none) &
    !$omp shared(V, VVdz) shared(nx0, nxN, nz, nzd, y_first, y_last, m, to) private(i, j, k, ib, kb, ii, ki)
    DO j = nx0, nxN
    DO ib = 0, (y_last - y_first + 4)/8
    DO kb = 0, (nzd - 1)/8
    DO ii = 0, 7
    DO ki = 0, 7
        i = y_first - 2 + ib*8 + ii
        k = kb*8 + ki + 1
        IF (i > y_last + 2 .OR. k > nzd) CYCLE
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
    ! x (j) innermost: rVVdx is x-fastest, so a wavefront reads contiguous lines.
    ! With y innermost every thread read three values 6 MB apart (20 ms -> ~1.5 ms).
    do i = max(y_first - 2, 1_C_INT), min(y_last + 2, ny - 1)
      do k = 1, nzB
        do j = 1, 2*nxd
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
    logical :: use_cb
    ! Local transposes: reading VVdz transposed in the RFT measured no faster
    ! than repack_zTOx_local, so the callbacks are for npxz > 1 only.
    use_cb = .not. fft_transpose_is_local
    if (use_cb) use_cb = xz_callbacks_enabled()
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
        if (use_cb .and. .not. fft_transpose_is_local) then
          ! IFT writes the send buffer itself (no pack kernel)
          call roctxPush("transform_to_physical IFT_pack")
          call IFT_pack(VVdz(:, :, :, to), c_loc(sendbuf(1, to)))
          call roctxPop("transform_to_physical IFT_pack")
          CALL alltoall(sendbuf(:, to), recvbuf(:, to), requests(m), "zTOx transform_to_physical", to)
        else
        call roctxPush("transform_to_physical IFT")
        CALL IFT(VVdz(:, :, :, to))
        call roctxPop("transform_to_physical IFT")
        if (use_cb) then
          ! local transpose: RFT_unpack reads VVdz directly
        else if (fft_transpose_is_local) then
          call roctxPush("transform_to_physical repack_zTOx_local")
          call repack_zTOx_local(VVdz(:, :, :, to), VVdx(:, :, :, to))
          call roctxPop("transform_to_physical repack_zTOx_local")
        else
          call roctxPush("transform_to_physical pack_zTOx")
          CALL pack_zTOx(VVdz(:, :, :, to), sendbuf(:, to))
          call roctxPop("transform_to_physical pack_zTOx")
          CALL alltoall(sendbuf(:, to), recvbuf(:, to), requests(m), "zTOx transform_to_physical", to)
        end if
        end if
      end if

      ! Step 2: wait, unpack, FFT (depending on overlap)
      if (MERGE(m > 1, .true., overlapping)) then
        if (.not. fft_transpose_is_local) then
          call roctxPush("MPI_Wait zTOx transform_to_physical")
          CALL alltoall_wait(requests(mm1), from)
          call roctxPop("MPI_Wait zTOx transform_to_physical")
        end if
        if (use_cb) then
          ! RFT reads the received blocks (or VVdz) and zero-pads x itself
          call roctxPush("transform_to_physical RFT_unpack")
          if (fft_transpose_is_local) then
            CALL RFT_unpack(c_loc(VVdz(1, 1, ny0 - 2, from)), VVdx(:, :, :, from), rVVdx(:, :, :, mm1))
          else
            CALL RFT_unpack(c_loc(recvbuf(1, from)), VVdx(:, :, :, from), rVVdx(:, :, :, mm1))
          end if
          call roctxPop("transform_to_physical RFT_unpack")
          cycle
        end if
        if (.not. fft_transpose_is_local) then
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
    integer(C_INT) ::  m, mm1, to, from, fa, fb
    real(C_DOUBLE), intent(in) :: ODE(1:3)
    type(MPI_Request), dimension(:) :: requests(6 + 3*nPhi)
    logical :: use_cb, use_cb_products
    type(c_ptr) :: xbuf
    use_cb = .not. fft_transpose_is_local
    if (use_cb) use_cb = xz_callbacks_enabled()
    use_cb_products = .false.
    if (use_cb) use_cb_products = xz_callbacks_products()
    ! Local transposes keep the compact x-pencil copy in products(:,:,:,to),
    ! which is then busy as the HFT input unless the callback forms the products.
    if (fft_transpose_is_local .and. .not. use_cb_products) use_cb = .false.
    ! Reverse pass to build rVVdx
    DO m = 1, MERGE(6 + 3*nPhi + 1, 6 + 3*nPhi, overlapping)

      ! Buffer indices depend on overlap mode

      mm1 = MERGE(m - 1, m, overlapping)
      to = MERGE(mod(m - 1, 2) + 1, 1, overlapping)
      from = MERGE(mod(m, 2) + 1, 1, overlapping)

      ! Step 1: Build, HFT, pack, and post alltoall
      if (m <= 6 + 3*nPhi .and. use_cb) then
        ! The HFT forms the product in its load callback and writes the send
        ! buffer (or, local, a compact copy in products(:,:,:,to)) itself.
        if (fft_transpose_is_local) then
          xbuf = c_loc(products(1, 1, ny0 - 2, to))
        else
          xbuf = c_loc(sendbuf(1, to))
        end if
        if (use_cb_products) then
          call product_factors(m, fa, fb)
          call roctxPush("transform_back HFT_products")
          call HFT_products(rVVdx, fa, fb, factor, VVdx(:, :, :, to), xbuf)
          call roctxPop("transform_back HFT_products")
        else
          call roctxPush("transform_back build_products")
          call build_products(m, to)
          call roctxPop("transform_back build_products")
          call roctxPush("transform_back HFT_pack")
          call HFT_pack(products(:, :, :, to), VVdx(:, :, :, to), xbuf)
          call roctxPop("transform_back HFT_pack")
        end if
        if (.not. fft_transpose_is_local) &
          call alltoall(sendbuf(:, to), recvbuf(:, to), requests(m), "xTOz transform_back_and_build_rhs", to)
      else if (m <= 6 + 3*nPhi) then
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
          call alltoall(sendbuf(:, to), recvbuf(:, to), requests(m), "xTOz transform_back_and_build_rhs", to)
        end if
      end if

      ! Step 2: Wait, unpack, FFT, and build RHS
      if (MERGE(m > 1, .true., overlapping)) then
        if (.not. fft_transpose_is_local) then
          call roctxPush("MPI_Wait xTOz transform_back_and_build_rhs")
          call alltoall_wait(requests(mm1), from)
          call roctxPop("MPI_Wait xTOz transform_back_and_build_rhs")
        end if
        if (use_cb) then
          call roctxPush("transform_back FFT_unpack")
          if (fft_transpose_is_local) then
            call FFT_unpack(c_loc(products(1, 1, ny0 - 2, from)), VVdz(:, :, :, from))
          else
            call FFT_unpack(c_loc(recvbuf(1, from)), VVdz(:, :, :, from))
          end if
          call roctxPop("transform_back FFT_unpack")
        else
        if (.not. fft_transpose_is_local) then
          call roctxPush("transform_back unpack_xTOz")
          call unpack_xTOz(recvbuf(:, from), VVdz(:, :, :, from))
          call roctxPop("transform_back unpack_xTOz")
        end if
        call roctxPush("transform_back FFT")
        call FFT(VVdz(:, :, :, from))
        call roctxPop("transform_back FFT")
        end if
        call roctxPush("transform_back buildrhs")
        call buildrhs(ODE, mm1, from)
        call roctxPop("transform_back buildrhs")
      end if
    END DO
    call roctxPush("transform_back buildrhs_finish")
    call buildrhs_finish(ODE)
    call roctxPop("transform_back buildrhs_finish")
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
