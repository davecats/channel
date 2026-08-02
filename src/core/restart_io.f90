#include "header.h"

! Reading and writing restart files.
!
! The on-disk format is a small header -- nx, ny, nz then alfa0, beta0, ni, a,
! ymin, ymax, time -- followed by the velocity and scalar field written through
! an MPI-IO subarray view.  The reader validates the header against the current
! configuration and stops on a mismatch.
!
! Everything the header needs beyond the grid (ni, time) is passed in rather
! than read from the DNS state module, so that this module sits below it and
! the dependency runs one way.
module restart_io

  use, intrinsic :: iso_c_binding
  use channel_grid
  use roctx, only: roctxPush, roctxPop
#ifdef HAVE_MPI
  use mpi_transpose, only: vel_read_type, vel_field_type, writeview_type, owned2write_type
  use mpi_f08
#endif

  implicit none
  private

  public :: restart_read, restart_write

contains

  ! Reads a restart file into R, or fills R with a perturbed laminar profile
  ! when the file is absent.  time is taken from the file header.
  SUBROUTINE restart_read(filename, R, ni, time, perturbation_amplitude)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: ni, perturbation_amplitude
    real(C_DOUBLE), intent(inout) :: time
    complex(C_DOUBLE_COMPLEX), intent(INOUT) :: R(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3 + nPhi)
    character(len=*), intent(IN) :: filename
    integer(C_SIZE_T) :: ix, iy, iz, io, iPhi
    integer(C_INT) :: r_nx, r_ny, r_nz
    real(C_DOUBLE) :: r_alfa0, r_beta0, r_ni, r_a, r_ymin, r_ymax
    real(C_DOUBLE) :: rn(1:3)
#ifdef HAVE_MPI
    INTEGER(MPI_OFFSET_KIND) :: disp = 3*C_INT + 7*C_DOUBLE
    TYPE(MPI_File) :: fh

    OPEN (UNIT=120, FILE=TRIM(filename), access="stream", status="old", action="read", iostat=io)
    IF (io == 0) THEN
      if (has_terminal) print *, "Reading from file "//filename
      READ (120, POS=1) r_nx, r_ny, r_nz, r_alfa0, r_beta0, r_ni, r_a, r_ymin, r_ymax, time
      call MPI_file_open(MPI_COMM_WORLD, TRIM(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh)
      call MPI_file_set_view(fh, disp, MPI_DOUBLE_COMPLEX, vel_read_type, 'native', MPI_INFO_NULL)
      call roctxPush("MPI_File_read_all restart")
      call MPI_file_read_all(fh, R, 1, vel_field_type, MPI_STATUS_IGNORE)
      call roctxPop("MPI_File_read_all restart")
      call MPI_file_close(fh)
      IF (r_nx /= nx .OR. r_ny /= ny .OR. r_nz /= nz .OR. r_alfa0 /= alfa0 .OR. r_beta0 /= beta0 .OR. r_ni /= ni .OR. r_a /= a .OR. r_ymin /= ymin .OR. r_ymax /= ymax) THEN
        IF (has_terminal) THEN
          PRINT *, "ERROR: mismatch in metadata between restart file and dns.in. Stopping."
          PRINT *, "From .out file:"
          PRINT *, r_nx, r_ny, r_nz, r_alfa0, r_beta0, r_ni, r_a, r_ymin, r_ymax
          PRINT *, "From dns.in:"
          PRINT *, nx, ny, nz, alfa0, beta0, ni, a, ymin, ymax
        END IF
        STOP
      END IF
    ELSE
#endif
      IF (has_terminal) PRINT *, "Restart file "//filename//" not found"
      R = 0
      IF (has_terminal) WRITE (*, *) "Generating initial field..."
      DO iy = ny0 - 2, nyN + 2; DO ix = nx0, nxN; DO iz = -nz, nz
          CALL RANDOM_NUMBER(rn)
          R(iy, iz, ix, 1) = perturbation_amplitude*EXP(dcmplx(0, rn(1) - 0.5))
          R(iy, iz, ix, 2) = perturbation_amplitude*EXP(dcmplx(0, rn(2) - 0.5))
          R(iy, iz, ix, 3) = perturbation_amplitude*EXP(dcmplx(0, rn(3) - 0.5))
          !!R(iy,iz,ix,1) = 0.0001*EXP(dcmplx(0,rn(1)-0.5));  R(iy,iz,ix,2) = 0.0001*EXP(dcmplx(0,rn(2)-0.5));  R(iy,iz,ix,3) = 0.0001*EXP(dcmplx(0,rn(3)-0.5));
        END DO; END DO; END DO
      IF (has_average) THEN
        DO iy = ny0 - 2, nyN + 2
          R(iy, 0, 0, 1) = 3*0.5*y(iy)*(2 - y(iy))
          !R(iy, 0, 0, 1) = 3*0.5*y(iy)*(2 - y(iy)) + 0.01*SIN(8*y(iy)*2*PI)/ni
          !R(iy, 0, 0, 1) = y(iy)*(2 - y(iy))*3.d0/2.d0 + 0.001*SIN(8*y(iy)*2*PI);
          !V(iy,0,0,1)=y(iy)-1
          DO iPhi = 1, nPhi
            R(iy, 0, 0, 3 + iPhi) = 3*0.5*y(iy)*(2 - y(iy))
          END DO
        END DO
      END IF
#ifdef HAVE_MPI
    END IF
#endif
    CLOSE (120)
  END SUBROUTINE restart_read


  ! Writes R plus the header the reader validates against.  The caller is
  ! responsible for having the field on the host.
  SUBROUTINE restart_write(filename, R, ni, time)
    IMPLICIT NONE
    real(C_DOUBLE), intent(in) :: ni, time
    complex(C_DOUBLE_COMPLEX), intent(in) :: R(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN, 1:3 + nPhi)
    character(len=*), intent(in) :: filename
    ! mpi stuff
#ifdef HAVE_MPI
    TYPE(MPI_File) :: fh
    INTEGER(MPI_OFFSET_KIND) :: disp
    TYPE(MPI_Status) :: status

    ! open file
    CALL MPI_File_open(MPI_COMM_WORLD, TRIM(filename), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)

    ! write header
    IF (has_terminal) THEN ! only one process does this
      CALL MPI_file_write(fh, [nx, ny, nz], 3, MPI_INTEGER, status)
      CALL MPI_file_write(fh, [alfa0, beta0, ni, a, ymin, ymax, time], 7, MPI_DOUBLE_PRECISION, status)
    END IF

    ! set view to subarray
    disp = 3*C_INT + 7*C_DOUBLE ! offset to skip header
    CALL MPI_File_set_view(fh, disp, MPI_DOUBLE_COMPLEX, writeview_type, 'native', MPI_INFO_NULL)

    ! finally write field
    call roctxPush("MPI_File_write_all restart")
    CALL MPI_File_write_all(fh, R, 1, owned2write_type, status)
    call roctxPop("MPI_File_write_all restart")

    ! close file
    call MPI_File_close(fh)
#endif
  END SUBROUTINE restart_write

end module restart_io
