! TODO: write into files as append + flag to remove average or 0 mode

#include "header.h"

module runtime
    use dnsdata
    use mpi_transpose
    use ffts
    implicit none
#include HEADER_RUNTIME
    integer :: llz ! largest z index (iz) belonging to large scales -> large scales are -llz:llz
    real(C_DOUBLE) :: lc ! threshold wavelength for filter
    real(C_DOUBLE) :: ys, yl ! requested y positions for large and small scale signals
    logical :: remove_spatial_average ! flag to remove spatial or actual average from ul
    real(C_DOUBLE) :: ubar_yl ! value of average flow in x direction at yl
    integer :: iys, iyl ! indexes of nearest y positions to the requested ones
    integer :: iv ! a counter
    integer, parameter :: fil=777, fis=778
    logical :: has_ys = .FALSE., has_yl = .FALSE.
#ifdef nonblockingXZ
     TYPE(MPI_REQUEST) :: req_rt(1:6)
#endif
contains



    subroutine runtime_setup()
    ! this subroutine is taken from a file in folder "runtime"
    ! the filename is contained in the macro RUNTIME_SETUP_SUBROUTINE
    ! such macro is set in header.h
#include RUNTIME_SETUP_SUBROUTINE

        open(15, file='instabudget.in')
            read(15, *) lc
            read(15, *) ys, yl
            read(15, *) remove_spatial_average, ubar_yl
        close(15)

        ! process input
        llz=FLOOR(2*PI/lc/beta0)
        call get_nearest_iy(ys, iys, has_ys)
        call get_nearest_iy(yl, iyl, has_yl)

        ! write out settings
        if (has_terminal) then
            print *
            print *, "CALCULATION OF STATISTICS AT RUNTIME"
            print *, "INSTANTANEOUS LES BUDGET FOR SMALL SCALES (AM LIKE)"
            print *
            write(*,"(A,F11.6)") "   requested ys:", ys
            write(*,"(A,F11.6)") "   requested yl:", yl
            print *
        end if
        if (has_ys .AND. has_average) write(*,"(A,I6,A,F11.6)") "From process", iproc, " - actual ys:", y(iys)
        if (has_yl .AND. has_average) write(*,"(A,I6,A,F11.6)") "From process", iproc, " - actual yl:", y(iyl)
        if (has_terminal) then
            print *
            write(*,"(A,F11.6)") "   requested cutoff spanwise wavelength:", lc
            write(*,"(A,I5)") "   cutoff spanwise index:", llz
            print *
        end if

        ! open file streams for output
        if (has_ys .AND. has_terminal) open(unit=fis,file="us.out",access="stream",action="write")
        if (has_yl .AND. has_terminal) open(unit=fil,file="ul.out",access="stream",action="write")

    end subroutine



    subroutine runtime_finalise()
    ! this subroutine is taken from a file in folder "runtime"
    ! the filename is contained in the macro RUNTIME_FINALISE_SUBROUTINE
    ! such macro is set in header.h
#include RUNTIME_FINALISE_SUBROUTINE
        if (has_ys .AND. has_terminal) close(fis)
        if (has_yl .AND. has_terminal) close(fil)
    end subroutine



    subroutine runtime_save()
    ! this subroutine is taken from a file in folder "runtime"
    ! the filename is contained in the macro RUNTIME_SAVE_SUBROUTINE
    ! such macro is set in header.h
#include RUNTIME_SAVE_SUBROUTINE
    integer, parameter :: comp = 1

    ! large scales
    if (has_yl) then

        VVdz(:,:,1,1) = 0 ! clean up VVdz first
        ! copy velocity to FFTW array
        VVdz(1:(llz+1),1:nxB,1,1)=V(iyl,0:llz,nx0:nxN,comp);
        VVdz((nzd+1-llz):nzd,1:nxB,1,1)=V(iyl,(-llz):(-1),nx0:nxN,comp);

        ! remove average
        if (remove_spatial_average .AND. has_average) VVdz(1,1,1,1) = 0

        ! do the transform
        call IFT(VVdz(1:nzd,1:nxB,1,1))  
#ifdef nonblockingXZ
        call MPI_IAlltoall(VVdz(:,:,1,1), 1, Mdz, VVdx(:,:,1,1), 1, Mdx, MPI_COMM_X, req_rt(1))
#endif
#ifndef nonblockingXZ
        call MPI_Alltoall(VVdz(:,:,1,1), 1, Mdz, VVdx(:,:,1,1), 1, Mdx, MPI_COMM_X)
        VVdx(nx+2:nxd+1,1:nzB,1,1)=0;    call RFT(VVdx(1:nxd+1,1:nzB,1,1),rVVdx(1:2*nxd+2,1:nzB,1,1))
#endif
#ifdef nonblockingXZ
        call MPI_wait(req_rt(1),MPI_STATUS_IGNORE); 
        VVdx(nx+2:nxd+1,1:nzB,1,1)=0; 
        call RFT(VVdx(1:nxd+1,1:nzB,1,1),rVVdx(1:2*nxd+2,1:nzB,1,1));
#endif

        ! save to disk
        if (has_average) write(fil) rVVdx(1,1,1,1)

    end if

        ! small scales
        if (has_ys) then

            VVdz(:,:,2,1) = 0 ! clean up VVdz first
            ! copy velocity
            VVdz((llz+2):(nz+1),1:nxB,2,1)=V(iys,(llz+1):nz,nx0:nxN,comp);
            VVdz((nzd+1-nz):(nzd-llz),1:nxB,2,1)=V(iys,(-nz):(-llz-1),nx0:nxN,comp);

            ! do the transform
            call IFT(VVdz(1:nzd,1:nxB,2,1))  
#ifdef nonblockingXZ
            call MPI_IAlltoall(VVdz(:,:,2,1), 1, Mdz, VVdx(:,:,2,1), 1, Mdx, MPI_COMM_X, req_rt(2))
#endif
#ifndef nonblockingXZ
            call MPI_Alltoall(VVdz(:,:,2,1), 1, Mdz, VVdx(:,:,2,1), 1, Mdx, MPI_COMM_X)
            VVdx(nx+2:nxd+1,1:nzB,2,1)=0;    CALL RFT(VVdx(1:nxd+1,1:nzB,2,1),rVVdx(1:2*nxd+2,1:nzB,2,1))
#endif
#ifdef nonblockingXZ
            call MPI_wait(req_rt(2),MPI_STATUS_IGNORE); 
            VVdx(nx+2:nxd+1,1:nzB,2,1)=0; 
            call RFT(VVdx(1:nxd+1,1:nzB,2,1),rVVdx(1:2*nxd+2,1:nzB,2,1));
#endif

            ! save to disk
            if (has_average) write(fis) rVVdx(1,1,2,1)*rVVdx(1,1,2,1)

        end if

        ! synchronise
        call MPI_Barrier(MPI_COMM_WORLD)
    
    end subroutine



#include RUNTIME_AUXILIARY_SUBROUTINES



    subroutine get_nearest_iy(reqy, iyy, has_yy)
    integer, intent(out) :: iyy
    logical, intent(out) :: has_yy
    real(C_DOUBLE), intent(in) :: reqy
    real(C_DOUBLE) :: cdiff
    integer :: ii

        iyy = ny0-2; cdiff = abs(reqy - y(iyy))
        do ii = ny0-2, nyN+2
            if (abs(reqy - y(ii)) < cdiff) then
                iyy = ii
                cdiff = abs(reqy - y(iyy))
            end if
        end do

        if ( (iyy >= miny) .AND. (iyy <= maxy) ) then
            has_yy = .TRUE.
        end if
        
    end subroutine



end module