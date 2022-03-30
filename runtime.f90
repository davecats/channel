! problem
! need to take into account that processes do not have all y
! this is relevant in searching for iys/iyl, calculating transforms, writing out
! maybe put barrier at the end


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
    integer :: iys, iyl ! indexes of nearest y positions to the requested ones
    integer :: iv ! a counter
    integer, parameter :: fil=777, fis=778
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
        close(15)

        ! process input
        llz=FLOOR(2*PI/lc/beta0)
        call get_nearest_idx(ys, y, iys)
        call get_nearest_idx(yl, y, iyl)

        ! write out settings
        if (has_terminal) then
            print *, "CALCULATION OF STATISTICS AT RUNTIME"
            print *, "INSTANTANEOUS LES BUDGET FOR SMALL SCALES (AM LIKE)"
            print *
            write(*,"(A,I5)") "   requested ys:", ys
            write(*,"(A,I5)") "   requested yl:", yl
            print *
            write(*,"(A,I5)") "   actual ys:", y(iys)
            write(*,"(A,I5)") "   actual yl:", y(iyl)
            print *
            write(*,"(A,I5)") "   requested cutoff spanwise wavelength:", lc
            write(*,"(A,I5)") "   cutoff spanwise index:", llz
        end if

        ! open file streams for output
        open(unit=fil,file="ul.out",access="stream",action="write")
        open(unit=fis,file="us.out",access="stream",action="write")

    end subroutine



    subroutine runtime_finalise()
    ! this subroutine is taken from a file in folder "runtime"
    ! the filename is contained in the macro RUNTIME_FINALISE_SUBROUTINE
    ! such macro is set in header.h
#include RUNTIME_FINALISE_SUBROUTINE
        close(fil)
        close(fis)
    end subroutine



    subroutine runtime_save()
    ! this subroutine is taken from a file in folder "runtime"
    ! the filename is contained in the macro RUNTIME_SAVE_SUBROUTINE
    ! such macro is set in header.h
#include RUNTIME_SAVE_SUBROUTINE
    integer, parameter :: comp = 1

    ! begin by transforming velocity
    VVdz(:,:,1:2,1) = 0 ! clean up VVdz first
    ! large scale
    VVdz(1:(llz+1),1:nxB,1,1)=V(iyl,0:llz,nx0:nxN,comp);
    VVdz((nzd+1-llz):nzd,1:nxB,1,1)=V(iyl,(-llz):(-1),nx0:nxN,comp);
    ! small scale
    VVdz((llz+2):(nz+1),1:nxB,2,1)=V(iys,(llz+1):nz,nx0:nxN,comp);
    VVdz((nzd+1-nz):(nzd-llz),1:nxB,2,1)=V(iys,(-nz):(-llz-1),nx0:nxN,comp);

    ! do the transform
    DO iV=1,2
       CALL IFT(VVdz(1:nzd,1:nxB,iV,1))  
#ifdef nonblockingXZ
       CALL MPI_IAlltoall(VVdz(:,:,iV,1), 1, Mdz, VVdx(:,:,iV,1), 1, Mdx, MPI_COMM_X, req_rt(iV))
#endif
#ifndef nonblockingXZ
       CALL MPI_Alltoall(VVdz(:,:,iV,1), 1, Mdz, VVdx(:,:,iV,1), 1, Mdx, MPI_COMM_X)
       VVdx(nx+2:nxd+1,1:nzB,iV,1)=0;    CALL RFT(VVdx(1:nxd+1,1:nzB,iV,1),rVVdx(1:2*nxd+2,1:nzB,iV,1))
#endif
     END DO
#ifdef nonblockingXZ
     DO iV=1,3
       CALL MPI_wait(req_rt(iV),MPI_STATUS_IGNORE); 
       VVdx(nx+2:nxd+1,1:nzB,iV,1)=0; 
       CALL RFT(VVdx(1:nxd+1,1:nzB,iV,1),rVVdx(1:2*nxd+2,1:nzB,iV,1));
     END DO
#endif

    write(fil) rVVdx(1,1,1,1)
    write(fis) rVVdx(1,1,2,1)

    end subroutine



#include RUNTIME_AUXILIARY_SUBROUTINES



    subroutine get_nearest_idx(val, arr, idx)
    integer, intent(out) :: idx
    real(C_DOUBLE), intent(in) :: val, arr(*)
    real(C_DOUBLE), allocatable :: diff(:)
    real(C_DOUBLE) :: cv
    integer :: hi, lo, ii 

        ! retrieve upper and lower bounds for array
        hi = ubound(arr, 1)
        lo = lbound(arr, 1)

        ! allocate diff and calculate its values
        allocate(diff(lo:hi))
        diff = abs(arr - val)

        ! scan for minimum
        idx=lo; cv=diff(lo) 
        do ii = (lo+1), hi
            if (diff(ii) < cv) then
                idx = ii
                cv = diff(ii)
            end if
        end do
        
    end subroutine



end module