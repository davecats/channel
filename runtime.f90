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
    integer :: iv, is=1, il=1 ! a counter
    integer, parameter :: ft=777
    logical :: has_ys = .FALSE., has_yl = .FALSE.
    integer(MPI_OFFSET_KIND) :: no_reads = 0, fsize, disp
    real(C_DOUBLE) :: targtime

    type(MPI_COMM) :: MPI_COMM_L, MPI_COMM_S
    logical, allocatable :: gather_yl(:), gather_ys(:)
    integer :: snp=0, lnp=0 ! no of procs that have pieces of planes y=ys and y=yl
    integer, allocatable :: s_list(:), l_list(:) ! list of procs who have pieces of planes y=ys and y=yl
    type(MPI_GROUP) :: world_group, s_group, l_group
    
    type(MPI_File) :: fh_l, fh_s
    type(MPI_Datatype) :: yconst_plane_type
    integer(C_INT) :: ierror
    type(MPI_Status) :: status


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

        ! create communicators for yl and ys
        allocate(gather_yl(nproc))
        allocate(gather_ys(nproc))
        call MPI_ALLGATHER(has_yl, 1, MPI_LOGICAL, gather_yl, 1, MPI_LOGICAL, MPI_COMM_WORLD)
        call MPI_ALLGATHER(has_ys, 1, MPI_LOGICAL, gather_ys, 1, MPI_LOGICAL, MPI_COMM_WORLD)
        ! get how many procs have pieces of planes y=ys and y=yl
        do iv = 1,nproc 
            if (gather_yl(iv)) lnp = lnp + 1
            if (gather_ys(iv)) snp = snp + 1
        end do
        ! get list of procs who have pieces of planes y=ys and y=yl
        allocate(l_list(lnp))
        allocate(s_list(snp))
        do iv = 1,nproc 
            if (gather_yl(iv)) then
                l_list(il) = iv - 1
                il = il + 1
            end if
            if (gather_ys(iv)) then
                s_list(is) = iv - 1
                is = is + 1
            end if
        end do
        ! create groups for small and large
        call MPI_COMM_GROUP(MPI_COMM_WORLD, world_group)
        call MPI_GROUP_INCL(world_group, lnp, l_list, l_group)
        call MPI_GROUP_INCL(world_group, snp, s_list, s_group)
        ! now create comm
        call MPI_COMM_CREATE(MPI_COMM_WORLD, l_group, MPI_COMM_L)
        call MPI_COMM_CREATE(MPI_COMM_WORLD, s_group, MPI_COMM_S)

        ! create datatype for writing
        call MPI_Type_create_subarray(1, [nzd], [nzB], [nz0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, yconst_plane_type, ierror)
        call MPI_Type_commit(yconst_plane_type, ierror)

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
        if (has_terminal) then ! look for correct time to start writing
            open(unit=ft,file="time.out",access="stream",action="readwrite",position="rewind")
            if (time_from_restart) then
                inquire(ft, size=fsize)
                if (fsize/C_DOUBLE > 0) then
                    read(ft) targtime
                    no_reads = no_reads + 1
                end if
                do while (fsize/C_DOUBLE - no_reads > 0 .AND. targtime < time)
                    read(ft) targtime
                    no_reads = no_reads + 1
                end do
                if (targtime == time) then
                    print *, "Runtime channel: restarting from time", time
                else
                    print *, "ERROR / RUNTIME CHANNEL, in runtime.f90 - Time of restart file not found. Exiting."
                    stop
                end if
            end if
        end if
        ! broadcast no_reads
        call MPI_bcast(no_reads,MPI_OFFSET_KIND,MPI_BYTE,0,MPI_COMM_WORLD) ! iproc == 0 has terminal
        ! get initial displacement for other files
        disp = no_reads*nzd*8_MPI_OFFSET_KIND
        ! open files for ul and us
        if (has_yl) call MPI_file_open(MPI_COMM_L, TRIM("ul.out"), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh_l)
        if (has_ys) call MPI_file_open(MPI_COMM_S, TRIM("us.out"), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh_s)

        ! write out stuff
        if (has_yl .AND. has_average) write(*,"(A,F11.6)") "Found actual yl:", y(iyl)
        if (has_ys .AND. has_average) write(*,"(A,F11.6)") "Found actual ys:", y(iys)
        if (has_terminal) then
            print *
            write(*,"(A,F11.6)") "   requested cutoff spanwise wavelength:", lc
            write(*,"(A,I5)") "   cutoff spanwise index:", llz
            print *
            open(unit=676767,file='runtime.nfo')
                write(676767,*) "Effective yl:", y(iyl)
                write(676767,*) "Effective ys:", y(iys)
                write(676767,*) "nzd", nzd
                write(676767,*) "Effective dz", 2*PI/beta0 / nzd
            close(676767)
        end if

        call MPI_Barrier(MPI_COMM_WORLD)        

    end subroutine



    subroutine runtime_finalise()
    ! this subroutine is taken from a file in folder "runtime"
    ! the filename is contained in the macro RUNTIME_FINALISE_SUBROUTINE
    ! such macro is set in header.h
#include RUNTIME_FINALISE_SUBROUTINE
        if (has_terminal) then
            endfile(ft) ! truncate file
            close(ft)
        end if
        if (has_yl) then
            call MPI_FILE_SET_SIZE(fh_l, disp) ! truncate file
            call MPI_File_close(fh_l)
        end if
        if (has_ys) then
            call MPI_FILE_SET_SIZE(fh_s, disp) ! truncate file
            call MPI_File_close(fh_s)
        end if
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
            call MPI_file_set_view(fh_l, disp, MPI_DOUBLE_PRECISION, yconst_plane_type, 'native', MPI_INFO_NULL)
            call MPI_File_write_all(fh_l, rVVdx(1,:,1,1), nzB, MPI_DOUBLE_PRECISION, status)

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
            rVVdx(1,:,2,1) = rVVdx(1,:,2,1)*rVVdx(1,:,2,1)
            call MPI_file_set_view(fh_s, disp, MPI_DOUBLE_PRECISION, yconst_plane_type, 'native', MPI_INFO_NULL)
            call MPI_File_write_all(fh_s, rVVdx(1,:,2,1), nzB, MPI_DOUBLE_PRECISION, status)

        end if

        if (has_terminal) then
            write(ft) time
        end if

        ! synchronise
        call MPI_Barrier(MPI_COMM_WORLD)

        ! update displacement
        disp = disp + nzd*8_MPI_OFFSET_KIND

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