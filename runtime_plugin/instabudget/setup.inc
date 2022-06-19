        allocate(dVdy(ny0-2:nyN+2,-nz:nz,nx0:nxN)); dVdy=0

        open(15, file='instabudget.in')
            read(15, *) lc
            read(15, *) ys, yl
            read(15, *) remove_spatial_average
            read(15, *) runtime_read_field_from_idx, runtime_tgt_file_idx
        close(15)
        if (runtime_read_field_from_idx) then
            write(rt_tgt_istr,*) runtime_tgt_file_idx
        else
            rt_tgt_istr = ""
        end if

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
            open(unit=ft,file="time."//TRIM(ADJUSTL(rt_tgt_istr))//".out",access="stream",action="readwrite",position="rewind")
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
        if (has_yl) call MPI_file_open(MPI_COMM_L, TRIM("ul."//TRIM(ADJUSTL(rt_tgt_istr))//".out"), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh_l)
        if (has_ys) then
            call MPI_file_open(MPI_COMM_S, TRIM("us."//TRIM(ADJUSTL(rt_tgt_istr))//".out"), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh_us)
            call MPI_file_open(MPI_COMM_S, TRIM("vs."//TRIM(ADJUSTL(rt_tgt_istr))//".out"), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh_vs)
            call MPI_file_open(MPI_COMM_S, TRIM("ws."//TRIM(ADJUSTL(rt_tgt_istr))//".out"), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh_ws)
            call MPI_file_open(MPI_COMM_S, TRIM("uvs."//TRIM(ADJUSTL(rt_tgt_istr))//".out"), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh_uvs)
            call MPI_file_open(MPI_COMM_S, TRIM("duldy."//TRIM(ADJUSTL(rt_tgt_istr))//".out"), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh_duldy)
            call MPI_file_open(MPI_COMM_S, TRIM("ul_at_ys."//TRIM(ADJUSTL(rt_tgt_istr))//".out"), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh_ul_at_ys)
        end if

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