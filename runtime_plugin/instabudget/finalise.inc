        if (has_terminal) then
            endfile(ft) ! truncate file
            close(ft)
        end if
        if (has_yl) then
            call MPI_FILE_SET_SIZE(fh_l, disp) ! truncate file
            call MPI_File_close(fh_l)
        end if
        if (has_ys) then
            ! truncate file
            call MPI_FILE_SET_SIZE(fh_us, disp) 
            call MPI_FILE_SET_SIZE(fh_vs, disp) 
            call MPI_FILE_SET_SIZE(fh_ws, disp) 
            call MPI_FILE_SET_SIZE(fh_uvs, disp) 
            call MPI_FILE_SET_SIZE(fh_duldy, disp)
            call MPI_FILE_SET_SIZE(fh_ul_at_ys, disp)
            ! close file handle
            call MPI_File_close(fh_us)
            call MPI_File_close(fh_vs)
            call MPI_File_close(fh_ws)
            call MPI_File_close(fh_uvs)
            call MPI_File_close(fh_duldy)
            call MPI_File_close(fh_ul_at_ys)
        end if