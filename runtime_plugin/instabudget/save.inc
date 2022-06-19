    integer :: comp

        ! get velocity gradient
#ifdef nonblockingY
        cnt_req = 1
        call MPI_Buffer_Attach(BUFveta, szBUFveta)
#endif
        do ix = nx0, nxN
            do iz = -nz, nz
#ifdef nonblockingY
                call COMPLEXderiv(V(:,iz,ix,1), dVdy(:,iz,ix),REQvetaTOuvw(cnt_req),iz+nz); cnt_req=cnt_req+1
#else
                call COMPLEXderiv(V(:,iz,ix,1), dVdy(:,iz,ix))
                call LeftLU5divStep2(D0mat, dVdy(:,iz,ix))
#endif
            end do
        end do
#ifdef nonblockingY
        CALL MPI_Buffer_Detach(ptrBUFveta,szBUFveta)
        cnt_req = 1
        call MPI_Buffer_Attach(BUFveta, szBUFveta)
        do ix = nx0, nxN
            do iz = -nz, nz
                call LeftLU5divStep2(D0mat,dVdy(:,iz,ix),REQvetaTOuvw(cnt_req),iz+nz); cnt_req=cnt_req+1
            end do
        end do
        CALL MPI_Buffer_Detach(ptrBUFveta,szBUFveta)
#endif

        ! large scales
        if (has_yl) then

            comp = 1

            VVdz(:,:,1,1) = 0 ! clean up VVdz first
            ! copy velocity to FFTW array
            call copy2vvdz_large(1,1,V(iyl,:,:,comp))

            ! remove average
            if (remove_spatial_average .AND. has_average) then
                VVdz(1,1,1,1) = 0
            end if

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

            VVdz(:,:,2,:) = 0 ! clean up VVdz first
            
            ! copy velocity
            do comp = 1,3 ! three components
                call copy2vvdz_small(2,comp,V(iys,:,:,comp))
            end do

            ! copy the velocity gradient
            ! ATTENTION: THIS IS LARGE SCALE GRADIENT AT YS
            call copy2vvdz_large(2,4,dVdy(iys,:,:))

            ! copy large scales at ys
            call copy2vvdz_large(2,5,V(iys,:,:,1))
    
            ! remove average
            if (remove_spatial_average .AND. has_average) then
                VVdz(1,1,2,1:5) = 0 ! remove mode 0,0 from all vel comp and from vel grad!
            end if

            ! do the transform
            do comp = 1,5
                call IFT(VVdz(1:nzd,1:nxB,2,comp))  
#ifdef nonblockingXZ
                call MPI_IAlltoall(VVdz(:,:,2,comp), 1, Mdz, VVdx(:,:,2,comp), 1, Mdx, MPI_COMM_X, req_rt(2))
#endif
#ifndef nonblockingXZ
                call MPI_Alltoall(VVdz(:,:,2,comp), 1, Mdz, VVdx(:,:,2,comp), 1, Mdx, MPI_COMM_X)
                VVdx(nx+2:nxd+1,1:nzB,2,comp)=0;    CALL RFT(VVdx(1:nxd+1,1:nzB,2,comp),rVVdx(1:2*nxd+2,1:nzB,2,comp))
#endif
#ifdef nonblockingXZ
                call MPI_wait(req_rt(2),MPI_STATUS_IGNORE); 
                VVdx(nx+2:nxd+1,1:nzB,2,comp)=0; 
                call RFT(VVdx(1:nxd+1,1:nzB,2,comp),rVVdx(1:2*nxd+2,1:nzB,2,comp));
#endif
            end do

            ! calculate stuff to save to disk
            rVVdx(1,:,2,6) = rVVdx(1,:,2,1)*rVVdx(1,:,2,2) ! us*vs --- MUST BE FIRST
            do comp = 1,3
                rVVdx(1,:,2,comp) = rVVdx(1,:,2,comp)*rVVdx(1,:,2,comp) ! uis*uis
            end do
            ! velocity gradient rVVdx(1,:,2,4) is already good as it is

            ! set view
            call MPI_file_set_view(fh_us,    disp, MPI_DOUBLE_PRECISION, yconst_plane_type, 'native', MPI_INFO_NULL)
            call MPI_file_set_view(fh_vs,    disp, MPI_DOUBLE_PRECISION, yconst_plane_type, 'native', MPI_INFO_NULL)
            call MPI_file_set_view(fh_ws,    disp, MPI_DOUBLE_PRECISION, yconst_plane_type, 'native', MPI_INFO_NULL)
            call MPI_file_set_view(fh_uvs,   disp, MPI_DOUBLE_PRECISION, yconst_plane_type, 'native', MPI_INFO_NULL)
            call MPI_file_set_view(fh_duldy, disp, MPI_DOUBLE_PRECISION, yconst_plane_type, 'native', MPI_INFO_NULL)
            call MPI_file_set_view(fh_ul_at_ys, disp, MPI_DOUBLE_PRECISION, yconst_plane_type, 'native', MPI_INFO_NULL)
            ! write
            call MPI_File_write_all(fh_us,    rVVdx(1,:,2,1), nzB, MPI_DOUBLE_PRECISION, status)
            call MPI_File_write_all(fh_vs,    rVVdx(1,:,2,2), nzB, MPI_DOUBLE_PRECISION, status)
            call MPI_File_write_all(fh_ws,    rVVdx(1,:,2,3), nzB, MPI_DOUBLE_PRECISION, status)
            call MPI_File_write_all(fh_uvs,   rVVdx(1,:,2,6), nzB, MPI_DOUBLE_PRECISION, status)
            call MPI_File_write_all(fh_duldy, rVVdx(1,:,2,4), nzB, MPI_DOUBLE_PRECISION, status)
            call MPI_File_write_all(fh_ul_at_ys, rVVdx(1,:,2,5), nzB, MPI_DOUBLE_PRECISION, status)
        end if

        if (has_terminal) then
            write(ft) time
        end if

        ! synchronise
        call MPI_Barrier(MPI_COMM_WORLD)

        ! update displacement
        disp = disp + nzd*8_MPI_OFFSET_KIND