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



    subroutine copy2vvdz_large(ivvdz1, ivvdz2, R)
    complex(C_DOUBLE_COMPLEX), intent(IN) :: R(-nz:nz,nx0:nxN)
    integer(C_INT), intent(IN) :: ivvdz1,ivvdz2
        VVdz(1:(llz+1),1:nxB,ivvdz1,ivvdz2)=R(0:llz,nx0:nxN);
        VVdz((nzd+1-llz):nzd,1:nxB,ivvdz1,ivvdz2)=R((-llz):(-1),nx0:nxN);
    end subroutine



    subroutine copy2vvdz_small(ivvdz1,ivvdz2, R)
    complex(C_DOUBLE_COMPLEX), intent(IN) :: R(-nz:nz,nx0:nxN)
    integer(C_INT), intent(IN) :: ivvdz1,ivvdz2
        VVdz((llz+2):(nz+1),1:nxB,ivvdz1,ivvdz2)=R((llz+1):nz,nx0:nxN);
        VVdz((nzd+1-nz):(nzd-llz),1:nxB,ivvdz1,ivvdz2)=R((-nz):(-llz-1),nx0:nxN);
    end subroutine