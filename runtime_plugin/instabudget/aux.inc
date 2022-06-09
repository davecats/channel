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