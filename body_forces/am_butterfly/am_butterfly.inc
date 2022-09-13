! Suppresses superposition of large scales at the wall
! Parameters are hardcoded for Poiseuille MSU, retau 1000

subroutine config_body_force()
    if (has_terminal) print *, "Using bodyforce, am_buttefly (suppression of superposition + small scales)"
    if (has_terminal) print *, ""
    iz_f = NINT((2.d0*PI/lambdaz_f)/(beta0/1000.d0))
    call set_body_force()
end subroutine

subroutine set_body_force()
integer(C_INT) :: ix,iz,iy,iV
real(C_DOUBLE) :: yp
    do iV=1,3
        do ix=nx0,nxN
            do iz=(-nz),nz
                do iy=ny0-2,nyn+2
                    yp = MERGE(ymax-y(iy), y(iy), y(iy)>1)*1000
                    ! lambdaz+ = 2.3 (y+)^2
                    if ( (abs(iz) <= iz_f) .AND. (yp <= 60) .AND. (.NOT. (ix==0 .AND. iz==0)) ) then 
                        F(iy,iz,ix,iV)=-amp*V(iy,iz,ix,iV)
                    else if ( (abs(iz) > iz_f) .AND. (yp > 60) ) then
                        F(iy,iz,ix,iV)=-amp*V(iy,iz,ix,iV)
                    end if
                end do
            end do
        end do
    end do
end subroutine