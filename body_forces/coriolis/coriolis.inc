! From the Bachelor's thesis of Mohit Karnvir Singh
! This bodyforce is the Coriolis force associated to a rotation omega in the spanwise direction (z, or 3)
! 2*omega is hardcoded in omega2

subroutine config_body_force()
    if (has_terminal) print *, "Using bodyforce, Coriolis force"
    if (has_terminal) print *, ""
    call set_body_force()
end subroutine

subroutine set_body_force()
integer :: c, r, g, t, pc, zeichen
real, parameter :: omega2 = -0.06
    do c=1,2
        pc=mod(c,2)+1 ! c=1/2 --> pc = 2/1
        zeichen= 2*pc-3
        do r=nx0,nxN
            do g=-nz,nz
                do t=ny0-2,nyn+2
                    F(t,g,r,pc) = zeichen * omega2 * V(t,g,r,c)             
                end do
            end do
        end do
    end do
end subroutine