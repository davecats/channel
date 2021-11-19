subroutine set_coriolis_READ_WALZ()
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