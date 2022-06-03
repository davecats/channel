! From the Bachelor's thesis of Mohit Karnvir Singh
! This bodyforce is the Coriolis force associated to a rotation omega in the spanwise direction (z, or 3)

subroutine config_body_force()
    if (has_terminal) print *, "Using bodyforce, Coriolis force"
    if (has_terminal) print *, "Reading Ro from coriolis.in"
    ! read coriolis.in
    open(15, file='camstar.in')
        read(15, *) omega2
        read(15, *) kz_cutoff
        read(15, *) y_threshold_bot
    close(15)
    ! process info from coriolis.in
    y_threshold_top = ymax - y_threshold_bot
    iz_thr = min(nz,floor(kz_cutoff/beta0))
    ! print out parameters
    if (has_terminal) print *, "Ro(tation number)", omega2/2
    if (has_terminal) print *, "kz_cutoff, iz_cutoff", kz_cutoff, iz_thr
    if (has_terminal) print *, "y_threshold_bottom", y_threshold_bot
    if (has_terminal) print *, ""
    ! sync and set body force
    call MPI_Barrier(MPI_COMM_WORLD)
    call set_body_force()
end subroutine

subroutine set_body_force()
integer :: c, ix, iz, iy, pc, zeichen
    do c=1,2
        pc=mod(c,2)+1 ! c=1/2 --> pc = 2/1
        zeichen= 2*pc-3
        do ix=nx0,nxN
            do iz=-iz_thr,iz_thr
                do iy=ny0-2,nyn+2
                    if (y(iy) <= y_threshold_bot .OR. y(iy) >= y_threshold_top) then
                        F(iy,iz,ix,pc) = zeichen * omega2 * V(iy,iz,ix,c)
                    end if           
                end do
            end do
        end do
    end do
end subroutine