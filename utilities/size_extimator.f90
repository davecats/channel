program size_extimator
use iso_c_binding, only: C_INT
implicit none

integer(C_INT) :: nx, ny, nz

    ! read from dns.in
    open(15, file='dns.in')
            read(15, *) nx, ny, nz;
    close(15)

    print *, "Extimated field size: ", 3*16*((nx+1.0)/1024)*((ny+3.0)/1024)*((2*nz+1.0)/1024), " GB"

end program