real(C_DOUBLE) :: omega2 ! minus one times rotation number (that's right: it's -Ro)
real(C_DOUBLE) :: kz_cutoff ! for large/small scale discrimination
real(C_DOUBLE) :: y_threshold_top, y_threshold_bot ! for confining Coriolis force to a wall-normal region
integer :: iz_thr