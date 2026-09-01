! The direct convection-velocity estimator against its exact answer.
!
! For a single travelling mode u(t) = A exp(-i omega t) with A constant in time,
! the three-level operator returns, exactly,
!
!   conj(u_n) du/dt|est = |A|^2 * G,   G = a e^(i omega h1) + b + c e^(-i omega h2)
!
! independently of the phase of A.  Only Im is accumulated, so the stored field
! must equal |A|^2 Im(G) to round-off.  h1 and h2 are deliberately unequal: the
! operator is second order either way, and this is what proves it.
program test_convvelo_direct
  use, intrinsic :: iso_c_binding
  use case_setup, only: nPhi, ny0, nyN, nz, nx0, nxN, V, sync_velocity_to_device, free_memory
  use convvelo, only: init_convvelo, reset_convvelo_stats, free_convvelo, &
                      convvelo_direct_enabled, convvelo_direct_output, &
                      convvelo_prev1, convvelo_prev2, &
                      save_convvelo_time_level, accumulate_direct_numerator, &
                      copy_convvelo_field_average
  use pressure_output, only: free_pressure_output
  use driver, only: initialize
  use mpi_f08
  implicit none

  character(len=*), parameter :: config_file = "tests/convvelo/dns.in"
  character(len=*), parameter :: restart_in = "tests/convvelo/Dati.cart.out"

  ! Unequal on purpose, and by much more than CFL control ever produces.
  real(C_DOUBLE), parameter :: h1 = 0.011d0, h2 = 0.017d0
  real(C_DOUBLE), parameter :: tol = 1.0d-12

  complex(C_DOUBLE_COMPLEX), allocatable :: field(:, :, :)
  complex(C_DOUBLE_COMPLEX) :: g_response, expected, phase
  real(C_DOUBLE) :: ca, cb, cc, omega, worst, err
  integer(C_INT) :: ic, ix, iy, iz, it, output_index, nfail
  real(C_DOUBLE) :: t_level(3)

  t_level = [0.0d0, h1, h1 + h2]
  ca = -h2/(h1*(h1 + h2))
  cb = (h2 - h1)/(h1*h2)
  cc = h1/(h2*(h1 + h2))

  nfail = 0
  worst = 0.0d0

  call initialize(config_file, restart_in)
  ! configure_convvelo has already run and left this off, since the test deck
  ! carries no [convvelo] section.  Turn it on before the layout is built.
  convvelo_direct_enabled = .true.
  call init_convvelo()
  call reset_convvelo_stats()

  do it = 1, 3
    do ic = 1, 3 + nPhi
      omega = mode_frequency(ic)
      phase = exp(cmplx(0.0d0, -omega*t_level(it), C_DOUBLE_COMPLEX))
      do ix = nx0, nxN
        do iz = -nz, nz
          do iy = ny0 - 2, nyN + 2
            V(iy, iz, ix, ic) = amplitude(iy, iz, ix, ic)*phase
          end do
        end do
      end do
    end do
    call sync_velocity_to_device()
    if (it == 1) call save_convvelo_time_level(convvelo_prev1)
    if (it == 2) call save_convvelo_time_level(convvelo_prev2)
    if (it == 3) call accumulate_direct_numerator(h1, h2)
  end do

  allocate (field(ny0 - 2:nyN + 2, -nz:nz, nx0:nxN))

  do ic = 1, 3 + nPhi
    output_index = convvelo_direct_output(ic)
    if (output_index == 0) then
      write (*, *) "FAIL: no dt output field for component ", ic
      nfail = nfail + 1
      cycle
    end if
    omega = mode_frequency(ic)
    g_response = ca*exp(cmplx(0.0d0, omega*h1, C_DOUBLE_COMPLEX)) + cb + &
                 cc*exp(cmplx(0.0d0, -omega*h2, C_DOUBLE_COMPLEX))
    call copy_convvelo_field_average(output_index, field)
    do ix = nx0, nxN
      do iz = -nz, nz
        do iy = ny0 - 2, nyN + 2
          expected = cmplx(0.0d0, abs(amplitude(iy, iz, ix, ic))**2*aimag(g_response), C_DOUBLE_COMPLEX)
          err = abs(field(iy, iz, ix) - expected)/max(abs(expected), 1.0d0)
          worst = max(worst, err)
          if (err > tol) nfail = nfail + 1
        end do
      end do
    end do
  end do

  write (*, '(A,ES12.5)') " worst relative error = ", worst
  if (nfail == 0) then
    write (*, *) "Convvelo direct estimator PASSED"
  else
    write (*, *) "Convvelo direct estimator FAILED, mismatching points = ", nfail
  end if

  deallocate (field)
  call free_pressure_output()
  call free_convvelo()
  call free_memory(.true.)
  call MPI_Finalize()
  if (nfail /= 0) error stop 1

contains

  ! Distinct per component so a swapped component shows up as a failure.
  real(C_DOUBLE) function mode_frequency(component)
    integer(C_INT), intent(in) :: component
    mode_frequency = 3.0d0 + 1.75d0*real(component, C_DOUBLE)
  end function mode_frequency

  ! Complex, bounded, nowhere zero, and varying in every index.
  complex(C_DOUBLE_COMPLEX) function amplitude(iy, iz, ix, component)
    integer(C_INT), intent(in) :: iy, iz, ix, component
    amplitude = cmplx(1.0d0 + 0.125d0*real(component, C_DOUBLE) + 0.0625d0*real(iy, C_DOUBLE), &
                      0.5d0 + 0.03125d0*real(iz, C_DOUBLE) - 0.015625d0*real(ix, C_DOUBLE), &
                      C_DOUBLE_COMPLEX)
  end function amplitude

end program test_convvelo_direct
