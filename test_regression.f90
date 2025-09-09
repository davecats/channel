PROGRAM regression_test
  USE dnsdata
  USE driver
  IMPLICIT NONE

  CHARACTER(len=256) :: config_file, restart_in, restart_expected
  COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE :: V_expected(:, :, :, :)
  REAL(C_DOUBLE) :: diffnorm, tol

  ! Filenames (adapt for your tests/data/ structure)
  config_file = "tests/data/dns_test.in"
  restart_in = "tests/data/start_field.out"
  restart_expected = "tests/data/end_field.out"

  ! Initialise with test input/restart
  CALL initialize(config_file, restart_in)

  CALL timeloop()

  ! Allocate and read expected result
  ALLOCATE (V_expected(SIZE(V, 1), SIZE(V, 2), SIZE(V, 3), SIZE(V, 4)))
  CALL read_restart_file(restart_expected, V_expected)

  ! Compare fields (L2 norm of difference)
  diffnorm = SQRT(SUM(ABS(V - V_expected)**2)/REAL(SIZE(V), C_DOUBLE))
  tol = 1.0e-12_C_DOUBLE

  IF (diffnorm < tol) THEN
    WRITE (*, *) "Regression test PASSED, diffnorm = ", diffnorm
  ELSE
    WRITE (*, *) "Regression test FAILED, diffnorm = ", diffnorm
    CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  END IF

  CALL free_fft()
  CALL free_memory(.TRUE.)
  CALL MPI_Finalize()
END PROGRAM regression_test
