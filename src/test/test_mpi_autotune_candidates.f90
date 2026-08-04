#include "build_options.h"

program test_mpi_autotune_candidates
  use, intrinsic :: iso_c_binding
  use mpi_autotune, only: mpi_autotune_has_pass_sequence
  implicit none

  integer :: nfail
  integer(C_INT), allocatable :: empty(:)

  nfail = 0
  allocate (empty(0))

  call expect(.not. mpi_autotune_has_pass_sequence(1_C_INT, [2_C_INT]), "npy=1 rejects nonempty pass list")
  call expect(mpi_autotune_has_pass_sequence(1_C_INT, empty), "npy=1 accepts empty pass list")

  call expect(mpi_autotune_has_pass_sequence(32_C_INT, [4_C_INT, 4_C_INT, 2_C_INT]), "32 includes 4,4,2")
  call expect(mpi_autotune_has_pass_sequence(32_C_INT, [4_C_INT, 2_C_INT, 2_C_INT, 2_C_INT]), "32 includes 4,2,2,2")
  call expect(mpi_autotune_has_pass_sequence(32_C_INT, [2_C_INT, 2_C_INT, 2_C_INT, 2_C_INT, 2_C_INT]), "32 includes 2,2,2,2,2")
  call expect(mpi_autotune_has_pass_sequence(32_C_INT, [8_C_INT, 4_C_INT]), "32 includes arity-8 option")

  call expect(mpi_autotune_has_pass_sequence(24_C_INT, [4_C_INT, 3_C_INT, 2_C_INT]), "24 includes 4,3,2")
  call expect(mpi_autotune_has_pass_sequence(24_C_INT, [4_C_INT, 6_C_INT]), "24 includes 4,6")
  call expect(mpi_autotune_has_pass_sequence(24_C_INT, [2_C_INT, 2_C_INT, 3_C_INT, 2_C_INT]), "24 includes 2,2,3,2")
  call expect(mpi_autotune_has_pass_sequence(24_C_INT, [8_C_INT, 3_C_INT]), "24 includes arity-8 option")

  call expect(.not. mpi_autotune_has_pass_sequence(24_C_INT, [5_C_INT, 4_C_INT]), "rejects unsupported arity")
  call expect(.not. mpi_autotune_has_pass_sequence(24_C_INT, [4_C_INT, 4_C_INT]), "rejects wrong product")

  if (nfail /= 0) error stop "MPI autotune candidate tests failed"

contains

  subroutine expect(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      nfail = nfail + 1
      write (*, *) "FAILED: ", trim(label)
    end if
  end subroutine expect

end program test_mpi_autotune_candidates
