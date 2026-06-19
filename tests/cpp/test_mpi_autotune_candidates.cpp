#include "channel/runtime.hpp"
#include "channel/types.hpp"

#include "test_common.hpp"

#include <iostream>

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();

  channel::test::require(channel::default_schur_pass_counts(1).empty(), "npy=1 has no Schur passes");
  channel::test::require((channel::default_schur_pass_counts(8) == std::vector<int>{4, 2}),
                         "npy=8 default factors prefer 4 then 2");
  channel::test::require(channel::valid_schur_pass_sequence(12, {3, 4}), "3x4 sequence valid for npy=12");
  channel::test::require(!channel::valid_schur_pass_sequence(12, {5, 2}), "unsupported arity rejected");
  channel::test::require(!channel::valid_schur_pass_sequence(12, {3, 3}), "wrong product rejected");

  if (runtime.rank() == 0) std::cout << "MPI autotune candidate helper test PASSED\n";
  return 0;
}
