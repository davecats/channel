#include "channel/decomposition.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <iostream>

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::DecompositionConfig cfg;
  cfg.nxpp = 8;
  cfg.nz = 1;
  cfg.ny = 16;
  cfg.nzd = 3;
  cfg.nphi = 1;
  cfg.npy = 1;

  channel::Decomposition decomp;
  decomp.initialize(cfg);
  channel::test::require(decomp.nproc == runtime.size(), "nproc matches runtime");
  channel::test::require(decomp.iproc == runtime.rank(), "iproc matches runtime");
  channel::test::require(decomp.npy == 1, "npy preserved");
  channel::test::require(decomp.npxz == runtime.size(), "npxz equals rank count for npy=1");
  channel::test::require(decomp.ny0 == 1, "ny0 for npy=1");
  channel::test::require(decomp.nyN == 15, "nyN for npy=1");
  channel::test::require(decomp.sendbuf.size() >= 1, "send buffer allocated");
  decomp.finalize();

  if (runtime.rank() == 0) std::cout << "Decomposition test PASSED\n";
  return 0;
}
