#include "channel/input.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1, "input parser test runs on one MPI rank");

  const std::string path = "test_channel_input.ini";
  {
    std::ofstream out(path);
    out << "[mesh]\n"
        << "nx = 7\n"
        << "ny = 9\n"
        << "nz = 5\n"
        << "alfa0 = 0.25\n"
        << "beta0 = 1.5\n"
        << "stretching = 1.2\n"
        << "ymin = 0\n"
        << "ymax = 2\n"
        << "[velocity]\n"
        << "ni = 1000\n"
        << "meanpx = 0.1\n"
        << "meanpz = -0.2\n"
        << "meanflowx = 2\n"
        << "meanflowz = 0.5\n"
        << "u0 = 0\n"
        << "uN = 0\n"
        << "[scalars]\n"
        << "nPhi = 2\n"
        << "pr = 2 4\n"
        << "meantx = 0\n"
        << "meantb = 0\n"
        << "t0 = 0\n"
        << "tN = 0\n"
        << "[parallel]\n"
        << "npy = 3\n"
        << "[timestepping]\n"
        << "deltat = 0.01\n"
        << "cflmax = 1\n"
        << "time = 4\n"
        << "dt_field = 1\n"
        << "dt_save = 1\n"
        << "t_max = 5\n"
        << "time_from_restart = false\n"
        << "nstep = 6\n";
  }

  const auto cfg = channel::read_channel_input(path);
  channel::test::require(cfg.mesh.nx == 7 && cfg.mesh.ny == 9 && cfg.mesh.nz == 5, "mesh dimensions parsed");
  channel::test::require(std::abs(cfg.velocity.viscosity - 0.001) < 1.0e-15, "velocity.ni inverted");
  channel::test::require(cfg.scalars.nphi == 2, "scalar count parsed");
  channel::test::require(std::abs(cfg.scalars.inverse_prandtl[0] - 0.5) < 1.0e-15, "Prandtl vector inverted");
  channel::test::require(cfg.parallel.npy == 3 && cfg.parallel.npy_was_set, "parallel npy parsed");
  channel::test::require(cfg.timestepping.nstep == 6 && std::abs(cfg.timestepping.time - 4.0) < 1.0e-15,
                         "timestepping parsed");

  if (runtime.rank() == 0) std::cout << "input parser test PASSED\n";
  return 0;
}
