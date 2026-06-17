#include "channel/runtime.hpp"
#include "channel/schur.hpp"

#include "test_common.hpp"

#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();

  channel::SchurSolver solver;
  channel::SchurSolverConfig cfg;
  cfg.npy = 8;
  cfg.ipy = 5;
  cfg.nlines = 96;
  cfg.pass_counts = {4, 2};
  cfg.exchange_mode = channel::ExchangeMode::Auto;
  cfg.comm_y = channel::world_comm();
  solver.prepare(cfg);

  const auto& levels = solver.levels();
  channel::test::require(levels.size() == 2, "one SchurLevel object per configured pass");
  channel::test::require(levels[0].arity() == 4, "level 0 arity");
  channel::test::require(levels[1].arity() == 2, "level 1 arity");
  channel::test::require(levels[0].prev_count() == 96, "level 0 previous line count");
  channel::test::require(levels[0].owned_count() == 24, "level 0 owned count");
  channel::test::require(levels[1].prev_count() == 24, "level 1 previous line count");
  channel::test::require(levels[1].owned_count() == 12, "level 1 owned count");
  channel::test::require(levels[1].config().exchange_mode == channel::ExchangeMode::AllGather,
                         "auto exchange resolves arity-2 root to allgather");

  auto bad_cfg = cfg;
  bad_cfg.ipy = -1;
  channel::test::require_throws([&] { solver.prepare(bad_cfg); }, "negative y-rank is rejected");

  bad_cfg = cfg;
  bad_cfg.ipy = bad_cfg.npy;
  channel::test::require_throws([&] { solver.prepare(bad_cfg); }, "out-of-range y-rank is rejected");

  bad_cfg = cfg;
  bad_cfg.nlines = 0;
  channel::test::require_throws([&] { solver.prepare(bad_cfg); }, "empty Schur line set is rejected");

  bad_cfg = cfg;
  bad_cfg.pass_counts = {3, 3};
  channel::test::require_throws([&] { solver.prepare(bad_cfg); }, "invalid Schur pass sequence is rejected");

  channel::SchurSolver single_rank_solver;
  channel::SchurSolverConfig single_cfg;
  single_cfg.npy = 1;
  single_cfg.ipy = 0;
  single_cfg.nlines = 3;
  single_cfg.pass_counts = {};
  single_cfg.exchange_mode = channel::ExchangeMode::Auto;
  single_cfg.comm_y = channel::world_comm();
  single_rank_solver.prepare(single_cfg);

  std::vector<channel::Complex> leaf_rows(
      static_cast<std::size_t>(single_cfg.nlines * channel::SchurLevel::row_width), channel::Complex(0.0, 0.0));
  for (int line = 0; line < single_cfg.nlines; ++line) {
    for (int k = 0; k < channel::SchurLevel::row_width; ++k) {
      leaf_rows[static_cast<std::size_t>(line * channel::SchurLevel::row_width + k)] =
          channel::Complex(0.25 * line + 0.01 * k, -0.02 * k);
    }
  }
  channel::DeviceVector<channel::Complex> leaf_rows_device("single_rank_schur_rows", leaf_rows.size());
  leaf_rows_device.copy_from_host(leaf_rows);
  channel::DeviceVector<channel::Complex> leaf_values_device;
  single_rank_solver.solve_from_leaf_rows(leaf_rows_device, leaf_values_device);
  const auto leaf_values = leaf_values_device.copy_to_host();
  for (int line = 0; line < single_cfg.nlines; ++line) {
    for (int k = 0; k < channel::SchurLevel::value_width; ++k) {
      const auto expected = k < 4 ? leaf_rows[static_cast<std::size_t>(line * channel::SchurLevel::row_width + k)]
                                  : channel::Complex(0.0, 0.0);
      channel::test::require_near(
          leaf_values[static_cast<std::size_t>(line * channel::SchurLevel::value_width + k)], expected, 1.0e-12,
          "single-rank Schur leaf value");
    }
  }

  solver.release();
  single_rank_solver.release();
  if (runtime.rank() == 0) std::cout << "SchurLevel object model test PASSED\n";
  return 0;
}
