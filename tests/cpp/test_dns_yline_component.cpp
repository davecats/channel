#include "channel/runtime.hpp"
#include "channel/yline.hpp"

#include "test_common.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

namespace {

int penta_index(int row, int line, int batch) {
  return line + row * batch;
}

channel::Complex exact_value(int global_row, int line) {
  return channel::Complex(0.25 + 0.07 * static_cast<double>(global_row) + 0.05 * static_cast<double>(line),
                          0.4 - 0.02 * static_cast<double>(global_row + 2 * line));
}

void coefficients(int global_row,
                  int global_n,
                  int line,
                  channel::Complex& ds,
                  channel::Complex& dl,
                  channel::Complex& d,
                  channel::Complex& du,
                  channel::Complex& dw) {
  ds = global_row >= 2 ? channel::Complex(-0.014 - 0.001 * line, 0.002) : channel::Complex(0.0, 0.0);
  dl = global_row >= 1 ? channel::Complex(-0.08, -0.004 - 0.001 * line) : channel::Complex(0.0, 0.0);
  d = channel::Complex(2.5 + 0.02 * static_cast<double>(global_row), 0.02 + 0.001 * line);
  du = global_row + 1 < global_n ? channel::Complex(-0.075 + 0.001 * line, 0.003) : channel::Complex(0.0, 0.0);
  dw = global_row + 2 < global_n ? channel::Complex(-0.011, -0.002 + 0.001 * line) : channel::Complex(0.0, 0.0);
}

channel::Complex apply_row(int global_row, int global_n, int line) {
  channel::Complex ds;
  channel::Complex dl;
  channel::Complex d;
  channel::Complex du;
  channel::Complex dw;
  coefficients(global_row, global_n, line, ds, dl, d, du, dw);

  channel::Complex rhs = d * exact_value(global_row, line);
  if (global_row >= 1) rhs += dl * exact_value(global_row - 1, line);
  if (global_row >= 2) rhs += ds * exact_value(global_row - 2, line);
  if (global_row + 1 < global_n) rhs += du * exact_value(global_row + 1, line);
  if (global_row + 2 < global_n) rhs += dw * exact_value(global_row + 2, line);
  return rhs;
}

void run_case(const channel::Runtime& runtime,
              channel::ExchangeMode exchange_mode,
              const std::string& label,
              int base_local_n,
              bool uneven_local_n,
              const std::vector<int>& pass_counts) {
  const int local_n = uneven_local_n ? base_local_n + runtime.rank() % 2 : base_local_n;
  constexpr int nx = 3;
  constexpr int nz = 2;
  constexpr int batch = nx * nz;
  constexpr int solve_component = 1;
  std::vector<int> local_counts(static_cast<std::size_t>(runtime.size()));
  const int err = MPI_Allgather(&local_n, 1, MPI_INT, local_counts.data(), 1, MPI_INT, channel::world_comm());
  channel::test::require(err == MPI_SUCCESS, label + " local_n allgather succeeds");

  int global_n = 0;
  int first_global_row = 0;
  for (int rank = 0; rank < runtime.size(); ++rank) {
    if (rank < runtime.rank()) first_global_row += local_counts[static_cast<std::size_t>(rank)];
    global_n += local_counts[static_cast<std::size_t>(rank)];
  }

  channel::DnsState state;
  state.resize({nx, local_n, nz, 3});
  const auto count = state.values_per_component();

  std::vector<channel::Complex> ds(count);
  std::vector<channel::Complex> dl(count);
  std::vector<channel::Complex> d(count);
  std::vector<channel::Complex> du(count);
  std::vector<channel::Complex> dw(count);
  std::vector<channel::Complex> rhs(count);
  std::vector<channel::Complex> sentinel_low(count, channel::Complex(-4.0, 1.0));
  std::vector<channel::Complex> sentinel_high(count, channel::Complex(8.0, -2.0));

  for (int row = 0; row < local_n; ++row) {
    const int global_row = first_global_row + row;
    for (int line = 0; line < batch; ++line) {
      const int p = penta_index(row, line, batch);
      coefficients(global_row, global_n, line, ds[p], dl[p], d[p], du[p], dw[p]);
      rhs[p] = apply_row(global_row, global_n, line);
    }
  }
  state.copy_component_from_host(0, sentinel_low);
  state.copy_component_from_host(solve_component, rhs);
  state.copy_component_from_host(2, sentinel_high);

  channel::DnsYLineComponentSolver solver;
  channel::DnsYLineComponentConfig cfg;
  cfg.component = solve_component;
  cfg.npy = runtime.size();
  cfg.ipy = runtime.rank();
  cfg.pass_counts = pass_counts;
  cfg.exchange_mode = exchange_mode;
  cfg.comm_y = channel::world_comm();
  solver.prepare(state, cfg);
  solver.copy_coefficients_from_host(ds, dl, d, du, dw);
  solver.solve(state);

  const auto solved = state.component_host(solve_component);
  for (int row = 0; row < local_n; ++row) {
    const int global_row = first_global_row + row;
    for (int line = 0; line < batch; ++line) {
      const int p = penta_index(row, line, batch);
      channel::test::require_near(solved[static_cast<std::size_t>(p)], exact_value(global_row, line), 7.0e-11,
                                  label + " DNS y-line component solution rank " +
                                      std::to_string(runtime.rank()) + " row " + std::to_string(row) +
                                      " line " + std::to_string(line));
    }
  }

  const auto low = state.component_host(0);
  const auto high = state.component_host(2);
  for (std::size_t i = 0; i < count; ++i) {
    channel::test::require_near(low[i], sentinel_low[i], 0.0, label + " low component untouched");
    channel::test::require_near(high[i], sentinel_high[i], 0.0, label + " high component untouched");
  }
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1 || runtime.size() == 2 || runtime.size() == 4 || runtime.size() == 6,
                         "DNS y-line component test must run with one, two, four, or six MPI ranks");
  const bool force_pass23 = argc > 1 && std::string(argv[1]) == "pass23";
  channel::test::require(!force_pass23 || runtime.size() == 6,
                         "DNS y-line component pass23 variant requires six MPI ranks");

  std::vector<int> pass_counts;
  if (runtime.size() == 2) {
    pass_counts = {2};
  } else if (runtime.size() == 4) {
    pass_counts = {2, 2};
  } else if (force_pass23) {
    pass_counts = {2, 3};
  }

  run_case(runtime, channel::ExchangeMode::AllToAll, "dns_alltoallv_uniform", 7, false, pass_counts);
  run_case(runtime, channel::ExchangeMode::Auto, "dns_auto_uniform", 7, false, pass_counts);
  run_case(runtime, channel::ExchangeMode::AllToAll, "dns_alltoallv_uneven_local_n", 6, true, pass_counts);
  run_case(runtime, channel::ExchangeMode::Auto, "dns_auto_uneven_local_n", 6, true, pass_counts);

  if (runtime.rank() == 0) std::cout << "DNS y-line component test PASSED\n";
  return 0;
}
