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

channel::Complex exact_value(int component, int global_row, int line) {
  return channel::Complex(0.2 + 0.11 * component + 0.05 * global_row + 0.03 * line,
                          -0.15 + 0.04 * component - 0.01 * (global_row + line));
}

void coefficients(int component,
                  int global_row,
                  int global_n,
                  int line,
                  channel::Complex& ds,
                  channel::Complex& dl,
                  channel::Complex& d,
                  channel::Complex& du,
                  channel::Complex& dw) {
  const double c = static_cast<double>(component + 1);
  ds = global_row >= 2 ? channel::Complex(-0.010 - 0.001 * c, 0.001 * (line + 1)) : channel::Complex(0.0, 0.0);
  dl = global_row >= 1 ? channel::Complex(-0.070 - 0.002 * c, -0.003) : channel::Complex(0.0, 0.0);
  d = channel::Complex(2.7 + 0.04 * c + 0.01 * global_row, 0.015 + 0.001 * line);
  du = global_row + 1 < global_n ? channel::Complex(-0.060 + 0.001 * line, 0.002 * c)
                                  : channel::Complex(0.0, 0.0);
  dw = global_row + 2 < global_n ? channel::Complex(-0.009, -0.001 * (line + component + 1))
                                  : channel::Complex(0.0, 0.0);
}

channel::Complex apply_row(int component, int global_row, int global_n, int line) {
  channel::Complex ds;
  channel::Complex dl;
  channel::Complex d;
  channel::Complex du;
  channel::Complex dw;
  coefficients(component, global_row, global_n, line, ds, dl, d, du, dw);

  channel::Complex rhs = d * exact_value(component, global_row, line);
  if (global_row >= 1) rhs += dl * exact_value(component, global_row - 1, line);
  if (global_row >= 2) rhs += ds * exact_value(component, global_row - 2, line);
  if (global_row + 1 < global_n) rhs += du * exact_value(component, global_row + 1, line);
  if (global_row + 2 < global_n) rhs += dw * exact_value(component, global_row + 2, line);
  return rhs;
}

void fill_component_system(int component,
                           int local_n,
                           int first_global_row,
                           int global_n,
                           int batch,
                           std::vector<channel::Complex>& ds,
                           std::vector<channel::Complex>& dl,
                           std::vector<channel::Complex>& d,
                           std::vector<channel::Complex>& du,
                           std::vector<channel::Complex>& dw,
                           std::vector<channel::Complex>& rhs) {
  for (int row = 0; row < local_n; ++row) {
    const int global_row = first_global_row + row;
    for (int line = 0; line < batch; ++line) {
      const int p = penta_index(row, line, batch);
      coefficients(component, global_row, global_n, line, ds[p], dl[p], d[p], du[p], dw[p]);
      rhs[p] = apply_row(component, global_row, global_n, line);
    }
  }
}

void run_case(const channel::Runtime& runtime,
              const std::string& label,
              channel::ExchangeMode exchange_mode,
              bool uneven_local_n,
              const std::vector<int>& pass_counts) {
  constexpr int nx = 2;
  constexpr int nz = 3;
  constexpr int batch = nx * nz;
  constexpr int component_count = 3;
  const int local_n = uneven_local_n ? 6 + runtime.rank() % 2 : 7;

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
  state.resize({nx, local_n, nz, component_count});
  const auto count = state.values_per_component();
  std::vector<channel::Complex> sentinel(count, channel::Complex(9.0, -3.0));
  state.copy_component_from_host(1, sentinel);

  channel::DnsImplicitYLineStep step;
  channel::DnsImplicitYLineStepConfig cfg;
  cfg.components = {0, 2};
  cfg.npy = runtime.size();
  cfg.ipy = runtime.rank();
  cfg.pass_counts = pass_counts;
  cfg.exchange_mode = exchange_mode;
  cfg.comm_y = channel::world_comm();
  step.prepare(state, cfg);

  for (const int component : cfg.components) {
    std::vector<channel::Complex> ds(count);
    std::vector<channel::Complex> dl(count);
    std::vector<channel::Complex> d(count);
    std::vector<channel::Complex> du(count);
    std::vector<channel::Complex> dw(count);
    std::vector<channel::Complex> rhs(count);
    fill_component_system(component, local_n, first_global_row, global_n, batch, ds, dl, d, du, dw, rhs);
    state.copy_component_from_host(component, rhs);
    step.copy_coefficients_from_host(component, ds, dl, d, du, dw);
  }

  step.solve(state);

  for (const int component : cfg.components) {
    const auto solved = state.component_host(component);
    for (int row = 0; row < local_n; ++row) {
      const int global_row = first_global_row + row;
      for (int line = 0; line < batch; ++line) {
        const int p = penta_index(row, line, batch);
        channel::test::require_near(solved[static_cast<std::size_t>(p)], exact_value(component, global_row, line),
                                    8.0e-11, label + " implicit DNS component solution");
      }
    }
  }

  const auto untouched = state.component_host(1);
  for (std::size_t i = 0; i < untouched.size(); ++i) {
    channel::test::require_near(untouched[i], sentinel[i], 0.0, label + " inactive component untouched");
  }

  channel::test::require(step.component_count() == 2, label + " prepared component count");
  channel::test::require(step.component_solver(0).config().component == 0, label + " component 0 lookup");
  channel::test::require(step.component_solver(2).config().component == 2, label + " component 2 lookup");
  channel::test::require_throws([&] { (void)step.component_solver(1); }, label + " missing component lookup rejected");
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1 || runtime.size() == 2 || runtime.size() == 4 || runtime.size() == 6,
                         "implicit DNS y-line step test must run with one, two, four, or six MPI ranks");
  const bool force_pass23 = argc > 1 && std::string(argv[1]) == "pass23";
  channel::test::require(!force_pass23 || runtime.size() == 6,
                         "implicit DNS y-line step pass23 variant requires six MPI ranks");

  std::vector<int> pass_counts;
  if (runtime.size() == 2) {
    pass_counts = {2};
  } else if (runtime.size() == 4) {
    pass_counts = {2, 2};
  } else if (force_pass23) {
    pass_counts = {2, 3};
  }

  channel::DnsState validation_state;
  validation_state.resize({2, 7, 2, 3});
  channel::DnsImplicitYLineStep validation_step;
  channel::DnsImplicitYLineStepConfig bad_cfg;
  bad_cfg.components = {0, 0};
  channel::test::require_throws([&] { validation_step.prepare(validation_state, bad_cfg); },
                                "duplicate implicit DNS components rejected");
  bad_cfg.components = {};
  channel::test::require_throws([&] { validation_step.prepare(validation_state, bad_cfg); },
                                "empty implicit DNS component list rejected");

  run_case(runtime, "implicit_dns_alltoallv_uniform", channel::ExchangeMode::AllToAll, false, pass_counts);
  run_case(runtime, "implicit_dns_auto_uniform", channel::ExchangeMode::Auto, false, pass_counts);
  run_case(runtime, "implicit_dns_alltoallv_uneven", channel::ExchangeMode::AllToAll, true, pass_counts);
  run_case(runtime, "implicit_dns_auto_uneven", channel::ExchangeMode::Auto, true, pass_counts);

  if (runtime.rank() == 0) std::cout << "DNS implicit y-line step test PASSED\n";
  return 0;
}
