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
  return channel::Complex(0.7 + 0.13 * static_cast<double>(global_row) + 0.04 * static_cast<double>(line),
                          -0.2 + 0.03 * static_cast<double>(global_row - line));
}

void coefficients(int global_row,
                  int global_n,
                  int line,
                  channel::Complex& ds,
                  channel::Complex& dl,
                  channel::Complex& d,
                  channel::Complex& du,
                  channel::Complex& dw) {
  ds = global_row >= 2 ? channel::Complex(-0.018 - 0.001 * line, 0.004) : channel::Complex(0.0, 0.0);
  dl = global_row >= 1 ? channel::Complex(-0.11, -0.006 - 0.001 * line) : channel::Complex(0.0, 0.0);
  d = channel::Complex(2.8 + 0.015 * static_cast<double>(global_row), 0.03 + 0.002 * line);
  du = global_row + 1 < global_n ? channel::Complex(-0.09 + 0.002 * line, 0.005) : channel::Complex(0.0, 0.0);
  dw = global_row + 2 < global_n ? channel::Complex(-0.013, -0.003 + 0.001 * line) : channel::Complex(0.0, 0.0);
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
  const int local_n = uneven_local_n ? base_local_n + runtime.rank() % 3 : base_local_n;
  constexpr int batch = 5;
  std::vector<int> local_counts(static_cast<std::size_t>(runtime.size()));
  const int err = MPI_Allgather(&local_n, 1, MPI_INT, local_counts.data(), 1, MPI_INT, channel::world_comm());
  channel::test::require(err == MPI_SUCCESS, label + " local_n allgather succeeds");

  int global_n = 0;
  int first_global_row = 0;
  for (int rank = 0; rank < runtime.size(); ++rank) {
    if (rank < runtime.rank()) first_global_row += local_counts[static_cast<std::size_t>(rank)];
    global_n += local_counts[static_cast<std::size_t>(rank)];
  }
  const auto count = static_cast<std::size_t>(local_n * batch);

  std::vector<channel::Complex> ds(count);
  std::vector<channel::Complex> dl(count);
  std::vector<channel::Complex> d(count);
  std::vector<channel::Complex> du(count);
  std::vector<channel::Complex> dw(count);
  std::vector<channel::Complex> rhs(count);

  for (int row = 0; row < local_n; ++row) {
    const int global_row = first_global_row + row;
    for (int line = 0; line < batch; ++line) {
      const int p = penta_index(row, line, batch);
      coefficients(global_row, global_n, line, ds[p], dl[p], d[p], du[p], dw[p]);
      rhs[p] = apply_row(global_row, global_n, line);
    }
  }

  channel::DeviceVector<channel::Complex> ds_device(label + "_ds", count);
  channel::DeviceVector<channel::Complex> dl_device(label + "_dl", count);
  channel::DeviceVector<channel::Complex> d_device(label + "_d", count);
  channel::DeviceVector<channel::Complex> du_device(label + "_du", count);
  channel::DeviceVector<channel::Complex> dw_device(label + "_dw", count);
  channel::DeviceVector<channel::Complex> x_device(label + "_x", count);
  ds_device.copy_from_host(ds);
  dl_device.copy_from_host(dl);
  d_device.copy_from_host(d);
  du_device.copy_from_host(du);
  dw_device.copy_from_host(dw);
  x_device.copy_from_host(rhs);

  channel::EndpointSchurYLineSolver solver;
  channel::EndpointSchurYLineConfig cfg;
  cfg.local_n = local_n;
  cfg.batch_count = batch;
  cfg.npy = runtime.size();
  cfg.ipy = runtime.rank();
  cfg.pass_counts = pass_counts;
  cfg.exchange_mode = exchange_mode;
  cfg.comm_y = channel::world_comm();
  solver.prepare(cfg);
  solver.solve(ds_device, dl_device, d_device, du_device, dw_device, x_device);

  const auto solved = x_device.copy_to_host();
  for (int row = 0; row < local_n; ++row) {
    const int global_row = first_global_row + row;
    for (int line = 0; line < batch; ++line) {
      const int p = penta_index(row, line, batch);
      channel::test::require_near(solved[static_cast<std::size_t>(p)], exact_value(global_row, line), 5.0e-11,
                                  label + " endpoint Schur y-line solution rank " +
                                      std::to_string(runtime.rank()) + " row " + std::to_string(row) +
                                      " line " + std::to_string(line));
    }
  }
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1 || runtime.size() == 2 || runtime.size() == 3 || runtime.size() == 4 ||
                             runtime.size() == 6,
                         "endpoint Schur y-line test must run with one, two, three, four, or six MPI ranks");
  const bool force_pass23 = argc > 1 && std::string(argv[1]) == "pass23";
  channel::test::require(!force_pass23 || runtime.size() == 6,
                         "endpoint Schur pass23 variant requires six MPI ranks");

  std::vector<int> pass_counts;
  if (runtime.size() == 2) {
    pass_counts = {2};
  } else if (runtime.size() == 4) {
    pass_counts = {2, 2};
  } else if (force_pass23) {
    pass_counts = {2, 3};
  }

  run_case(runtime, channel::ExchangeMode::AllToAll, "alltoallv_uniform", 7, false, pass_counts);
  run_case(runtime, channel::ExchangeMode::Auto, "auto_exchange_uniform", 7, false, pass_counts);
  run_case(runtime, channel::ExchangeMode::AllToAll, "alltoallv_uneven_local_n", 7, true, pass_counts);
  run_case(runtime, channel::ExchangeMode::Auto, "auto_exchange_uneven_local_n", 7, true, pass_counts);
  run_case(runtime, channel::ExchangeMode::AllToAll, "alltoallv_small_interior", 5, false, pass_counts);
  run_case(runtime, channel::ExchangeMode::Auto, "auto_exchange_small_interior", 5, false, pass_counts);

  if (runtime.rank() == 0) std::cout << "Endpoint Schur y-line test PASSED\n";
  return 0;
}
