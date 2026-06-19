#include "channel/dns_solver.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>

namespace {

int dns_index(int y, int line, int batch) {
  return line + y * batch;
}

channel::Complex fft_value(int y, int line) {
  return channel::Complex(0.17 * y + 0.09 * line, -0.21 * y + 0.13 * line);
}

channel::Complex payload(int source, int target, int value) {
  return {500.0 * source + 50.0 * target + value, -3.0 * source - 0.5 * target - 0.125 * value};
}

channel::Complex exact_value(int global_row, int line) {
  return channel::Complex(0.4 + 0.06 * global_row + 0.02 * line,
                          -0.3 + 0.015 * (global_row - line));
}

void coefficients(int global_row,
                  int global_n,
                  int line,
                  channel::Complex& ds,
                  channel::Complex& dl,
                  channel::Complex& d,
                  channel::Complex& du,
                  channel::Complex& dw) {
  ds = global_row >= 2 ? channel::Complex(-0.012, 0.002 + 0.0005 * line) : channel::Complex(0.0, 0.0);
  dl = global_row >= 1 ? channel::Complex(-0.082, -0.004) : channel::Complex(0.0, 0.0);
  d = channel::Complex(2.6 + 0.012 * global_row, 0.018 + 0.001 * line);
  du = global_row + 1 < global_n ? channel::Complex(-0.071 + 0.001 * line, 0.003) : channel::Complex(0.0, 0.0);
  dw = global_row + 2 < global_n ? channel::Complex(-0.010, -0.002) : channel::Complex(0.0, 0.0);
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

void run_case(const channel::Runtime& runtime, channel::ExchangeMode transpose_mode, const std::string& label) {
  constexpr int nx = 4;
  constexpr int nz = 1;
  constexpr int batch = nx * nz;
  constexpr int ny = 7;
  constexpr int components = 3;
  const int total = ny * batch;
  channel::test::require(total % runtime.size() == 0, label + " values divide across ranks");
  const int values_per_peer = total / runtime.size();

  std::vector<int> local_counts(static_cast<std::size_t>(runtime.size()), ny);
  int global_n = 0;
  int first_global_row = 0;
  for (int rank = 0; rank < runtime.size(); ++rank) {
    if (rank < runtime.rank()) first_global_row += local_counts[static_cast<std::size_t>(rank)];
    global_n += local_counts[static_cast<std::size_t>(rank)];
  }

  channel::DnsState state;
  state.resize({nx, ny, nz, components});

  std::vector<channel::Complex> fft_component(state.values_per_component());
  for (int y = 0; y < ny; ++y) {
    for (int line = 0; line < batch; ++line) {
      fft_component[static_cast<std::size_t>(dns_index(y, line, batch))] = fft_value(y, line);
    }
  }

  std::vector<channel::Complex> transpose_component(state.values_per_component(), channel::Complex(-1.0, -1.0));
  if (transpose_mode == channel::ExchangeMode::AllGather) {
    for (int i = 0; i < values_per_peer; ++i) {
      transpose_component[static_cast<std::size_t>(i)] = payload(runtime.rank(), 0, i);
    }
  } else {
    for (int target = 0; target < runtime.size(); ++target) {
      for (int i = 0; i < values_per_peer; ++i) {
        transpose_component[static_cast<std::size_t>(target * values_per_peer + i)] =
            payload(runtime.rank(), target, i);
      }
    }
  }

  std::vector<channel::Complex> ds(state.values_per_component());
  std::vector<channel::Complex> dl(state.values_per_component());
  std::vector<channel::Complex> d(state.values_per_component());
  std::vector<channel::Complex> du(state.values_per_component());
  std::vector<channel::Complex> dw(state.values_per_component());
  std::vector<channel::Complex> rhs(state.values_per_component());
  for (int row = 0; row < ny; ++row) {
    const int global_row = first_global_row + row;
    for (int line = 0; line < batch; ++line) {
      const int p = dns_index(row, line, batch);
      coefficients(global_row, global_n, line, ds[p], dl[p], d[p], du[p], dw[p]);
      rhs[p] = apply_row(global_row, global_n, line);
    }
  }

  state.copy_component_from_host(0, fft_component);
  state.copy_component_from_host(1, transpose_component);
  state.copy_component_from_host(2, rhs);

  channel::DnsLinearStep step;
  channel::DnsLinearStepConfig cfg;
  cfg.forward_ffts = {{0, channel::FftNormalization::InverseLength}};
  cfg.transposes = {{1, values_per_peer, transpose_mode, channel::world_comm(), label + "_transpose"}};
  cfg.enable_implicit_yline = true;
  cfg.implicit_yline.components = {2};
  cfg.implicit_yline.npy = runtime.size();
  cfg.implicit_yline.ipy = runtime.rank();
  if (runtime.size() == 2) {
    cfg.implicit_yline.pass_counts = {2};
  }
  cfg.implicit_yline.exchange_mode = channel::ExchangeMode::Auto;
  cfg.implicit_yline.comm_y = channel::world_comm();
  cfg.inverse_ffts = {{0, channel::FftNormalization::InverseLength}};
  step.prepare(state, cfg);
  step.copy_yline_coefficients_from_host(2, ds, dl, d, du, dw);
  step.advance(state);

  const auto got_fft = state.component_host(0);
  for (std::size_t i = 0; i < fft_component.size(); ++i) {
    channel::test::require_near(got_fft[i], fft_component[i], 1.0e-11, label + " FFT roundtrip component");
  }

  const auto got_transpose = state.component_host(1);
  for (int source = 0; source < runtime.size(); ++source) {
    for (int i = 0; i < values_per_peer; ++i) {
      const auto want = transpose_mode == channel::ExchangeMode::AllGather ? payload(source, 0, i)
                                                                           : payload(source, runtime.rank(), i);
      channel::test::require_near(got_transpose[static_cast<std::size_t>(source * values_per_peer + i)], want, 0.0,
                                  label + " transposed component");
    }
  }

  const auto solved = state.component_host(2);
  for (int row = 0; row < ny; ++row) {
    const int global_row = first_global_row + row;
    for (int line = 0; line < batch; ++line) {
      const int p = dns_index(row, line, batch);
      channel::test::require_near(solved[static_cast<std::size_t>(p)], exact_value(global_row, line), 8.0e-11,
                                  label + " implicit y-line component");
    }
  }

  channel::test::require(step.forward_fft_count() == 1, label + " forward FFT stage count");
  channel::test::require(step.transpose_count() == 1, label + " transpose stage count");
  channel::test::require(step.inverse_fft_count() == 1, label + " inverse FFT stage count");
  channel::test::require(step.has_implicit_yline(), label + " implicit y-line stage enabled");
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1 || runtime.size() == 2,
                         "DNS linear step test must run with one or two MPI ranks");

  channel::DnsState validation_state;
  validation_state.resize({2, 7, 2, 3});
  channel::DnsLinearStep validation_step;
  channel::DnsLinearStepConfig empty_cfg;
  validation_step.prepare(validation_state, empty_cfg);
  channel::test::require(!validation_step.has_implicit_yline(), "empty linear step has no implicit y-line");
  channel::test::require_throws(
      [&] {
        std::vector<channel::Complex> values(validation_state.values_per_component(), channel::Complex(0.0, 0.0));
        validation_step.copy_yline_coefficients_from_host(0, values, values, values, values, values);
      },
      "coefficient upload without implicit y-line stage rejected");

  run_case(runtime, channel::ExchangeMode::AllToAll, "linear_alltoall");
  run_case(runtime, channel::ExchangeMode::AllGather, "linear_allgather");

  if (runtime.rank() == 0) std::cout << "DNS linear step test PASSED\n";
  return 0;
}
