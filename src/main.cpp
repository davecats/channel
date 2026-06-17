#include "channel/dns_solver.hpp"
#include "channel/input.hpp"
#include "channel/runge_kutta.hpp"
#include "channel/runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
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

channel::Complex mean_value(int global_y, int line) {
  return channel::Complex(0.25 + 0.04 * global_y + 0.01 * line,
                          -0.12 + 0.02 * global_y - 0.015 * line);
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

double magnitude(const channel::Complex& value) {
  return std::hypot(static_cast<double>(value.real()), static_cast<double>(value.imag()));
}

double max_allreduce(double local) {
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global;
}

channel::ExchangeMode parse_fixture_mode(const std::string& mode) {
  if (mode == "alltoall") return channel::ExchangeMode::AllToAll;
  if (mode == "allgather") return channel::ExchangeMode::AllGather;
  throw std::runtime_error("unknown linear fixture mode: " + mode);
}

bool has_arg(int argc, char** argv, const std::string& needle) {
  for (int i = 1; i < argc; ++i) {
    if (needle == argv[i]) return true;
  }
  return false;
}

int int_arg(int argc, char** argv, const std::string& needle, int fallback) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (needle == argv[i]) return std::stoi(argv[i + 1]);
  }
  return fallback;
}

std::vector<double> benchmark_derivatives(int ny) {
  std::vector<double> derivatives(static_cast<std::size_t>(ny) *
                                  channel::DnsNonlinearVelocityRhsStage::derivative_orders *
                                  channel::DnsNonlinearVelocityRhsStage::derivative_stencil,
                                  0.0);
  auto index = [](int y, int order, int offset) {
    return (y * channel::DnsNonlinearVelocityRhsStage::derivative_orders + order) *
               channel::DnsNonlinearVelocityRhsStage::derivative_stencil +
           (offset + 2);
  };
  for (int y = 0; y < ny; ++y) {
    derivatives[static_cast<std::size_t>(index(y, 0, 0))] = 1.0;
    derivatives[static_cast<std::size_t>(index(y, 1, -1))] = -0.5;
    derivatives[static_cast<std::size_t>(index(y, 1, 1))] = 0.5;
    derivatives[static_cast<std::size_t>(index(y, 2, -1))] = 1.0;
    derivatives[static_cast<std::size_t>(index(y, 2, 0))] = -2.0;
    derivatives[static_cast<std::size_t>(index(y, 2, 1))] = 1.0;
    derivatives[static_cast<std::size_t>(index(y, 3, -2))] = 1.0;
    derivatives[static_cast<std::size_t>(index(y, 3, -1))] = -4.0;
    derivatives[static_cast<std::size_t>(index(y, 3, 0))] = 6.0;
    derivatives[static_cast<std::size_t>(index(y, 3, 1))] = -4.0;
    derivatives[static_cast<std::size_t>(index(y, 3, 2))] = 1.0;
  }
  return derivatives;
}

void fill_benchmark_velocity(channel::DnsState& state, const channel::ChannelInput& input, int first_global_y) {
  const auto& grid = state.grid();
  const int lines = static_cast<int>(grid.line_count());
  std::vector<channel::Complex> u(grid.values_per_component());
  std::vector<channel::Complex> v(grid.values_per_component());
  std::vector<channel::Complex> w(grid.values_per_component());
  for (int y = 0; y < grid.ny; ++y) {
    const double gy = static_cast<double>(first_global_y + y);
    for (int line = 0; line < lines; ++line) {
      const int p = dns_index(y, line, lines);
      const double phase = 0.013 * gy + 0.007 * static_cast<double>(line);
      u[static_cast<std::size_t>(p)] = channel::Complex(input.velocity.meanflowx + 0.01 * std::sin(phase), 0.0);
      v[static_cast<std::size_t>(p)] = channel::Complex(0.005 * std::cos(0.5 * phase), 0.0);
      w[static_cast<std::size_t>(p)] = channel::Complex(input.velocity.meanflowz + 0.01 * std::sin(0.75 * phase), 0.0);
    }
  }
  state.copy_component_from_host(0, u);
  state.copy_component_from_host(1, v);
  state.copy_component_from_host(2, w);
}

void make_line_wavenumbers(const channel::ChannelInput& input,
                           std::vector<channel::Complex>& ialfa,
                           std::vector<channel::Complex>& ibeta,
                           std::vector<double>& k2) {
  const int nx = input.mesh.nx;
  const int nz = input.mesh.nz;
  const int lines = nx * nz;
  ialfa.assign(static_cast<std::size_t>(lines), channel::Complex(0.0, 0.0));
  ibeta.assign(static_cast<std::size_t>(lines), channel::Complex(0.0, 0.0));
  k2.assign(static_cast<std::size_t>(lines), 0.0);
  for (int iz = 0; iz < nz; ++iz) {
    const int signed_iz = iz <= nz / 2 ? iz : iz - nz;
    for (int ix = 0; ix < nx; ++ix) {
      const int p = iz * nx + ix;
      ialfa[static_cast<std::size_t>(p)] = channel::Complex(0.0, input.mesh.alfa0 * ix);
      ibeta[static_cast<std::size_t>(p)] = channel::Complex(0.0, input.mesh.beta0 * signed_iz);
      k2[static_cast<std::size_t>(p)] =
          (input.mesh.alfa0 * ix) * (input.mesh.alfa0 * ix) +
          (input.mesh.beta0 * signed_iz) * (input.mesh.beta0 * signed_iz);
    }
  }
}

void add_identity_yline(channel::DnsRungeKuttaTimestepper& stepper,
                        const channel::DnsState& state,
                        int stage_count) {
  std::vector<channel::Complex> zero(state.values_per_component(), channel::Complex(0.0, 0.0));
  std::vector<channel::Complex> one(state.values_per_component(), channel::Complex(1.0, 0.0));
  for (int stage = 0; stage < stage_count; ++stage) {
    stepper.copy_yline_coefficients_from_host(stage, 0, zero, zero, one, zero, zero);
    stepper.copy_yline_coefficients_from_host(stage, 1, zero, zero, one, zero, zero);
  }
}

struct BandedCoefficients {
  std::vector<channel::Complex> ds;
  std::vector<channel::Complex> dl;
  std::vector<channel::Complex> d;
  std::vector<channel::Complex> du;
  std::vector<channel::Complex> dw;
};

int derivative_index(int y, int order, int offset) {
  return (y * channel::DnsNonlinearVelocityRhsStage::derivative_orders + order) *
             channel::DnsNonlinearVelocityRhsStage::derivative_stencil +
         (offset + 2);
}

double derivative_value(const std::vector<double>& derivatives, int y, int order, int offset) {
  return derivatives[static_cast<std::size_t>(derivative_index(y, order, offset))];
}

BandedCoefficients assemble_velocity_operator(const channel::DnsState& state,
                                              const channel::LineRange& local_y,
                                              int global_ny,
                                              const std::vector<double>& k2,
                                              const std::vector<double>& derivatives,
                                              double lambda,
                                              double viscosity,
                                              bool biharmonic) {
  const auto& grid = state.grid();
  const int lines = static_cast<int>(grid.line_count());
  BandedCoefficients coeffs;
  coeffs.ds.assign(grid.values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.dl.assign(grid.values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.d.assign(grid.values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.du.assign(grid.values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.dw.assign(grid.values_per_component(), channel::Complex(0.0, 0.0));

  for (int y = 0; y < grid.ny; ++y) {
    const int global_y = local_y.first + y;
    for (int line = 0; line < lines; ++line) {
      const int p = dns_index(y, line, lines);
      if (global_y == 0 || global_y == global_ny - 1) {
        coeffs.d[static_cast<std::size_t>(p)] = channel::Complex(1.0, 0.0);
        continue;
      }

      std::array<channel::Complex, 5> row{};
      const double k2_line = k2[static_cast<std::size_t>(line)];
      for (int offset = -2; offset <= 2; ++offset) {
        if (global_y + offset < 0 || global_y + offset >= global_ny) continue;
        const double d0 = derivative_value(derivatives, y, 0, offset);
        const double d2 = derivative_value(derivatives, y, 2, offset);
        const double d4 = derivative_value(derivatives, y, 3, offset);
        double value = 0.0;
        if (biharmonic) {
          value = lambda * (d2 - k2_line * d0) -
                  viscosity * (d4 - 2.0 * k2_line * d2 + k2_line * k2_line * d0);
        } else {
          value = lambda * d0 - viscosity * (d2 - k2_line * d0);
        }
        row[static_cast<std::size_t>(offset + 2)] = channel::Complex(value, 0.0);
      }
      coeffs.ds[static_cast<std::size_t>(p)] = row[0];
      coeffs.dl[static_cast<std::size_t>(p)] = row[1];
      coeffs.d[static_cast<std::size_t>(p)] = row[2];
      coeffs.du[static_cast<std::size_t>(p)] = row[3];
      coeffs.dw[static_cast<std::size_t>(p)] = row[4];
    }
  }
  return coeffs;
}

void add_velocity_yline(channel::DnsRungeKuttaTimestepper& stepper,
                        const channel::DnsState& state,
                        const channel::LineRange& local_y,
                        int global_ny,
                        const std::vector<double>& k2,
                        const std::vector<double>& derivatives,
                        const std::array<channel::RungeKuttaWeights, 3>& weights,
                        double dt,
                        double viscosity) {
  for (int stage = 0; stage < static_cast<int>(weights.size()); ++stage) {
    const double lambda = weights[static_cast<std::size_t>(stage)][0] / dt;
    auto eta = assemble_velocity_operator(state, local_y, global_ny, k2, derivatives, lambda, viscosity, false);
    auto v = assemble_velocity_operator(state, local_y, global_ny, k2, derivatives, lambda, viscosity, true);
    stepper.copy_yline_coefficients_from_host(stage, 0, eta.ds, eta.dl, eta.d, eta.du, eta.dw);
    stepper.copy_yline_coefficients_from_host(stage, 1, v.ds, v.dl, v.d, v.du, v.dw);
  }
}

int run_parsed_benchmark(const channel::Runtime& runtime, int argc, char** argv) {
  const std::string path = argc > 2 ? argv[2] : "dns.in";
  const bool identity_yline = has_arg(argc, argv, "--identity-yline");
  const auto input = channel::read_channel_input(path);
  if (input.parallel.npy_was_set && input.parallel.npy != runtime.size()) {
    throw std::runtime_error("this C++ runner currently uses all MPI ranks as y ranks; parallel.npy must match mpirun -np");
  }
  const auto local_y = channel::split_range(runtime.rank(), input.mesh.ny, runtime.size());
  constexpr int components = 11;
  channel::DnsState state;
  state.resize({input.mesh.nx, local_y.count, input.mesh.nz, components});
  state.fill(channel::Complex(0.0, 0.0));
  fill_benchmark_velocity(state, input, local_y.first);

  std::vector<channel::Complex> ialfa;
  std::vector<channel::Complex> ibeta;
  std::vector<double> k2;
  make_line_wavenumbers(input, ialfa, ibeta, k2);
  const auto derivatives = benchmark_derivatives(local_y.count);

  channel::DnsNonlinearVelocityRhsConfig nonlinear_cfg;
  nonlinear_cfg.u_component = 0;
  nonlinear_cfg.v_component = 1;
  nonlinear_cfg.w_component = 2;
  nonlinear_cfg.eta_rhs_component = 0;
  nonlinear_cfg.d2v_rhs_component = 1;
  nonlinear_cfg.product_components = {5, 6, 7, 8, 9, 10};
  const int nxd = 3 * (input.mesh.nx + 1) / 2;
  const int nzd = 3 * input.mesh.nz;
  nonlinear_cfg.viscosity = input.velocity.viscosity;
  nonlinear_cfg.mean_pressure = channel::Complex(input.velocity.meanpx, input.velocity.meanpz);

  const auto weights = channel::channel_rk3_weights();
  channel::DnsRungeKuttaTimestepperConfig cfg;
  const double run_dt = input.timestepping.dt > 0.0 ? input.timestepping.dt : 1.0;
  cfg.dt = run_dt;
  cfg.time = input.timestepping.time;
  cfg.enable_nonlinear_product_transform = true;
  cfg.nonlinear_product_transform.u_component = 0;
  cfg.nonlinear_product_transform.v_component = 1;
  cfg.nonlinear_product_transform.w_component = 2;
  cfg.nonlinear_product_transform.product_components = nonlinear_cfg.product_components;
  cfg.nonlinear_product_transform.product_factor = 1.0 / (2.0 * static_cast<double>(nxd) * static_cast<double>(nzd));
  cfg.nonlinear_product_transform.enable_velocity_inverse_ffts = true;
  cfg.nonlinear_product_transform.enable_product_forward_ffts = true;
  cfg.enable_nonlinear_velocity_rhs = true;
  cfg.nonlinear_velocity_rhs = nonlinear_cfg;
  cfg.enable_velocity_recovery = true;
  cfg.velocity_recovery.u_component = 0;
  cfg.velocity_recovery.v_component = 1;
  cfg.velocity_recovery.w_component = 2;
  cfg.velocity_recovery.eta_component = 0;
  cfg.velocity_recovery.dvdy_component = 2;
  cfg.stages.resize(weights.size());
  for (std::size_t stage = 0; stage < weights.size(); ++stage) {
    cfg.stages[stage].weights = weights[stage];
    cfg.stages[stage].linear.enable_implicit_yline = true;
    cfg.stages[stage].linear.implicit_yline.components = {0, 1};
    cfg.stages[stage].linear.implicit_yline.npy = runtime.size();
    cfg.stages[stage].linear.implicit_yline.ipy = runtime.rank();
    cfg.stages[stage].linear.implicit_yline.exchange_mode = channel::ExchangeMode::Auto;
    cfg.stages[stage].linear.implicit_yline.comm_y = channel::world_comm();
  }

  channel::DnsRungeKuttaTimestepper stepper;
  stepper.prepare(state, std::move(cfg));
  stepper.copy_nonlinear_line_wavenumbers_from_host(ialfa, ibeta, k2);
  stepper.copy_nonlinear_y_derivatives_from_host(derivatives);
  stepper.copy_recovery_line_wavenumbers_from_host(ialfa, ibeta, k2);
  stepper.copy_recovery_y_derivatives_from_host(derivatives);
  if (identity_yline) {
    add_identity_yline(stepper, state, static_cast<int>(weights.size()));
  } else {
    add_velocity_yline(stepper, state, local_y, input.mesh.ny, k2, derivatives, weights, run_dt, input.velocity.viscosity);
  }

  const int requested_steps = int_arg(argc, argv, "--steps", input.timestepping.nstep);
  const int steps = std::max(1, requested_steps);
  MPI_Barrier(MPI_COMM_WORLD);
  const double total_start = MPI_Wtime();
  double local_step_sum = 0.0;
  for (int step = 0; step < steps; ++step) {
    Kokkos::Profiling::pushRegion("timestep");
    MPI_Barrier(MPI_COMM_WORLD);
    const double step_start = MPI_Wtime();
    stepper.advance_one_step(state);
    Kokkos::fence();
    const double step_elapsed = MPI_Wtime() - step_start;
    double max_step_elapsed = 0.0;
    MPI_Reduce(&step_elapsed, &max_step_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    local_step_sum += step_elapsed;
    if (runtime.rank() == 0) {
      std::cout << "step " << (step + 1) << "/" << steps
                << " elapsed_s=" << max_step_elapsed
                << " time=" << stepper.time() << '\n';
    }
     Kokkos::Profiling::popRegion();
  }
  const double total_elapsed = MPI_Wtime() - total_start;
  double max_total_elapsed = 0.0;
  double max_step_sum = 0.0;
  MPI_Reduce(&total_elapsed, &max_total_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_step_sum, &max_step_sum, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (runtime.rank() == 0) {
    const auto values = state.values_per_component();
    const double gib = static_cast<double>(state.total_values() * sizeof(channel::Complex)) / (1024.0 * 1024.0 * 1024.0);
    std::cout << "channel C++ parsed run complete\n"
              << "  mode: RK nonlinear velocity"
              << (identity_yline ? " + identity implicit y-line" : " + assembled implicit velocity y-line") << '\n'
              << "  note: pressure/projection and production x/z dealiasing layout are not ported yet\n"
              << "  ranks=" << runtime.size()
              << " global_grid=" << input.mesh.nx << "x" << input.mesh.ny << "x" << input.mesh.nz
              << " local_y=" << local_y.count
              << " components=" << components
              << " values_per_component_per_rank=" << values
              << " state_GiB_per_rank=" << gib << '\n'
              << "  steps=" << steps
              << " avg_step_s=" << (max_step_sum / static_cast<double>(steps))
              << " total_s=" << max_total_elapsed << '\n';
  }
  return 0;
}

int run_linear_fixture(const channel::Runtime& runtime, channel::ExchangeMode transpose_mode) {
  constexpr int nx = 4;
  constexpr int nz = 1;
  constexpr int batch = nx * nz;
  constexpr int ny = 8;
  constexpr int components = 4;
  const int total = ny * batch;
  if (total % runtime.size() != 0) {
    throw std::runtime_error("linear fixture requires local values divisible by MPI ranks");
  }
  const int values_per_peer = total / runtime.size();
  const int first_global_row = runtime.rank() * ny;
  const int global_n = runtime.size() * ny;
  const int global_mean_ny = global_n + 1;

  channel::DnsState state;
  state.resize({nx, ny, nz, components});

  std::vector<channel::Complex> fft_component(state.values_per_component());
  std::vector<channel::Complex> transpose_component(state.values_per_component(), channel::Complex(-1.0, -1.0));
  std::vector<channel::Complex> mean_component(state.values_per_component(), channel::Complex(11.0, -11.0));
  std::vector<channel::Complex> ds(state.values_per_component());
  std::vector<channel::Complex> dl(state.values_per_component());
  std::vector<channel::Complex> d(state.values_per_component());
  std::vector<channel::Complex> du(state.values_per_component());
  std::vector<channel::Complex> dw(state.values_per_component());
  std::vector<channel::Complex> rhs(state.values_per_component());

  for (int y = 0; y < ny; ++y) {
    for (int line = 0; line < batch; ++line) {
      const int p = dns_index(y, line, batch);
      const int global_row = first_global_row + y;
      fft_component[static_cast<std::size_t>(p)] = fft_value(y, line);
      mean_component[static_cast<std::size_t>(p)] = line == 0 ? mean_value(global_row + 1, line)
                                                             : channel::Complex(-13.0 - line, 13.0 + y);
      coefficients(global_row, global_n, line, ds[p], dl[p], d[p], du[p], dw[p]);
      rhs[static_cast<std::size_t>(p)] = apply_row(global_row, global_n, line);
    }
  }

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

  state.copy_component_from_host(0, fft_component);
  state.copy_component_from_host(1, transpose_component);
  state.copy_component_from_host(2, rhs);
  state.copy_component_from_host(3, mean_component);

  std::vector<double> mean_matrix(static_cast<std::size_t>(global_mean_ny + 1) * 5, 0.0);
  for (int iy = 1; iy <= global_mean_ny - 1; ++iy) {
    mean_matrix[static_cast<std::size_t>(iy * 5 + 2)] = 1.0;
  }
  channel::CompactBoundaryRows mean_boundaries;
  mean_boundaries.lower_ghost = {1.0, 0.0, 0.0, 0.0, 0.0};
  mean_boundaries.lower = {0.0, 1.0, 0.0, 0.0, 0.0};
  mean_boundaries.upper = {0.0, 0.0, 0.0, 1.0, 0.0};
  mean_boundaries.upper_ghost = {0.0, 0.0, 0.0, 0.0, 1.0};

  channel::DnsLinearStepConfig cfg;
  cfg.forward_ffts = {{0, channel::FftNormalization::InverseLength}};
  cfg.transposes = {{1, values_per_peer, transpose_mode, channel::world_comm(), "driver_linear_transpose"}};
  cfg.enable_implicit_yline = true;
  cfg.implicit_yline.components = {2};
  cfg.implicit_yline.npy = runtime.size();
  cfg.implicit_yline.ipy = runtime.rank();
  cfg.implicit_yline.exchange_mode = channel::ExchangeMode::Auto;
  cfg.implicit_yline.comm_y = channel::world_comm();
  channel::DnsMeanCorrectionConfig mean_cfg;
  mean_cfg.component = 3;
  mean_cfg.line = 0;
  mean_cfg.root = 0;
  mean_cfg.comm_y = channel::world_comm();
  mean_cfg.boundaries = mean_boundaries;
  mean_cfg.label = "driver_mean_correction";
  cfg.mean_corrections = {mean_cfg};
  cfg.inverse_ffts = {{0, channel::FftNormalization::InverseLength}};

  channel::DnsLinearStep step;
  step.prepare(state, cfg);
  step.copy_yline_coefficients_from_host(2, ds, dl, d, du, dw);
  step.copy_mean_correction_matrix_from_host(0, mean_matrix);
  step.advance(state);

  double local_fft_error = 0.0;
  const auto got_fft = state.component_host(0);
  for (std::size_t i = 0; i < fft_component.size(); ++i) {
    local_fft_error = std::max(local_fft_error, magnitude(got_fft[i] - fft_component[i]));
  }

  double local_transpose_error = 0.0;
  const auto got_transpose = state.component_host(1);
  for (int source = 0; source < runtime.size(); ++source) {
    for (int i = 0; i < values_per_peer; ++i) {
      const auto want = transpose_mode == channel::ExchangeMode::AllGather ? payload(source, 0, i)
                                                                           : payload(source, runtime.rank(), i);
      const auto got = got_transpose[static_cast<std::size_t>(source * values_per_peer + i)];
      local_transpose_error = std::max(local_transpose_error, magnitude(got - want));
    }
  }

  double local_yline_error = 0.0;
  const auto solved = state.component_host(2);
  for (int row = 0; row < ny; ++row) {
    const int global_row = first_global_row + row;
    for (int line = 0; line < batch; ++line) {
      const int p = dns_index(row, line, batch);
      local_yline_error =
          std::max(local_yline_error, magnitude(solved[static_cast<std::size_t>(p)] - exact_value(global_row, line)));
    }
  }

  double local_mean_error = 0.0;
  const auto got_mean = state.component_host(3);
  for (int row = 0; row < ny; ++row) {
    const int global_y = first_global_row + row + 1;
    for (int line = 0; line < batch; ++line) {
      const int p = dns_index(row, line, batch);
      const auto want = line == 0 ? mean_value(global_y, line) : mean_component[static_cast<std::size_t>(p)];
      local_mean_error = std::max(local_mean_error, magnitude(got_mean[static_cast<std::size_t>(p)] - want));
    }
  }

  const double fft_error = max_allreduce(local_fft_error);
  const double transpose_error = max_allreduce(local_transpose_error);
  const double yline_error = max_allreduce(local_yline_error);
  const double mean_error = max_allreduce(local_mean_error);
  if (fft_error > 1.0e-11 || transpose_error > 0.0 || yline_error > 1.0e-10 || mean_error > 1.0e-12) {
    if (runtime.rank() == 0) {
      std::cerr << "linear fixture mismatch: fft=" << fft_error << " transpose=" << transpose_error
                << " yline=" << yline_error << " mean=" << mean_error << '\n';
    }
    return 1;
  }

  if (runtime.rank() == 0) {
    std::cout << "channel linear fixture PASSED"
              << " ranks=" << runtime.size()
              << " fft_error=" << fft_error
              << " transpose_error=" << transpose_error
              << " yline_error=" << yline_error
              << " mean_error=" << mean_error << '\n';
  }
  return 0;
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  const std::string command = argc > 1 ? argv[1] : "";
  if (command == "--linear-fixture") {
    const std::string mode = argc > 2 ? argv[2] : "alltoall";
    return run_linear_fixture(runtime, parse_fixture_mode(mode));
  }
  if (command == "--run") {
    return run_parsed_benchmark(runtime, argc, argv);
  }

  if (runtime.rank() == 0) {
    std::cout << "channel C++/Kokkos port bootstrap\n";
    std::cout << "usage: channel --linear-fixture [alltoall|allgather]\n";
    std::cout << "       channel --run [dns.in] [--steps N] [--identity-yline]\n";
  }
  return 0;
}
