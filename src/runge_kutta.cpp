#include "channel/runge_kutta.hpp"

#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

#include <mpi.h>

namespace channel {
namespace {

bool debug_timestep_enabled() {
  const char* value = std::getenv("CHANNEL_DEBUG_TIMESTEP");
  return value != nullptr && std::string(value) != "0";
}

bool debug_ghosts_enabled() {
  const char* value = std::getenv("CHANNEL_DEBUG_TIMESTEP");
  return value != nullptr && std::string(value).find("ghost") != std::string::npos;
}

bool debug_line_enabled() {
  const char* value = std::getenv("CHANNEL_DEBUG_TIMESTEP");
  return value != nullptr && std::string(value).find("line") != std::string::npos;
}

bool debug_worst_line_enabled() {
  const char* value = std::getenv("CHANNEL_DEBUG_TIMESTEP");
  return value != nullptr && std::string(value).find("worstline") != std::string::npos;
}

int debug_int_env(const char* name, int fallback) {
  const char* value = std::getenv(name);
  if (value == nullptr || *value == '\0') return fallback;
  try {
    return std::stoi(value);
  } catch (...) {
    return fallback;
  }
}

void debug_state_summary(const char* label,
                         int stage,
                         double time,
                         double dt,
                         const DnsState& state,
                         int max_components = 6) {
  if (!debug_timestep_enabled()) return;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const bool print_ghosts = debug_ghosts_enabled();
  const bool print_line = debug_line_enabled();
  const bool print_worst_line = debug_worst_line_enabled();
  const bool print_probe = std::getenv("CHANNEL_DEBUG_TIMESTEP") != nullptr &&
                           std::string(std::getenv("CHANNEL_DEBUG_TIMESTEP")).find("probe") != std::string::npos;
  if (rank != 0 && !print_ghosts && !print_line && !print_worst_line) return;

  if (rank == 0) {
    std::cout << "CPP_DEBUG " << label
              << " stage=" << stage
              << " time=" << time
              << " dt=" << dt;
    const int components = std::min(state.grid().components, max_components);
    for (int component = 0; component < components; ++component) {
      const auto values = state.component_host(component);
      double max_abs = 0.0;
      int nonfinite = 0;
      for (const auto& value : values) {
        const double re = static_cast<double>(value.real());
        const double im = static_cast<double>(value.imag());
        const double mag = std::hypot(re, im);
        if (!std::isfinite(re) || !std::isfinite(im) || !std::isfinite(mag)) {
          ++nonfinite;
        } else {
          max_abs = std::max(max_abs, mag);
        }
      }
      std::cout << " c" << component << "_max=" << max_abs
                << " c" << component << "_bad=" << nonfinite;
    }
    std::cout << '\n';
    if (print_probe && state.grid().components >= 3) {
      const auto& grid = state.grid();
      const int y_storage = 12 - grid.y_first;
      const int z_storage = 14;
      const int x_storage = 5;
      if (y_storage >= 0 && y_storage < grid.ny && z_storage >= 0 && z_storage < grid.nz &&
          x_storage >= 0 && x_storage < grid.nx) {
        const auto c0 = state.component_host(0);
        const auto c1 = state.component_host(1);
        const auto c2 = state.component_host(2);
        const auto p = static_cast<std::size_t>((y_storage * grid.nz + z_storage) * grid.nx + x_storage);
        std::cout << std::scientific << std::setprecision(16)
                  << "CPP_PROBE " << label
                  << " stage=" << stage
                  << " c0=" << c0[p].real() << " " << c0[p].imag()
                  << " c1=" << c1[p].real() << " " << c1[p].imag()
                  << " c2=" << c2[p].real() << " " << c2[p].imag()
                  << '\n';
      }
    }
  }

  int size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  const auto& grid = state.grid();
  if (grid.components < 3 || grid.nx < 1 || grid.nz < 1) return;

  std::array<std::vector<Complex>, 3> components{
      state.component_host(0),
      state.component_host(1),
      state.component_host(2),
  };
  const int line_count = static_cast<int>(state.line_count());
  if (print_worst_line) {
    int worst_component = -1;
    int worst_y = -1;
    int worst_line = -1;
    int worst_bad = 0;
    double worst_mag = -1.0;
    for (int component = 0; component < std::min(grid.components, 3); ++component) {
      const auto& values = components[static_cast<std::size_t>(component)];
      for (int y = 0; y < grid.ny; ++y) {
        for (int line = 0; line < line_count; ++line) {
          const auto& value = values[static_cast<std::size_t>(y * line_count + line)];
          const double re = static_cast<double>(value.real());
          const double im = static_cast<double>(value.imag());
          const bool bad = !std::isfinite(re) || !std::isfinite(im);
          const double mag = bad ? std::numeric_limits<double>::infinity() : std::hypot(re, im);
          if ((bad && worst_bad == 0) || (bad == (worst_bad != 0) && mag > worst_mag)) {
            worst_component = component;
            worst_y = y;
            worst_line = line;
            worst_bad = bad ? 1 : 0;
            worst_mag = mag;
          }
        }
      }
    }
    for (int rank_to_print = 0; rank_to_print < size; ++rank_to_print) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (rank != rank_to_print) continue;
      const int worst_ix = worst_line >= 0 ? worst_line % grid.nx : -1;
      const int worst_iz_storage = worst_line >= 0 ? worst_line / grid.nx : -1;
      const int worst_iz = worst_iz_storage > grid.nz / 2 ? worst_iz_storage - grid.nz : worst_iz_storage;
      std::cout << std::scientific << std::setprecision(16)
                << "CPP_WORST_LINE " << label
                << " stage=" << stage
                << " rank=" << rank
                << " component=" << worst_component
                << " y=" << (worst_y >= 0 ? grid.y_first + worst_y : -999)
                << " ix=" << worst_ix
                << " iz=" << worst_iz
                << " bad=" << worst_bad
                << " mag=" << worst_mag
                << '\n';
      if (worst_line >= 0) {
        for (int y = 0; y < grid.ny; ++y) {
          const int global_y = grid.y_first + y;
          const auto p = static_cast<std::size_t>(y * line_count + worst_line);
          std::cout << std::scientific << std::setprecision(16)
                    << "CPP_WORST_VALUES " << label
                    << " stage=" << stage
                    << " rank=" << rank
                    << " y=" << global_y
                    << " ix=" << worst_ix
                    << " iz=" << worst_iz
                    << " (" << components[0][p].real() << "," << components[0][p].imag() << ")"
                    << " (" << components[1][p].real() << "," << components[1][p].imag() << ")"
                    << " (" << components[2][p].real() << "," << components[2][p].imag() << ")"
                    << '\n';
        }
      }
      std::cout.flush();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  if (print_line) {
    const int requested_z = debug_int_env("CHANNEL_DEBUG_LINE_Z", -3);
    const int x_storage = debug_int_env("CHANNEL_DEBUG_LINE_X", 5);
    const int z_storage = requested_z < 0 ? requested_z + grid.nz : requested_z;
    if (z_storage >= 0 && z_storage < grid.nz && x_storage >= 0 && x_storage < grid.nx) {
      const int line = z_storage * grid.nx + x_storage;
      for (int rank_to_print = 0; rank_to_print < size; ++rank_to_print) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank != rank_to_print) continue;
        for (int y = 0; y < grid.ny; ++y) {
          const int global_y = grid.y_first + y;
          const auto p = static_cast<std::size_t>(y * line_count + line);
          std::cout << std::scientific << std::setprecision(16)
                    << "CPP_LINE " << label
                    << " stage=" << stage
                    << " rank=" << rank
                    << " y=" << global_y
                    << " ix=" << x_storage
                    << " iz=" << requested_z
                    << " (" << components[0][p].real() << "," << components[0][p].imag() << ")"
                    << " (" << components[1][p].real() << "," << components[1][p].imag() << ")"
                    << " (" << components[2][p].real() << "," << components[2][p].imag() << ")"
                    << '\n';
        }
        std::cout.flush();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  if (!print_ghosts) return;

  const int line = 0;
  const int z = 0;
  const int x = 0;
  for (int rank_to_print = 0; rank_to_print < size; ++rank_to_print) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank != rank_to_print) continue;
    for (int y = 0; y < grid.ny; ++y) {
      const int global_y = grid.y_first + y;
      const auto p = static_cast<std::size_t>(y * line_count + line);
      std::cout << std::scientific << std::setprecision(16)
                << "CPP_GHOST " << label
                << " stage=" << stage
                << " rank=" << rank
                << " y=" << global_y
                << " ix=" << x
                << " iz=" << z
                << " (" << components[0][p].real() << "," << components[0][p].imag() << ")"
                << " (" << components[1][p].real() << "," << components[1][p].imag() << ")"
                << " (" << components[2][p].real() << "," << components[2][p].imag() << ")"
                << '\n';
    }
    std::cout.flush();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace

std::array<RungeKuttaWeights, 3> channel_rk3_weights() {
  return {{{120.0 / 32.0, 2.0, 0.0},
           {120.0 / 8.0, 50.0 / 8.0, 34.0 / 8.0},
           {120.0 / 20.0, 90.0 / 20.0, 50.0 / 20.0}}};
}

void DnsRungeKuttaTimestepper::prepare(const DnsState& state, DnsRungeKuttaTimestepperConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  if (cfg.dt <= 0.0) {
    throw std::runtime_error("DnsRungeKuttaTimestepper requires dt > 0");
  }
  if (cfg.stages.empty()) {
    throw std::runtime_error("DnsRungeKuttaTimestepper requires at least one stage");
  }

  grid_ = grid;
  cfg_ = std::move(cfg);
  time_ = cfg_.time;

  nonlinear_product_transform_.reset();
  if (cfg_.enable_nonlinear_product_transform) {
    nonlinear_product_transform_ = std::make_unique<DnsNonlinearProductTransformStage>();
    nonlinear_product_transform_->prepare(state, cfg_.nonlinear_product_transform);
  }

  nonlinear_velocity_rhs_.reset();
  if (cfg_.enable_nonlinear_velocity_rhs) {
    nonlinear_velocity_rhs_ = std::make_unique<DnsNonlinearVelocityRhsStage>();
    auto nonlinear_cfg = cfg_.nonlinear_velocity_rhs;
    nonlinear_cfg.dt = cfg_.dt;
    nonlinear_cfg.implicit_weight = cfg_.stages.front().weights[0];
    nonlinear_cfg.explicit_weight = cfg_.stages.front().weights[1];
    nonlinear_cfg.history_weight = cfg_.stages.front().weights[2];
    nonlinear_velocity_rhs_->prepare(state, std::move(nonlinear_cfg));
  }

  velocity_recovery_.reset();
  if (cfg_.enable_velocity_recovery) {
    velocity_recovery_ = std::make_unique<DnsVelocityRecoveryStage>();
    velocity_recovery_->prepare(state, cfg_.velocity_recovery);
  }

  linear_steps_.clear();
  velocity_mean_corrections_.clear();
  reused_stage_cache_.clear();
  if (cfg_.reuse_linear_stage_storage) {
    linear_steps_.resize(1);
    linear_steps_[0].prepare(state, cfg_.stages.front().linear);
    velocity_mean_corrections_.resize(1);
    if (cfg_.stages.front().enable_velocity_mean_correction) {
      velocity_mean_corrections_[0] = std::make_unique<DnsVelocityMeanCorrectionStep>();
      velocity_mean_corrections_[0]->prepare(state, cfg_.stages.front().velocity_mean_correction);
    }
    reused_stage_cache_.resize(cfg_.stages.size());
    for (auto& cache : reused_stage_cache_) {
      cache.component_coefficients.resize(grid_.components);
      cache.component_ready.assign(static_cast<std::size_t>(grid_.components), false);
    }
  } else {
    linear_steps_.resize(cfg_.stages.size());
    velocity_mean_corrections_.resize(cfg_.stages.size());
    for (std::size_t stage = 0; stage < cfg_.stages.size(); ++stage) {
      linear_steps_[stage].prepare(state, cfg_.stages[stage].linear);
      if (cfg_.stages[stage].enable_velocity_mean_correction) {
        velocity_mean_corrections_[stage] = std::make_unique<DnsVelocityMeanCorrectionStep>();
        velocity_mean_corrections_[stage]->prepare(state, cfg_.stages[stage].velocity_mean_correction);
      }
    }
  }
}

void DnsRungeKuttaTimestepper::copy_nonlinear_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                                         std::span<const Complex> ibeta,
                                                                         std::span<const double> k2) {
  if (!nonlinear_velocity_rhs_) {
    throw std::runtime_error("DnsRungeKuttaTimestepper has no nonlinear velocity RHS stage");
  }
  nonlinear_velocity_rhs_->copy_line_wavenumbers_from_host(ialfa, ibeta, k2);
}

void DnsRungeKuttaTimestepper::copy_nonlinear_y_derivatives_from_host(std::span<const double> derivatives) {
  if (!nonlinear_velocity_rhs_) {
    throw std::runtime_error("DnsRungeKuttaTimestepper has no nonlinear velocity RHS stage");
  }
  nonlinear_velocity_rhs_->copy_y_derivatives_from_host(derivatives);
}

void DnsRungeKuttaTimestepper::copy_recovery_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                                        std::span<const Complex> ibeta,
                                                                        std::span<const double> k2) {
  if (!velocity_recovery_) {
    throw std::runtime_error("DnsRungeKuttaTimestepper has no velocity recovery stage");
  }
  velocity_recovery_->copy_line_wavenumbers_from_host(ialfa, ibeta, k2);
}

void DnsRungeKuttaTimestepper::copy_recovery_y_derivatives_from_host(std::span<const double> derivatives) {
  if (!velocity_recovery_) {
    throw std::runtime_error("DnsRungeKuttaTimestepper has no velocity recovery stage");
  }
  velocity_recovery_->copy_y_derivatives_from_host(derivatives);
}

void DnsRungeKuttaTimestepper::copy_recovery_compact_dvdy_coefficients_from_host(std::span<const Complex> ds,
                                                                                 std::span<const Complex> dl,
                                                                                 std::span<const Complex> d,
                                                                                 std::span<const Complex> du,
                                                                                 std::span<const Complex> dw) {
  if (!velocity_recovery_) {
    throw std::runtime_error("DnsRungeKuttaTimestepper has no velocity recovery stage");
  }
  velocity_recovery_->copy_compact_dvdy_coefficients_from_host(ds, dl, d, du, dw);
}

void DnsRungeKuttaTimestepper::copy_recovery_compact_dvdy_boundary_data_from_host(const DnsYLineBoundaryData& boundary) {
  if (!velocity_recovery_) {
    throw std::runtime_error("DnsRungeKuttaTimestepper has no velocity recovery stage");
  }
  velocity_recovery_->copy_compact_dvdy_boundary_data_from_host(boundary);
}

void DnsRungeKuttaTimestepper::copy_yline_coefficients_from_host(int stage,
                                                                 int component,
                                                                 std::span<const Complex> ds,
                                                                 std::span<const Complex> dl,
                                                                 std::span<const Complex> d,
                                                                 std::span<const Complex> du,
                                                                 std::span<const Complex> dw) {
  check_stage(stage, "DnsRungeKuttaTimestepper::copy_yline_coefficients_from_host");
  if (cfg_.reuse_linear_stage_storage) {
    if (component < 0 || component >= grid_.components) {
      throw std::runtime_error("DnsRungeKuttaTimestepper component index out of range");
    }
    auto& cache = reused_stage_cache_[static_cast<std::size_t>(stage)];
    auto& coefficients = cache.component_coefficients[static_cast<std::size_t>(component)];
    coefficients[0].assign(ds.begin(), ds.end());
    coefficients[1].assign(dl.begin(), dl.end());
    coefficients[2].assign(d.begin(), d.end());
    coefficients[3].assign(du.begin(), du.end());
    coefficients[4].assign(dw.begin(), dw.end());
    cache.component_ready[static_cast<std::size_t>(component)] = true;
    if (stage == 0) {
      linear_steps_[0].copy_yline_coefficients_from_host(component,
                                                         coefficients[0],
                                                         coefficients[1],
                                                         coefficients[2],
                                                         coefficients[3],
                                                         coefficients[4]);
    }
    return;
  }
  linear_steps_[static_cast<std::size_t>(stage)].copy_yline_coefficients_from_host(component, ds, dl, d, du, dw);
}

void DnsRungeKuttaTimestepper::copy_yline_boundary_data_from_host(int stage,
                                                                  int component,
                                                                  const DnsYLineBoundaryData& boundary) {
  check_stage(stage, "DnsRungeKuttaTimestepper::copy_yline_boundary_data_from_host");
  linear_steps_[static_cast<std::size_t>(linear_storage_index(stage))].copy_yline_boundary_data_from_host(component,
                                                                                                         boundary);
}

void DnsRungeKuttaTimestepper::copy_mean_correction_matrix_from_host(int stage,
                                                                     int mean_stage,
                                                                     std::span<const double> matrix_rows) {
  linear_steps_[static_cast<std::size_t>(linear_storage_index(stage))]
      .copy_mean_correction_matrix_from_host(mean_stage, matrix_rows);
}

void DnsRungeKuttaTimestepper::copy_velocity_mean_correction_matrix_from_host(int stage,
                                                                              std::span<const double> matrix_rows) {
  check_stage(stage, "DnsRungeKuttaTimestepper::copy_velocity_mean_correction_matrix_from_host");
  if (cfg_.reuse_linear_stage_storage) {
    auto& cache = reused_stage_cache_[static_cast<std::size_t>(stage)];
    cache.velocity_mean_correction_matrix.assign(matrix_rows.begin(), matrix_rows.end());
    cache.velocity_mean_correction_ready = true;
    if (stage == 0) {
      velocity_mean_corrections_[0]->copy_matrix_from_host(cache.velocity_mean_correction_matrix);
    }
    return;
  }
  auto& correction = velocity_mean_corrections_[static_cast<std::size_t>(stage)];
  if (!correction) {
    throw std::runtime_error("DnsRungeKuttaTimestepper stage has no velocity mean-correction step");
  }
  correction->copy_matrix_from_host(matrix_rows);
}

void DnsRungeKuttaTimestepper::configure_device_velocity_yline_coefficients(std::span<const double> k2,
                                                                            std::span<const double> derivatives,
                                                                            int active_y_global_first,
                                                                            int global_y_count,
                                                                            double viscosity,
                                                                            int eta_component,
                                                                            int v_component) {
  if (!cfg_.reuse_linear_stage_storage) {
    throw std::runtime_error("device velocity y-line coefficient assembly requires reused linear stage storage");
  }
  const auto line_count = grid_.line_count();
  const int active_count = grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count;
  if (k2.size() != line_count ||
      derivatives.size() != static_cast<std::size_t>(active_count) *
                                DnsNonlinearVelocityRhsStage::derivative_orders *
                                DnsNonlinearVelocityRhsStage::derivative_stencil) {
    throw std::runtime_error("DnsRungeKuttaTimestepper::configure_device_velocity_yline_coefficients size mismatch");
  }
  velocity_yline_k2_.resize("velocity_yline_k2", k2.size());
  velocity_yline_derivatives_.resize("velocity_yline_derivatives", derivatives.size());
  velocity_yline_k2_.copy_from_host(k2);
  velocity_yline_derivatives_.copy_from_host(derivatives);
  velocity_yline_active_y_global_first_ = active_y_global_first;
  velocity_yline_global_y_count_ = global_y_count;
  velocity_yline_viscosity_ = viscosity;
  velocity_yline_eta_component_ = eta_component;
  velocity_yline_v_component_ = v_component;
  use_device_velocity_yline_coefficients_ = true;
}

void DnsRungeKuttaTimestepper::set_dt(double dt) {
  if (dt <= 0.0) {
    throw std::runtime_error("DnsRungeKuttaTimestepper requires dt > 0");
  }
  cfg_.dt = dt;
}

void DnsRungeKuttaTimestepper::advance_one_step(DnsState& state) {
  check_state(state, "DnsRungeKuttaTimestepper::advance_one_step");
  for (std::size_t stage = 0; stage < cfg_.stages.size(); ++stage) {
    Kokkos::Profiling::ScopedRegion f90_stage_region("rk_substep");
    const std::string stage_label = "rk_substep_" + std::to_string(stage);
    Kokkos::Profiling::ScopedRegion stage_region(stage_label);
    const auto weights = cfg_.stages[stage].weights;
    if (weights[0] == 0.0) {
      throw std::runtime_error("DnsRungeKuttaTimestepper stage implicit weight must be nonzero");
    }
    time_ += 2.0 / weights[0] * cfg_.dt;
    debug_state_summary("stage_start", static_cast<int>(stage), time_, cfg_.dt, state);
    if (nonlinear_product_transform_ || nonlinear_velocity_rhs_) {
      Kokkos::Profiling::ScopedRegion f90_rhs_region("transform_back_and_build_rhs");
      if (nonlinear_product_transform_) {
        nonlinear_product_transform_->apply(state);
        debug_state_summary("after_product_transform", static_cast<int>(stage), time_, cfg_.dt, state);
        if (cfg_.after_product_transform_observer) {
          cfg_.after_product_transform_observer(static_cast<int>(stage), state);
        }
      }
      if (nonlinear_velocity_rhs_) {
        nonlinear_velocity_rhs_->set_time_scheme(cfg_.dt, weights[0], weights[1], weights[2]);
        nonlinear_velocity_rhs_->apply(state);
        debug_state_summary("after_velocity_rhs", static_cast<int>(stage), time_, cfg_.dt, state);
      }
    }
    if (cfg_.reuse_linear_stage_storage) {
      Kokkos::Profiling::ScopedRegion region("rk_yline_coefficients");
      activate_reused_linear_stage(static_cast<int>(stage));
    }
    {
      Kokkos::Profiling::ScopedRegion region("linsolve_velocity");
      {
        Kokkos::Profiling::ScopedRegion linear_region("cpp_linsolve_linear_step");
        linear_steps_[static_cast<std::size_t>(linear_storage_index(static_cast<int>(stage)))].advance(state);
        debug_state_summary("after_linear_step", static_cast<int>(stage), time_, cfg_.dt, state);
      }
      if (velocity_recovery_) {
        Kokkos::Profiling::ScopedRegion recovery_region("velocity_recovery");
        velocity_recovery_->apply(state);
        debug_state_summary("after_velocity_recovery", static_cast<int>(stage), time_, cfg_.dt, state);
      }
      auto& correction =
          velocity_mean_corrections_[static_cast<std::size_t>(linear_storage_index(static_cast<int>(stage)))];
      if (correction) {
        Kokkos::Profiling::ScopedRegion correction_region("velocity_mean_correction");
        correction->apply(state);
        debug_state_summary("after_velocity_mean_correction", static_cast<int>(stage), time_, cfg_.dt, state);
      }
    }
  }
}

DnsLinearStep& DnsRungeKuttaTimestepper::linear_step(int stage) {
  check_stage(stage, "DnsRungeKuttaTimestepper::linear_step");
  return linear_steps_[static_cast<std::size_t>(linear_storage_index(stage))];
}

const DnsLinearStep& DnsRungeKuttaTimestepper::linear_step(int stage) const {
  check_stage(stage, "DnsRungeKuttaTimestepper::linear_step");
  return linear_steps_[static_cast<std::size_t>(linear_storage_index(stage))];
}

void DnsRungeKuttaTimestepper::check_state(const DnsState& state, const char* caller) const {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DnsRungeKuttaTimestepper::check_stage(int stage, const char* caller) const {
  if (stage < 0 || stage >= static_cast<int>(cfg_.stages.size())) {
    throw std::runtime_error(std::string(caller) + " stage index out of range");
  }
}

int DnsRungeKuttaTimestepper::linear_storage_index(int stage) const {
  check_stage(stage, "DnsRungeKuttaTimestepper::linear_storage_index");
  return cfg_.reuse_linear_stage_storage ? 0 : stage;
}

void DnsRungeKuttaTimestepper::activate_reused_linear_stage(int stage) {
  if (!cfg_.reuse_linear_stage_storage) return;
  auto& cache = reused_stage_cache_[static_cast<std::size_t>(stage)];
  if (use_device_velocity_yline_coefficients_) {
    const auto weights = cfg_.stages[static_cast<std::size_t>(stage)].weights;
    const double lambda = weights[0] / cfg_.dt;
    linear_steps_[0].assemble_velocity_yline_coefficients(velocity_yline_eta_component_,
                                                          velocity_yline_k2_,
                                                          velocity_yline_derivatives_,
                                                          velocity_yline_active_y_global_first_,
                                                          velocity_yline_global_y_count_,
                                                          lambda,
                                                          velocity_yline_viscosity_,
                                                          false);
    linear_steps_[0].assemble_velocity_yline_coefficients(velocity_yline_v_component_,
                                                          velocity_yline_k2_,
                                                          velocity_yline_derivatives_,
                                                          velocity_yline_active_y_global_first_,
                                                          velocity_yline_global_y_count_,
                                                          lambda,
                                                          velocity_yline_viscosity_,
                                                          true);
  } else {
    for (int component : cfg_.stages[static_cast<std::size_t>(stage)].linear.implicit_yline.components) {
      if (component < 0 || component >= grid_.components ||
          !cache.component_ready[static_cast<std::size_t>(component)]) {
        throw std::runtime_error("DnsRungeKuttaTimestepper missing reused y-line coefficients");
      }
      const auto& coefficients = cache.component_coefficients[static_cast<std::size_t>(component)];
      linear_steps_[0].copy_yline_coefficients_from_host(component,
                                                         coefficients[0],
                                                         coefficients[1],
                                                         coefficients[2],
                                                         coefficients[3],
                                                         coefficients[4]);
    }
  }
  if (cfg_.stages[static_cast<std::size_t>(stage)].enable_velocity_mean_correction) {
    if (!cache.velocity_mean_correction_ready) {
      throw std::runtime_error("DnsRungeKuttaTimestepper missing reused velocity mean-correction matrix");
    }
    velocity_mean_corrections_[0]->copy_matrix_from_host(cache.velocity_mean_correction_matrix);
  }
}

} // namespace channel
