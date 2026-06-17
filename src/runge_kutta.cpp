#include "channel/runge_kutta.hpp"

#include <stdexcept>
#include <string>
#include <utility>

namespace channel {

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
  linear_steps_.resize(cfg_.stages.size());
  for (std::size_t stage = 0; stage < cfg_.stages.size(); ++stage) {
    linear_steps_[stage].prepare(state, cfg_.stages[stage].linear);
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

void DnsRungeKuttaTimestepper::copy_yline_coefficients_from_host(int stage,
                                                                 int component,
                                                                 std::span<const Complex> ds,
                                                                 std::span<const Complex> dl,
                                                                 std::span<const Complex> d,
                                                                 std::span<const Complex> du,
                                                                 std::span<const Complex> dw) {
  linear_step(stage).copy_yline_coefficients_from_host(component, ds, dl, d, du, dw);
}

void DnsRungeKuttaTimestepper::copy_mean_correction_matrix_from_host(int stage,
                                                                     int mean_stage,
                                                                     std::span<const double> matrix_rows) {
  linear_step(stage).copy_mean_correction_matrix_from_host(mean_stage, matrix_rows);
}

void DnsRungeKuttaTimestepper::advance_one_step(DnsState& state) {
  check_state(state, "DnsRungeKuttaTimestepper::advance_one_step");
  for (std::size_t stage = 0; stage < cfg_.stages.size(); ++stage) {
    const auto weights = cfg_.stages[stage].weights;
    if (weights[0] == 0.0) {
      throw std::runtime_error("DnsRungeKuttaTimestepper stage implicit weight must be nonzero");
    }
    time_ += 2.0 / weights[0] * cfg_.dt;
    if (nonlinear_product_transform_) {
      nonlinear_product_transform_->apply(state);
    }
    if (nonlinear_velocity_rhs_) {
      nonlinear_velocity_rhs_->set_time_scheme(cfg_.dt, weights[0], weights[1], weights[2]);
      nonlinear_velocity_rhs_->apply(state);
    }
    linear_steps_[stage].advance(state);
    if (velocity_recovery_) {
      velocity_recovery_->apply(state);
    }
  }
}

DnsLinearStep& DnsRungeKuttaTimestepper::linear_step(int stage) {
  check_stage(stage, "DnsRungeKuttaTimestepper::linear_step");
  return linear_steps_[static_cast<std::size_t>(stage)];
}

const DnsLinearStep& DnsRungeKuttaTimestepper::linear_step(int stage) const {
  check_stage(stage, "DnsRungeKuttaTimestepper::linear_step");
  return linear_steps_[static_cast<std::size_t>(stage)];
}

void DnsRungeKuttaTimestepper::check_state(const DnsState& state, const char* caller) const {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DnsRungeKuttaTimestepper::check_stage(int stage, const char* caller) const {
  if (stage < 0 || stage >= static_cast<int>(linear_steps_.size())) {
    throw std::runtime_error(std::string(caller) + " stage index out of range");
  }
}

} // namespace channel
