#include "channel/dns_solver.hpp"

#include <string>
#include <stdexcept>
#include <utility>

namespace channel {

void DnsLinearStep::prepare(const DnsState& state, DnsLinearStepConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  grid_ = grid;
  cfg_ = std::move(cfg);

  forward_ffts_.clear();
  forward_ffts_.reserve(cfg_.forward_ffts.size());
  for (const auto& fft_cfg : cfg_.forward_ffts) {
    auto plan = std::make_unique<DnsComponentFftPlan>();
    plan->configure(state, fft_cfg);
    forward_ffts_.push_back(std::move(plan));
  }

  transposes_.clear();
  transposes_.reserve(cfg_.transposes.size());
  for (const auto& transpose_cfg : cfg_.transposes) {
    auto plan = std::make_unique<DnsComponentTransposePlan>();
    plan->configure(state, transpose_cfg);
    transposes_.push_back(std::move(plan));
  }

  implicit_yline_.reset();
  if (cfg_.enable_implicit_yline) {
    implicit_yline_ = std::make_unique<DnsImplicitYLineStep>();
    implicit_yline_->prepare(state, cfg_.implicit_yline);
  }

  mean_corrections_.clear();
  mean_corrections_.reserve(cfg_.mean_corrections.size());
  for (const auto& mean_cfg : cfg_.mean_corrections) {
    auto step = std::make_unique<DnsMeanCorrectionStep>();
    step->prepare(state, mean_cfg);
    mean_corrections_.push_back(std::move(step));
  }

  inverse_ffts_.clear();
  inverse_ffts_.reserve(cfg_.inverse_ffts.size());
  for (const auto& fft_cfg : cfg_.inverse_ffts) {
    auto plan = std::make_unique<DnsComponentFftPlan>();
    plan->configure(state, fft_cfg);
    inverse_ffts_.push_back(std::move(plan));
  }
}

void DnsLinearStep::copy_yline_coefficients_from_host(int component,
                                                      std::span<const Complex> ds,
                                                      std::span<const Complex> dl,
                                                      std::span<const Complex> d,
                                                      std::span<const Complex> du,
                                                      std::span<const Complex> dw) {
  if (!implicit_yline_) {
    throw std::runtime_error("DnsLinearStep has no implicit y-line stage");
  }
  implicit_yline_->copy_coefficients_from_host(component, ds, dl, d, du, dw);
}

void DnsLinearStep::copy_mean_correction_matrix_from_host(int stage, std::span<const double> matrix_rows) {
  if (stage < 0 || stage >= static_cast<int>(mean_corrections_.size())) {
    throw std::runtime_error("DnsLinearStep mean-correction stage index out of range");
  }
  mean_corrections_[static_cast<std::size_t>(stage)]->copy_matrix_from_host(matrix_rows);
}

void DnsLinearStep::advance(DnsState& state) {
  check_grid(state, "DnsLinearStep::advance");
  for (auto& plan : forward_ffts_) {
    plan->execute(state, FftDirection::Forward, "dns_linear_forward_fft");
  }
  for (auto& plan : transposes_) {
    plan->execute(state);
  }
  if (implicit_yline_) {
    implicit_yline_->solve(state);
  }
  for (auto& step : mean_corrections_) {
    step->apply(state);
  }
  for (auto& plan : inverse_ffts_) {
    plan->execute(state, FftDirection::Inverse, "dns_linear_inverse_fft");
  }
}

void DnsLinearStep::check_grid(const DnsState& state, const char* caller) const {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DnsVelocityStep::prepare(const DnsState& state, DnsVelocityStepConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  grid_ = grid;
  cfg_ = std::move(cfg);

  nonlinear_product_transform_.reset();
  if (cfg_.enable_nonlinear_product_transform) {
    nonlinear_product_transform_ = std::make_unique<DnsNonlinearProductTransformStage>();
    nonlinear_product_transform_->prepare(state, cfg_.nonlinear_product_transform);
  }

  nonlinear_velocity_rhs_.reset();
  if (cfg_.enable_nonlinear_velocity_rhs) {
    nonlinear_velocity_rhs_ = std::make_unique<DnsNonlinearVelocityRhsStage>();
    nonlinear_velocity_rhs_->prepare(state, cfg_.nonlinear_velocity_rhs);
  }
  linear_.prepare(state, cfg_.linear);

  velocity_recovery_.reset();
  if (cfg_.enable_velocity_recovery) {
    velocity_recovery_ = std::make_unique<DnsVelocityRecoveryStage>();
    velocity_recovery_->prepare(state, cfg_.velocity_recovery);
  }
}

void DnsVelocityStep::copy_nonlinear_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                                std::span<const Complex> ibeta,
                                                                std::span<const double> k2) {
  if (!nonlinear_velocity_rhs_) {
    throw std::runtime_error("DnsVelocityStep has no nonlinear velocity RHS stage");
  }
  nonlinear_velocity_rhs_->copy_line_wavenumbers_from_host(ialfa, ibeta, k2);
}

void DnsVelocityStep::copy_nonlinear_y_derivatives_from_host(std::span<const double> derivatives) {
  if (!nonlinear_velocity_rhs_) {
    throw std::runtime_error("DnsVelocityStep has no nonlinear velocity RHS stage");
  }
  nonlinear_velocity_rhs_->copy_y_derivatives_from_host(derivatives);
}

void DnsVelocityStep::copy_recovery_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                               std::span<const Complex> ibeta,
                                                               std::span<const double> k2) {
  if (!velocity_recovery_) {
    throw std::runtime_error("DnsVelocityStep has no velocity recovery stage");
  }
  velocity_recovery_->copy_line_wavenumbers_from_host(ialfa, ibeta, k2);
}

void DnsVelocityStep::copy_recovery_y_derivatives_from_host(std::span<const double> derivatives) {
  if (!velocity_recovery_) {
    throw std::runtime_error("DnsVelocityStep has no velocity recovery stage");
  }
  velocity_recovery_->copy_y_derivatives_from_host(derivatives);
}

void DnsVelocityStep::copy_yline_coefficients_from_host(int component,
                                                        std::span<const Complex> ds,
                                                        std::span<const Complex> dl,
                                                        std::span<const Complex> d,
                                                        std::span<const Complex> du,
                                                        std::span<const Complex> dw) {
  linear_.copy_yline_coefficients_from_host(component, ds, dl, d, du, dw);
}

void DnsVelocityStep::copy_mean_correction_matrix_from_host(int stage, std::span<const double> matrix_rows) {
  linear_.copy_mean_correction_matrix_from_host(stage, matrix_rows);
}

void DnsVelocityStep::advance(DnsState& state) {
  check_grid(state, "DnsVelocityStep::advance");
  if (nonlinear_product_transform_) {
    nonlinear_product_transform_->apply(state);
  }
  if (nonlinear_velocity_rhs_) {
    nonlinear_velocity_rhs_->apply(state);
  }
  linear_.advance(state);
  if (velocity_recovery_) {
    velocity_recovery_->apply(state);
  }
}

void DnsVelocityStep::check_grid(const DnsState& state, const char* caller) const {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

} // namespace channel
