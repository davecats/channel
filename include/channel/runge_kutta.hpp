#pragma once

#include "channel/dns_solver.hpp"
#include "channel/nonlinear.hpp"
#include "channel/velocity_recovery.hpp"

#include <array>
#include <memory>
#include <span>
#include <vector>

namespace channel {

using RungeKuttaWeights = std::array<double, 3>;

[[nodiscard]] std::array<RungeKuttaWeights, 3> channel_rk3_weights();

struct DnsRungeKuttaStageConfig {
  RungeKuttaWeights weights = {1.0, 1.0, 0.0};
  DnsLinearStepConfig linear;
};

struct DnsRungeKuttaTimestepperConfig {
  double dt = 1.0;
  double time = 0.0;
  bool enable_nonlinear_product_transform = false;
  DnsNonlinearProductTransformConfig nonlinear_product_transform;
  bool enable_nonlinear_velocity_rhs = true;
  DnsNonlinearVelocityRhsConfig nonlinear_velocity_rhs;
  bool enable_velocity_recovery = false;
  DnsVelocityRecoveryConfig velocity_recovery;
  std::vector<DnsRungeKuttaStageConfig> stages;
};

class DnsRungeKuttaTimestepper {
public:
  void prepare(const DnsState& state, DnsRungeKuttaTimestepperConfig cfg);
  void copy_nonlinear_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                 std::span<const Complex> ibeta,
                                                 std::span<const double> k2);
  void copy_nonlinear_y_derivatives_from_host(std::span<const double> derivatives);
  void copy_recovery_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                std::span<const Complex> ibeta,
                                                std::span<const double> k2);
  void copy_recovery_y_derivatives_from_host(std::span<const double> derivatives);
  void copy_yline_coefficients_from_host(int stage,
                                         int component,
                                         std::span<const Complex> ds,
                                         std::span<const Complex> dl,
                                         std::span<const Complex> d,
                                         std::span<const Complex> du,
                                         std::span<const Complex> dw);
  void copy_mean_correction_matrix_from_host(int stage, int mean_stage, std::span<const double> matrix_rows);
  void advance_one_step(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] double time() const { return time_; }
  [[nodiscard]] double dt() const { return cfg_.dt; }
  [[nodiscard]] int stage_count() const { return static_cast<int>(linear_steps_.size()); }
  [[nodiscard]] bool has_nonlinear_product_transform() const { return nonlinear_product_transform_ != nullptr; }
  [[nodiscard]] bool has_nonlinear_velocity_rhs() const { return nonlinear_velocity_rhs_ != nullptr; }
  [[nodiscard]] bool has_velocity_recovery() const { return velocity_recovery_ != nullptr; }
  [[nodiscard]] DnsLinearStep& linear_step(int stage);
  [[nodiscard]] const DnsLinearStep& linear_step(int stage) const;

private:
  void check_state(const DnsState& state, const char* caller) const;
  void check_stage(int stage, const char* caller) const;

  DnsGrid grid_;
  DnsRungeKuttaTimestepperConfig cfg_;
  double time_ = 0.0;
  std::unique_ptr<DnsNonlinearProductTransformStage> nonlinear_product_transform_;
  std::unique_ptr<DnsNonlinearVelocityRhsStage> nonlinear_velocity_rhs_;
  std::unique_ptr<DnsVelocityRecoveryStage> velocity_recovery_;
  std::vector<DnsLinearStep> linear_steps_;
};

} // namespace channel
