#pragma once

#include "channel/dns_solver.hpp"
#include "channel/nonlinear.hpp"
#include "channel/velocity_recovery.hpp"

#include <array>
#include <functional>
#include <memory>
#include <span>
#include <vector>

namespace channel {

using RungeKuttaWeights = std::array<double, 3>;

[[nodiscard]] std::array<RungeKuttaWeights, 3> channel_rk3_weights();

struct DnsRungeKuttaStageConfig {
  RungeKuttaWeights weights = {1.0, 1.0, 0.0};
  DnsLinearStepConfig linear;
  bool enable_velocity_mean_correction = false;
  DnsVelocityMeanCorrectionConfig velocity_mean_correction;
};

struct DnsRungeKuttaTimestepperConfig {
  double dt = 1.0;
  double time = 0.0;
  bool reuse_linear_stage_storage = false;
  bool enable_nonlinear_product_transform = false;
  DnsNonlinearProductTransformConfig nonlinear_product_transform;
  bool enable_nonlinear_velocity_rhs = true;
  DnsNonlinearVelocityRhsConfig nonlinear_velocity_rhs;
  bool enable_velocity_recovery = false;
  DnsVelocityRecoveryConfig velocity_recovery;
  std::vector<DnsRungeKuttaStageConfig> stages;
  std::function<void(int, const DnsState&)> after_product_transform_observer;
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
  void copy_recovery_compact_dvdy_coefficients_from_host(std::span<const Complex> ds,
                                                         std::span<const Complex> dl,
                                                         std::span<const Complex> d,
                                                         std::span<const Complex> du,
                                                         std::span<const Complex> dw);
  void copy_recovery_compact_dvdy_boundary_data_from_host(const DnsYLineBoundaryData& boundary);
  void copy_yline_coefficients_from_host(int stage,
                                         int component,
                                         std::span<const Complex> ds,
                                         std::span<const Complex> dl,
                                         std::span<const Complex> d,
                                         std::span<const Complex> du,
                                         std::span<const Complex> dw);
  void copy_yline_boundary_data_from_host(int stage, int component, const DnsYLineBoundaryData& boundary);
  void copy_mean_correction_matrix_from_host(int stage, int mean_stage, std::span<const double> matrix_rows);
  void copy_velocity_mean_correction_matrix_from_host(int stage, std::span<const double> matrix_rows);
  void configure_device_velocity_yline_coefficients(std::span<const double> k2,
                                                    std::span<const double> derivatives,
                                                    int active_y_global_first,
                                                    int global_y_count,
                                                    double viscosity,
                                                    int eta_component,
                                                    int v_component);
  void set_dt(double dt);
  void advance_one_step(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] double time() const { return time_; }
  [[nodiscard]] double dt() const { return cfg_.dt; }
  [[nodiscard]] int stage_count() const { return static_cast<int>(cfg_.stages.size()); }
  [[nodiscard]] bool has_nonlinear_product_transform() const { return nonlinear_product_transform_ != nullptr; }
  [[nodiscard]] bool has_nonlinear_velocity_rhs() const { return nonlinear_velocity_rhs_ != nullptr; }
  [[nodiscard]] bool has_velocity_recovery() const { return velocity_recovery_ != nullptr; }
  [[nodiscard]] DnsLinearStep& linear_step(int stage);
  [[nodiscard]] const DnsLinearStep& linear_step(int stage) const;

private:
  struct StageCoefficientCache {
    std::vector<std::array<std::vector<Complex>, 5>> component_coefficients;
    std::vector<double> velocity_mean_correction_matrix;
    std::vector<bool> component_ready;
    bool velocity_mean_correction_ready = false;
  };

  void check_state(const DnsState& state, const char* caller) const;
  void check_stage(int stage, const char* caller) const;
  [[nodiscard]] int linear_storage_index(int stage) const;
  void activate_reused_linear_stage(int stage);

  DnsGrid grid_;
  DnsRungeKuttaTimestepperConfig cfg_;
  double time_ = 0.0;
  std::unique_ptr<DnsNonlinearProductTransformStage> nonlinear_product_transform_;
  std::unique_ptr<DnsNonlinearVelocityRhsStage> nonlinear_velocity_rhs_;
  std::unique_ptr<DnsVelocityRecoveryStage> velocity_recovery_;
  std::vector<DnsLinearStep> linear_steps_;
  std::vector<std::unique_ptr<DnsVelocityMeanCorrectionStep>> velocity_mean_corrections_;
  std::vector<StageCoefficientCache> reused_stage_cache_;
  DeviceVector<double> velocity_yline_k2_;
  DeviceVector<double> velocity_yline_derivatives_;
  int velocity_yline_active_y_global_first_ = 0;
  int velocity_yline_global_y_count_ = 0;
  int velocity_yline_eta_component_ = -1;
  int velocity_yline_v_component_ = -1;
  double velocity_yline_viscosity_ = 0.0;
  bool use_device_velocity_yline_coefficients_ = false;
};

} // namespace channel
