#pragma once

#include "channel/dnsdata.hpp"
#include "channel/ffts.hpp"
#include "channel/mean_correction.hpp"
#include "channel/mpi_transpose.hpp"
#include "channel/nonlinear.hpp"
#include "channel/velocity_recovery.hpp"
#include "channel/yline.hpp"

#include <memory>
#include <span>
#include <vector>

namespace channel {

struct DnsLinearStepConfig {
  std::vector<DnsComponentFftConfig> forward_ffts;
  std::vector<DnsComponentTransposeConfig> transposes;
  bool enable_implicit_yline = false;
  DnsImplicitYLineStepConfig implicit_yline;
  std::vector<DnsMeanCorrectionConfig> mean_corrections;
  std::vector<DnsComponentFftConfig> inverse_ffts;
};

class DnsLinearStep {
public:
  void prepare(const DnsState& state, DnsLinearStepConfig cfg);
  void copy_yline_coefficients_from_host(int component,
                                         std::span<const Complex> ds,
                                         std::span<const Complex> dl,
                                         std::span<const Complex> d,
                                         std::span<const Complex> du,
                                         std::span<const Complex> dw);
  void copy_mean_correction_matrix_from_host(int stage, std::span<const double> matrix_rows);
  void advance(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsLinearStepConfig& config() const { return cfg_; }
  [[nodiscard]] int forward_fft_count() const { return static_cast<int>(forward_ffts_.size()); }
  [[nodiscard]] int transpose_count() const { return static_cast<int>(transposes_.size()); }
  [[nodiscard]] int mean_correction_count() const { return static_cast<int>(mean_corrections_.size()); }
  [[nodiscard]] int inverse_fft_count() const { return static_cast<int>(inverse_ffts_.size()); }
  [[nodiscard]] bool has_implicit_yline() const { return implicit_yline_ != nullptr; }

private:
  void check_grid(const DnsState& state, const char* caller) const;

  DnsGrid grid_;
  DnsLinearStepConfig cfg_;
  std::vector<std::unique_ptr<DnsComponentFftPlan>> forward_ffts_;
  std::vector<std::unique_ptr<DnsComponentTransposePlan>> transposes_;
  std::unique_ptr<DnsImplicitYLineStep> implicit_yline_;
  std::vector<std::unique_ptr<DnsMeanCorrectionStep>> mean_corrections_;
  std::vector<std::unique_ptr<DnsComponentFftPlan>> inverse_ffts_;
};

struct DnsVelocityStepConfig {
  bool enable_nonlinear_product_transform = false;
  DnsNonlinearProductTransformConfig nonlinear_product_transform;
  bool enable_nonlinear_velocity_rhs = false;
  DnsNonlinearVelocityRhsConfig nonlinear_velocity_rhs;
  DnsLinearStepConfig linear;
  bool enable_velocity_recovery = false;
  DnsVelocityRecoveryConfig velocity_recovery;
};

class DnsVelocityStep {
public:
  void prepare(const DnsState& state, DnsVelocityStepConfig cfg);
  void copy_nonlinear_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                 std::span<const Complex> ibeta,
                                                 std::span<const double> k2);
  void copy_nonlinear_y_derivatives_from_host(std::span<const double> derivatives);
  void copy_recovery_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                std::span<const Complex> ibeta,
                                                std::span<const double> k2);
  void copy_recovery_y_derivatives_from_host(std::span<const double> derivatives);
  void copy_yline_coefficients_from_host(int component,
                                         std::span<const Complex> ds,
                                         std::span<const Complex> dl,
                                         std::span<const Complex> d,
                                         std::span<const Complex> du,
                                         std::span<const Complex> dw);
  void copy_mean_correction_matrix_from_host(int stage, std::span<const double> matrix_rows);
  void advance(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] bool has_nonlinear_product_transform() const { return nonlinear_product_transform_ != nullptr; }
  [[nodiscard]] bool has_nonlinear_velocity_rhs() const { return nonlinear_velocity_rhs_ != nullptr; }
  [[nodiscard]] DnsLinearStep& linear_step() { return linear_; }
  [[nodiscard]] const DnsLinearStep& linear_step() const { return linear_; }

private:
  void check_grid(const DnsState& state, const char* caller) const;

  DnsGrid grid_;
  DnsVelocityStepConfig cfg_;
  std::unique_ptr<DnsNonlinearProductTransformStage> nonlinear_product_transform_;
  std::unique_ptr<DnsNonlinearVelocityRhsStage> nonlinear_velocity_rhs_;
  DnsLinearStep linear_;
  std::unique_ptr<DnsVelocityRecoveryStage> velocity_recovery_;
};

} // namespace channel
