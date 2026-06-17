#pragma once

#include "channel/device_vector.hpp"
#include "channel/dnsdata.hpp"
#include "channel/ffts.hpp"
#include "channel/mpi_transpose.hpp"
#include "channel/types.hpp"

#include <array>
#include <memory>
#include <span>
#include <vector>

namespace channel {

enum class VelocityProduct : int {
  UU = 0,
  VV = 1,
  WW = 2,
  UV = 3,
  VW = 4,
  UW = 5,
};

struct DnsNonlinearProductTransformConfig {
  int u_component = 0;
  int v_component = 1;
  int w_component = 2;
  std::array<int, 6> product_components = {5, 6, 7, 8, 9, 10};
  bool enable_velocity_inverse_ffts = true;
  FftNormalization velocity_inverse_normalization = FftNormalization::InverseLength;
  std::vector<DnsComponentTransposeConfig> velocity_transposes;
  double product_factor = 1.0;
  bool enable_product_forward_ffts = true;
  FftNormalization product_forward_normalization = FftNormalization::InverseLength;
  std::vector<DnsComponentTransposeConfig> product_transposes;
};

class DnsNonlinearProductTransformStage {
public:
  void prepare(const DnsState& state, DnsNonlinearProductTransformConfig cfg);
  void apply(DnsState& state);
  void build_products(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsNonlinearProductTransformConfig& config() const { return cfg_; }
  [[nodiscard]] int velocity_fft_count() const { return static_cast<int>(velocity_inverse_ffts_.size()); }
  [[nodiscard]] int velocity_transpose_count() const { return static_cast<int>(velocity_transposes_.size()); }
  [[nodiscard]] int product_fft_count() const { return static_cast<int>(product_forward_ffts_.size()); }
  [[nodiscard]] int product_transpose_count() const { return static_cast<int>(product_transposes_.size()); }

private:
  void check_state(const DnsState& state, const char* caller) const;
  void check_component(int component, const char* name) const;

  DnsGrid grid_;
  DnsNonlinearProductTransformConfig cfg_;
  std::vector<std::unique_ptr<DnsComponentFftPlan>> velocity_inverse_ffts_;
  std::vector<std::unique_ptr<DnsComponentTransposePlan>> velocity_transposes_;
  std::vector<std::unique_ptr<DnsComponentFftPlan>> product_forward_ffts_;
  std::vector<std::unique_ptr<DnsComponentTransposePlan>> product_transposes_;
  bool prepared_ = false;
};

struct DnsNonlinearVelocityRhsConfig {
  int u_component = 0;
  int v_component = 1;
  int w_component = 2;
  int eta_rhs_component = 3;
  int d2v_rhs_component = 4;
  std::array<int, 6> product_components = {5, 6, 7, 8, 9, 10};
  bool build_products_from_velocity = false;
  double product_factor = 1.0;
  double viscosity = 0.0;
  double dt = 1.0;
  double implicit_weight = 1.0;
  double explicit_weight = 1.0;
  double history_weight = 0.0;
  Complex mean_pressure = Complex(0.0, 0.0);
};

class DnsNonlinearVelocityRhsStage {
public:
  void prepare(const DnsState& state, DnsNonlinearVelocityRhsConfig cfg);
  void copy_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                       std::span<const Complex> ibeta,
                                       std::span<const double> k2);
  void copy_y_derivatives_from_host(std::span<const double> derivatives);
  void set_time_scheme(double dt, double implicit_weight, double explicit_weight, double history_weight);
  void reset_history();
  void build_velocity_products(DnsState& state);
  void initialize_rhs(DnsState& state);
  void accumulate_product(DnsState& state, VelocityProduct product);
  void apply(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsNonlinearVelocityRhsConfig& config() const { return cfg_; }
  [[nodiscard]] std::vector<Complex> eta_history_host() const { return old_eta_rhs_.copy_to_host(); }
  [[nodiscard]] std::vector<Complex> d2v_history_host() const { return old_d2v_rhs_.copy_to_host(); }

  static constexpr int derivative_orders = 4;
  static constexpr int derivative_stencil = 5;

private:
  void check_state(const DnsState& state, const char* caller) const;
  void check_prepared(const char* caller) const;
  void check_component(int component, const char* name) const;

  DnsGrid grid_;
  DnsNonlinearVelocityRhsConfig cfg_;
  DeviceVector<Complex> ialfa_;
  DeviceVector<Complex> ibeta_;
  DeviceVector<double> k2_;
  DeviceVector<double> derivatives_;
  DeviceVector<Complex> old_eta_rhs_;
  DeviceVector<Complex> old_d2v_rhs_;
  bool prepared_ = false;
  bool wavenumbers_ready_ = false;
  bool derivatives_ready_ = false;
};

} // namespace channel
