#pragma once

#include "channel/device_vector.hpp"
#include "channel/dnsdata.hpp"
#include "channel/types.hpp"
#include "channel/yline.hpp"

#include <array>
#include <span>
#include <vector>

namespace channel {

struct DnsVelocityRecoveryConfig {
  int u_component = 0;
  int v_component = 1;
  int w_component = 2;
  int eta_component = 3;
  int dvdy_component = 2;
  bool preserve_zero_mode = true;
  bool enable_compact_dvdy = false;
  int npy = 1;
  int ipy = 0;
  std::vector<int> pass_counts;
  ExchangeMode exchange_mode = ExchangeMode::Auto;
  MpiComm comm_y = world_comm();
  std::array<double, 5> lower_boundary_rhs_coeff = {};
  std::array<double, 5> lower_ghost_rhs_coeff = {};
  std::array<double, 5> upper_boundary_rhs_coeff = {};
  std::array<double, 5> upper_ghost_rhs_coeff = {};
};

class DnsVelocityRecoveryStage {
public:
  void prepare(const DnsState& state, DnsVelocityRecoveryConfig cfg);
  void copy_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                       std::span<const Complex> ibeta,
                                       std::span<const double> k2);
  void copy_y_derivatives_from_host(std::span<const double> derivatives);
  void copy_compact_dvdy_coefficients_from_host(std::span<const Complex> ds,
                                                std::span<const Complex> dl,
                                                std::span<const Complex> d,
                                                std::span<const Complex> du,
                                                std::span<const Complex> dw);
  void copy_compact_dvdy_boundary_data_from_host(const DnsYLineBoundaryData& boundary);
  void apply(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsVelocityRecoveryConfig& config() const { return cfg_; }

  static constexpr int derivative_orders = 4;
  static constexpr int derivative_stencil = 5;

private:
  void check_state(const DnsState& state, const char* caller) const;
  void check_prepared(const char* caller) const;
  void check_component(int component, const char* name) const;

  DnsGrid grid_;
  DnsVelocityRecoveryConfig cfg_;
  DeviceVector<Complex> ialfa_;
  DeviceVector<Complex> ibeta_;
  DeviceVector<double> k2_;
  DeviceVector<double> derivatives_;
  DnsYLineComponentSolver compact_dvdy_solver_;
  DnsYLineBoundaryData compact_dvdy_boundary_;
  DeviceVector<Complex> compact_lower_boundary_rhs_;
  DeviceVector<Complex> compact_lower_ghost_rhs_;
  DeviceVector<Complex> compact_upper_boundary_rhs_;
  DeviceVector<Complex> compact_upper_ghost_rhs_;
  bool prepared_ = false;
  bool wavenumbers_ready_ = false;
  bool derivatives_ready_ = false;
  bool compact_dvdy_coefficients_ready_ = false;
  bool compact_dvdy_boundary_ready_ = false;
};

} // namespace channel
