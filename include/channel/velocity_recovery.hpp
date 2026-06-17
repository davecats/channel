#pragma once

#include "channel/device_vector.hpp"
#include "channel/dnsdata.hpp"
#include "channel/types.hpp"

#include <span>

namespace channel {

struct DnsVelocityRecoveryConfig {
  int u_component = 0;
  int v_component = 1;
  int w_component = 2;
  int eta_component = 3;
  int dvdy_component = 2;
  bool preserve_zero_mode = true;
};

class DnsVelocityRecoveryStage {
public:
  void prepare(const DnsState& state, DnsVelocityRecoveryConfig cfg);
  void copy_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                       std::span<const Complex> ibeta,
                                       std::span<const double> k2);
  void copy_y_derivatives_from_host(std::span<const double> derivatives);
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
  bool prepared_ = false;
  bool wavenumbers_ready_ = false;
  bool derivatives_ready_ = false;
};

} // namespace channel
