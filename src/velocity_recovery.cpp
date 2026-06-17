#include "channel/velocity_recovery.hpp"

#include "channel/runtime.hpp"

#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace channel {
namespace {

KOKKOS_INLINE_FUNCTION constexpr int derivative_index(int y, int order, int offset) {
  return (y * DnsVelocityRecoveryStage::derivative_orders + order) *
             DnsVelocityRecoveryStage::derivative_stencil +
         (offset + 2);
}

void check_span_size(std::size_t got, std::size_t expected, const char* name) {
  if (got != expected) {
    throw std::runtime_error(std::string(name) + " size mismatch");
  }
}

} // namespace

void DnsVelocityRecoveryStage::prepare(const DnsState& state, DnsVelocityRecoveryConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  if (grid.values_per_component() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsVelocityRecoveryStage component size exceeds int range");
  }
  grid_ = grid;
  cfg_ = std::move(cfg);

  check_component(cfg_.u_component, "u_component");
  check_component(cfg_.v_component, "v_component");
  check_component(cfg_.w_component, "w_component");
  check_component(cfg_.eta_component, "eta_component");
  check_component(cfg_.dvdy_component, "dvdy_component");

  const auto line_count = grid_.line_count();
  ialfa_.resize("velocity_recovery_ialfa", line_count);
  ibeta_.resize("velocity_recovery_ibeta", line_count);
  k2_.resize("velocity_recovery_k2", line_count);
  derivatives_.resize("velocity_recovery_y_derivatives",
                      static_cast<std::size_t>(grid_.ny) * derivative_orders * derivative_stencil);
  prepared_ = true;
  wavenumbers_ready_ = false;
  derivatives_ready_ = false;
}

void DnsVelocityRecoveryStage::copy_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                               std::span<const Complex> ibeta,
                                                               std::span<const double> k2) {
  check_prepared("DnsVelocityRecoveryStage::copy_line_wavenumbers_from_host");
  check_span_size(ialfa.size(), grid_.line_count(), "ialfa");
  check_span_size(ibeta.size(), grid_.line_count(), "ibeta");
  check_span_size(k2.size(), grid_.line_count(), "k2");
  ialfa_.copy_from_host(ialfa);
  ibeta_.copy_from_host(ibeta);
  k2_.copy_from_host(k2);
  wavenumbers_ready_ = true;
}

void DnsVelocityRecoveryStage::copy_y_derivatives_from_host(std::span<const double> derivatives) {
  check_prepared("DnsVelocityRecoveryStage::copy_y_derivatives_from_host");
  check_span_size(derivatives.size(), derivatives_.size(), "y derivatives");
  derivatives_.copy_from_host(derivatives);
  derivatives_ready_ = true;
}

void DnsVelocityRecoveryStage::apply(DnsState& state) {
  check_state(state, "DnsVelocityRecoveryStage::apply");
  if (!wavenumbers_ready_ || !derivatives_ready_) {
    throw std::runtime_error("DnsVelocityRecoveryStage metadata was not uploaded");
  }

  const int ny = grid_.ny;
  const int lines = static_cast<int>(grid_.line_count());
  const int count = static_cast<int>(grid_.values_per_component());
  const auto u = state.component_view(cfg_.u_component);
  const auto v = state.component_view(cfg_.v_component);
  const auto w = state.component_view(cfg_.w_component);
  const auto eta = state.component_view(cfg_.eta_component);
  const auto dvdy_out = state.component_view(cfg_.dvdy_component);
  const auto ialfa = ialfa_.view();
  const auto ibeta = ibeta_.view();
  const auto k2 = k2_.view();
  const auto der = derivatives_.view();
  const bool preserve_zero_mode = cfg_.preserve_zero_mode;

  Kokkos::parallel_for(
      "dns_velocity_recovery",
      Kokkos::RangePolicy<ExecutionSpace>(0, count),
      KOKKOS_LAMBDA(const int p) {
        const int y = p / lines;
        const int line = p - y * lines;
        Complex dvdy(0.0, 0.0);
        for (int offset = -2; offset <= 2; ++offset) {
          const int yy = y + offset;
          if (yy < 0 || yy >= ny) continue;
          dvdy += der(derivative_index(y, 1, offset)) * v(yy * lines + line);
        }
        dvdy_out(p) = dvdy;

        const double k2_line = k2(line);
        if (k2_line == 0.0 && preserve_zero_mode) return;
        if (k2_line == 0.0) {
          u(p) = Complex(0.0, 0.0);
          w(p) = Complex(0.0, 0.0);
          return;
        }

        const Complex recovered_u = (ialfa(line) * dvdy - ibeta(line) * eta(p)) / k2_line;
        const Complex recovered_w = (ibeta(line) * dvdy + ialfa(line) * eta(p)) / k2_line;
        u(p) = recovered_u;
        w(p) = recovered_w;
      });
  fence("dns_velocity_recovery");
}

void DnsVelocityRecoveryStage::check_state(const DnsState& state, const char* caller) const {
  check_prepared(caller);
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DnsVelocityRecoveryStage::check_prepared(const char* caller) const {
  if (!prepared_) {
    throw std::runtime_error(std::string(caller) + " called before prepare");
  }
}

void DnsVelocityRecoveryStage::check_component(int component, const char* name) const {
  if (component < 0 || component >= grid_.components) {
    throw std::runtime_error(std::string("DnsVelocityRecoveryStage ") + name + " out of range");
  }
}

} // namespace channel
