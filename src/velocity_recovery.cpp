#include "channel/velocity_recovery.hpp"

#include "channel/runtime.hpp"

#include <Kokkos_Profiling_ScopedRegion.hpp>

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
                      static_cast<std::size_t>(grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count) *
                          derivative_orders * derivative_stencil);
  if (cfg_.enable_compact_dvdy) {
    compact_lower_boundary_rhs_.resize("velocity_recovery_lower_boundary_rhs", line_count);
    compact_lower_ghost_rhs_.resize("velocity_recovery_lower_ghost_rhs", line_count);
    compact_upper_boundary_rhs_.resize("velocity_recovery_upper_boundary_rhs", line_count);
    compact_upper_ghost_rhs_.resize("velocity_recovery_upper_ghost_rhs", line_count);
    DnsYLineComponentConfig solver_cfg;
    solver_cfg.component = cfg_.dvdy_component;
    solver_cfg.npy = cfg_.npy;
    solver_cfg.ipy = cfg_.ipy;
    solver_cfg.pass_counts = cfg_.pass_counts;
    solver_cfg.exchange_mode = cfg_.exchange_mode;
    solver_cfg.comm_y = cfg_.comm_y;
    solver_cfg.fill_interface_ghosts = true;
    compact_dvdy_solver_.prepare(state, std::move(solver_cfg));
  }
  prepared_ = true;
  wavenumbers_ready_ = false;
  derivatives_ready_ = false;
  compact_dvdy_coefficients_ready_ = !cfg_.enable_compact_dvdy;
  compact_dvdy_boundary_ready_ = !cfg_.enable_compact_dvdy;
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

void DnsVelocityRecoveryStage::copy_compact_dvdy_coefficients_from_host(std::span<const Complex> ds,
                                                                        std::span<const Complex> dl,
                                                                        std::span<const Complex> d,
                                                                        std::span<const Complex> du,
                                                                        std::span<const Complex> dw) {
  check_prepared("DnsVelocityRecoveryStage::copy_compact_dvdy_coefficients_from_host");
  if (!cfg_.enable_compact_dvdy) {
    throw std::runtime_error("DnsVelocityRecoveryStage compact dvdy is not enabled");
  }
  compact_dvdy_solver_.copy_coefficients_from_host(ds, dl, d, du, dw);
  compact_dvdy_coefficients_ready_ = true;
}

void DnsVelocityRecoveryStage::copy_compact_dvdy_boundary_data_from_host(const DnsYLineBoundaryData& boundary) {
  check_prepared("DnsVelocityRecoveryStage::copy_compact_dvdy_boundary_data_from_host");
  if (!cfg_.enable_compact_dvdy) {
    throw std::runtime_error("DnsVelocityRecoveryStage compact dvdy is not enabled");
  }
  compact_dvdy_boundary_ = boundary;
  compact_dvdy_solver_.copy_boundary_data_from_host(boundary);
  compact_dvdy_boundary_ready_ = true;
}

void DnsVelocityRecoveryStage::apply(DnsState& state) {
  check_state(state, "DnsVelocityRecoveryStage::apply");
  if (!wavenumbers_ready_ || !derivatives_ready_) {
    throw std::runtime_error("DnsVelocityRecoveryStage metadata was not uploaded");
  }
  if (!compact_dvdy_coefficients_ready_ || !compact_dvdy_boundary_ready_) {
    throw std::runtime_error("DnsVelocityRecoveryStage compact dvdy metadata was not uploaded");
  }

  const int ny = grid_.ny;
  const int nx = grid_.nx;
  const int active_first = grid_.active_y_storage_first();
  const int active_last = grid_.active_y_storage_last_exclusive();
  const auto v = state.component_view_3d(cfg_.v_component);
  const auto dvdy_out = state.component_view_3d(cfg_.dvdy_component);
  const auto der = derivatives_.view();

  if (cfg_.enable_compact_dvdy) {
    {
      Kokkos::Profiling::ScopedRegion region("velocity_recovery compact_dvdy_rhs");
      Kokkos::parallel_for(
          "dns_velocity_recovery_compact_dvdy_rhs",
          Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({active_first, 0, 0}, {active_last, grid_.nz, grid_.nx}),
          KOKKOS_LAMBDA(const int y, const int z, const int x) {
            const int dy = y - active_first;
            Complex rhs(0.0, 0.0);
            for (int offset = -2; offset <= 2; ++offset) {
              const int yy = y + offset;
              rhs += der(derivative_index(dy, 1, offset)) * v(yy, z, x);
            }
            dvdy_out(y, z, x) = rhs;
          });
      fence("dns_velocity_recovery_compact_dvdy_rhs");
    }
    const int line_count = static_cast<int>(grid_.line_count());
    DnsYLineDeviceBoundaryData boundary;
    boundary.lower_enabled = compact_dvdy_boundary_.lower_enabled;
    boundary.upper_enabled = compact_dvdy_boundary_.upper_enabled;
    boundary.lower_ghost_eq = compact_dvdy_boundary_.lower_ghost_eq;
    boundary.lower_boundary_eq = compact_dvdy_boundary_.lower_boundary_eq;
    boundary.upper_boundary_eq = compact_dvdy_boundary_.upper_boundary_eq;
    boundary.upper_ghost_eq = compact_dvdy_boundary_.upper_ghost_eq;
    boundary.lower_boundary_rhs = &compact_lower_boundary_rhs_;
    boundary.lower_ghost_rhs = &compact_lower_ghost_rhs_;
    boundary.upper_boundary_rhs = &compact_upper_boundary_rhs_;
    boundary.upper_ghost_rhs = &compact_upper_ghost_rhs_;

    if (boundary.lower_enabled) {
      Kokkos::Profiling::ScopedRegion region("velocity_recovery lower_boundary_rhs");
      auto lower_boundary = compact_lower_boundary_rhs_.view();
      auto lower_ghost = compact_lower_ghost_rhs_.view();
      Kokkos::Array<double, 5> lower_boundary_coeff{};
      Kokkos::Array<double, 5> lower_ghost_coeff{};
      for (int i = 0; i < 5; ++i) {
        lower_boundary_coeff[static_cast<std::size_t>(i)] = cfg_.lower_boundary_rhs_coeff[static_cast<std::size_t>(i)];
        lower_ghost_coeff[static_cast<std::size_t>(i)] = cfg_.lower_ghost_rhs_coeff[static_cast<std::size_t>(i)];
      }
      const int y_first = grid_.y_first;
      const int nx_local = grid_.nx;
      Kokkos::parallel_for(
          "dns_velocity_recovery_lower_boundary_rhs",
          Kokkos::RangePolicy<ExecutionSpace>(0, line_count),
          KOKKOS_LAMBDA(const int line) {
            const int z = line / nx_local;
            const int x = line - z * nx_local;
            Complex boundary_value(0.0, 0.0);
            Complex ghost_value(0.0, 0.0);
            for (int i = 0; i < 5; ++i) {
              const int y_storage = (-1 + i) - y_first;
              boundary_value += lower_boundary_coeff[static_cast<std::size_t>(i)] * v(y_storage, z, x);
              ghost_value += lower_ghost_coeff[static_cast<std::size_t>(i)] * v(y_storage, z, x);
            }
            lower_boundary(line) = boundary_value;
            lower_ghost(line) = ghost_value;
          });
      fence("dns_velocity_recovery_lower_boundary_rhs");
    }
    if (boundary.upper_enabled) {
      Kokkos::Profiling::ScopedRegion region("velocity_recovery upper_boundary_rhs");
      auto upper_boundary = compact_upper_boundary_rhs_.view();
      auto upper_ghost = compact_upper_ghost_rhs_.view();
      Kokkos::Array<double, 5> upper_boundary_coeff{};
      Kokkos::Array<double, 5> upper_ghost_coeff{};
      for (int i = 0; i < 5; ++i) {
        upper_boundary_coeff[static_cast<std::size_t>(i)] = cfg_.upper_boundary_rhs_coeff[static_cast<std::size_t>(i)];
        upper_ghost_coeff[static_cast<std::size_t>(i)] = cfg_.upper_ghost_rhs_coeff[static_cast<std::size_t>(i)];
      }
      const int first = grid_.active_y_first + (grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count) - 3;
      const int y_first = grid_.y_first;
      const int nx_local = grid_.nx;
      Kokkos::parallel_for(
          "dns_velocity_recovery_upper_boundary_rhs",
          Kokkos::RangePolicy<ExecutionSpace>(0, line_count),
          KOKKOS_LAMBDA(const int line) {
            const int z = line / nx_local;
            const int x = line - z * nx_local;
            Complex boundary_value(0.0, 0.0);
            Complex ghost_value(0.0, 0.0);
            for (int i = 0; i < 5; ++i) {
              const int y_storage = (first + i) - y_first;
              boundary_value += upper_boundary_coeff[static_cast<std::size_t>(i)] * v(y_storage, z, x);
              ghost_value += upper_ghost_coeff[static_cast<std::size_t>(i)] * v(y_storage, z, x);
            }
            upper_boundary(line) = boundary_value;
            upper_ghost(line) = ghost_value;
          });
      fence("dns_velocity_recovery_upper_boundary_rhs");
    }
    {
      Kokkos::Profiling::ScopedRegion region("velocity_recovery compact_boundary_copy");
      compact_dvdy_solver_.copy_boundary_data_from_device(boundary);
    }
    {
      Kokkos::Profiling::ScopedRegion f90_region("linsolve d_v_dy");
      Kokkos::Profiling::ScopedRegion region("velocity_recovery compact_dvdy_solve");
      compact_dvdy_solver_.solve(state);
    }
  } else {
    Kokkos::Profiling::ScopedRegion f90_region("linsolve d_v_dy");
    Kokkos::Profiling::ScopedRegion region("velocity_recovery explicit_dvdy");
    Kokkos::parallel_for(
        "dns_velocity_recovery_explicit_dvdy",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({active_first, 0, 0}, {active_last, grid_.nz, grid_.nx}),
        KOKKOS_LAMBDA(const int y, const int z, const int x) {
          const int dy = y - active_first;
          Complex dvdy(0.0, 0.0);
          for (int offset = -2; offset <= 2; ++offset) {
            const int yy = y + offset;
            if (yy < 0 || yy >= ny) continue;
            dvdy += der(derivative_index(dy, 1, offset)) * v(yy, z, x);
          }
          dvdy_out(y, z, x) = dvdy;
        });
    fence("dns_velocity_recovery_explicit_dvdy");
  }

  const auto u = state.component_view_3d(cfg_.u_component);
  const auto w = state.component_view_3d(cfg_.w_component);
  const auto eta = state.component_view_3d(cfg_.eta_component);
  const auto ialfa = ialfa_.view();
  const auto ibeta = ibeta_.view();
  const auto k2 = k2_.view();
  const bool preserve_zero_mode = cfg_.preserve_zero_mode;
  const int recovery_first = cfg_.enable_compact_dvdy ? 0 : active_first;
  const int recovery_last = cfg_.enable_compact_dvdy ? grid_.ny : active_last;

  {
    Kokkos::Profiling::ScopedRegion f90_region("linsolve recover_u_w");
    Kokkos::Profiling::ScopedRegion region("velocity_recovery recover_u_w");
    Kokkos::parallel_for(
        "dns_velocity_recovery",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({recovery_first, 0, 0}, {recovery_last, grid_.nz, grid_.nx}),
        KOKKOS_LAMBDA(const int y, const int z, const int x) {
          const int line = z * nx + x;
          const Complex dvdy = dvdy_out(y, z, x);

          const double k2_line = k2(line);
          if (k2_line == 0.0 && preserve_zero_mode) return;
          if (k2_line == 0.0) {
            u(y, z, x) = Complex(0.0, 0.0);
            w(y, z, x) = Complex(0.0, 0.0);
            return;
          }

          const Complex recovered_u = (ialfa(line) * dvdy - ibeta(line) * eta(y, z, x)) / k2_line;
          const Complex recovered_w = (ibeta(line) * dvdy + ialfa(line) * eta(y, z, x)) / k2_line;
          u(y, z, x) = recovered_u;
          w(y, z, x) = recovered_w;
        });
    fence("dns_velocity_recovery");
  }
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
