#include "channel/nonlinear.hpp"

#include "channel/runtime.hpp"

#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace channel {
namespace {

KOKKOS_INLINE_FUNCTION constexpr int derivative_index(int y, int order, int offset) {
  return (y * DnsNonlinearVelocityRhsStage::derivative_orders + order) *
             DnsNonlinearVelocityRhsStage::derivative_stencil +
         (offset + 2);
}

void check_span_size(std::size_t got, std::size_t expected, const char* name) {
  if (got != expected) {
    throw std::runtime_error(std::string(name) + " size mismatch");
  }
}

} // namespace

void DnsNonlinearProductTransformStage::prepare(const DnsState& state, DnsNonlinearProductTransformConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  if (grid.values_per_component() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsNonlinearProductTransformStage component size exceeds int range");
  }
  grid_ = grid;
  cfg_ = std::move(cfg);

  check_component(cfg_.u_component, "u_component");
  check_component(cfg_.v_component, "v_component");
  check_component(cfg_.w_component, "w_component");
  for (const int component : cfg_.product_components) {
    check_component(component, "product_component");
  }

  velocity_inverse_ffts_.clear();
  if (cfg_.enable_velocity_inverse_ffts) {
    for (const int component : {cfg_.u_component, cfg_.v_component, cfg_.w_component}) {
      auto plan = std::make_unique<DnsComponentFftPlan>();
      plan->configure(state, {component, cfg_.velocity_inverse_normalization});
      velocity_inverse_ffts_.push_back(std::move(plan));
    }
  }

  velocity_transposes_.clear();
  velocity_transposes_.reserve(cfg_.velocity_transposes.size());
  for (const auto& transpose_cfg : cfg_.velocity_transposes) {
    check_component(transpose_cfg.component, "velocity_transpose_component");
    auto plan = std::make_unique<DnsComponentTransposePlan>();
    plan->configure(state, transpose_cfg);
    velocity_transposes_.push_back(std::move(plan));
  }

  product_forward_ffts_.clear();
  if (cfg_.enable_product_forward_ffts) {
    for (const int component : cfg_.product_components) {
      auto plan = std::make_unique<DnsComponentFftPlan>();
      plan->configure(state, {component, cfg_.product_forward_normalization});
      product_forward_ffts_.push_back(std::move(plan));
    }
  }

  product_transposes_.clear();
  product_transposes_.reserve(cfg_.product_transposes.size());
  for (const auto& transpose_cfg : cfg_.product_transposes) {
    check_component(transpose_cfg.component, "product_transpose_component");
    auto plan = std::make_unique<DnsComponentTransposePlan>();
    plan->configure(state, transpose_cfg);
    product_transposes_.push_back(std::move(plan));
  }

  prepared_ = true;
}

void DnsNonlinearProductTransformStage::apply(DnsState& state) {
  check_state(state, "DnsNonlinearProductTransformStage::apply");
  for (auto& plan : velocity_inverse_ffts_) {
    plan->execute(state, FftDirection::Inverse, "dns_nonlinear_velocity_inverse_fft");
  }
  for (auto& plan : velocity_transposes_) {
    plan->execute(state);
  }
  build_products(state);
  for (auto& plan : product_forward_ffts_) {
    plan->execute(state, FftDirection::Forward, "dns_nonlinear_product_forward_fft");
  }
  for (auto& plan : product_transposes_) {
    plan->execute(state);
  }
}

void DnsNonlinearProductTransformStage::check_state(const DnsState& state, const char* caller) const {
  if (!prepared_) {
    throw std::runtime_error(std::string(caller) + " called before prepare");
  }
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DnsNonlinearProductTransformStage::check_component(int component, const char* name) const {
  if (component < 0 || component >= grid_.components) {
    throw std::runtime_error(std::string("DnsNonlinearProductTransformStage ") + name + " out of range");
  }
}

void DnsNonlinearProductTransformStage::build_products(DnsState& state) {
  const int count = static_cast<int>(grid_.values_per_component());
  const auto u = state.component_view(cfg_.u_component);
  const auto v = state.component_view(cfg_.v_component);
  const auto w = state.component_view(cfg_.w_component);
  const auto uu = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::UU)]);
  const auto vv = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::VV)]);
  const auto ww = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::WW)]);
  const auto uv = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::UV)]);
  const auto vw = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::VW)]);
  const auto uw = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::UW)]);
  const double factor = cfg_.product_factor;

  Kokkos::parallel_for(
      "dns_nonlinear_product_transform_build_products",
      Kokkos::RangePolicy<ExecutionSpace>(0, count),
      KOKKOS_LAMBDA(const int i) {
        const double ur = u(i).real();
        const double vr = v(i).real();
        const double wr = w(i).real();
        uu(i) = Complex(factor * ur * ur, 0.0);
        vv(i) = Complex(factor * vr * vr, 0.0);
        ww(i) = Complex(factor * wr * wr, 0.0);
        uv(i) = Complex(factor * ur * vr, 0.0);
        vw(i) = Complex(factor * vr * wr, 0.0);
        uw(i) = Complex(factor * ur * wr, 0.0);
      });
  fence("dns_nonlinear_product_transform_build_products");
}

void DnsNonlinearVelocityRhsStage::prepare(const DnsState& state, DnsNonlinearVelocityRhsConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  if (grid.values_per_component() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage component size exceeds int range");
  }
  grid_ = grid;
  cfg_ = std::move(cfg);

  check_component(cfg_.u_component, "u_component");
  check_component(cfg_.v_component, "v_component");
  check_component(cfg_.w_component, "w_component");
  check_component(cfg_.eta_rhs_component, "eta_rhs_component");
  check_component(cfg_.d2v_rhs_component, "d2v_rhs_component");
  for (const int component : cfg_.product_components) {
    check_component(component, "product_component");
  }
  if (cfg_.dt <= 0.0) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage requires dt > 0");
  }

  const auto line_count = grid_.line_count();
  ialfa_.resize("nonlinear_ialfa", line_count);
  ibeta_.resize("nonlinear_ibeta", line_count);
  k2_.resize("nonlinear_k2", line_count);
  derivatives_.resize("nonlinear_y_derivatives",
                      static_cast<std::size_t>(grid_.ny) * derivative_orders * derivative_stencil);
  old_eta_rhs_.resize("nonlinear_old_eta_rhs", grid_.values_per_component());
  old_d2v_rhs_.resize("nonlinear_old_d2v_rhs", grid_.values_per_component());
  reset_history();
  prepared_ = true;
  wavenumbers_ready_ = false;
  derivatives_ready_ = false;
}

void DnsNonlinearVelocityRhsStage::copy_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                                   std::span<const Complex> ibeta,
                                                                   std::span<const double> k2) {
  check_prepared("DnsNonlinearVelocityRhsStage::copy_line_wavenumbers_from_host");
  check_span_size(ialfa.size(), grid_.line_count(), "ialfa");
  check_span_size(ibeta.size(), grid_.line_count(), "ibeta");
  check_span_size(k2.size(), grid_.line_count(), "k2");
  ialfa_.copy_from_host(ialfa);
  ibeta_.copy_from_host(ibeta);
  k2_.copy_from_host(k2);
  wavenumbers_ready_ = true;
}

void DnsNonlinearVelocityRhsStage::copy_y_derivatives_from_host(std::span<const double> derivatives) {
  check_prepared("DnsNonlinearVelocityRhsStage::copy_y_derivatives_from_host");
  check_span_size(derivatives.size(), derivatives_.size(), "y derivatives");
  derivatives_.copy_from_host(derivatives);
  derivatives_ready_ = true;
}

void DnsNonlinearVelocityRhsStage::set_time_scheme(double dt,
                                                   double implicit_weight,
                                                   double explicit_weight,
                                                   double history_weight) {
  check_prepared("DnsNonlinearVelocityRhsStage::set_time_scheme");
  if (dt <= 0.0) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage requires dt > 0");
  }
  cfg_.dt = dt;
  cfg_.implicit_weight = implicit_weight;
  cfg_.explicit_weight = explicit_weight;
  cfg_.history_weight = history_weight;
}

void DnsNonlinearVelocityRhsStage::reset_history() {
  old_eta_rhs_.fill(Complex(0.0, 0.0));
  old_d2v_rhs_.fill(Complex(0.0, 0.0));
}

void DnsNonlinearVelocityRhsStage::build_velocity_products(DnsState& state) {
  check_state(state, "DnsNonlinearVelocityRhsStage::build_velocity_products");
  const int count = static_cast<int>(grid_.values_per_component());
  const auto u = state.component_view(cfg_.u_component);
  const auto v = state.component_view(cfg_.v_component);
  const auto w = state.component_view(cfg_.w_component);
  const auto uu = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::UU)]);
  const auto vv = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::VV)]);
  const auto ww = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::WW)]);
  const auto uv = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::UV)]);
  const auto vw = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::VW)]);
  const auto uw = state.component_view(cfg_.product_components[static_cast<int>(VelocityProduct::UW)]);
  const double factor = cfg_.product_factor;

  Kokkos::parallel_for(
      "dns_nonlinear_build_velocity_products",
      Kokkos::RangePolicy<ExecutionSpace>(0, count),
      KOKKOS_LAMBDA(const int i) {
        const double ur = u(i).real();
        const double vr = v(i).real();
        const double wr = w(i).real();
        uu(i) = Complex(factor * ur * ur, 0.0);
        vv(i) = Complex(factor * vr * vr, 0.0);
        ww(i) = Complex(factor * wr * wr, 0.0);
        uv(i) = Complex(factor * ur * vr, 0.0);
        vw(i) = Complex(factor * vr * wr, 0.0);
        uw(i) = Complex(factor * ur * wr, 0.0);
      });
  fence("dns_nonlinear_build_velocity_products");
}

void DnsNonlinearVelocityRhsStage::initialize_rhs(DnsState& state) {
  check_state(state, "DnsNonlinearVelocityRhsStage::initialize_rhs");
  if (!wavenumbers_ready_ || !derivatives_ready_) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage metadata was not uploaded");
  }

  const int ny = grid_.ny;
  const int lines = static_cast<int>(grid_.line_count());
  const int count = static_cast<int>(grid_.values_per_component());
  const auto u = state.component_view(cfg_.u_component);
  const auto v = state.component_view(cfg_.v_component);
  const auto w = state.component_view(cfg_.w_component);
  const auto eta = state.component_view(cfg_.eta_rhs_component);
  const auto d2v = state.component_view(cfg_.d2v_rhs_component);
  const auto old_eta = old_eta_rhs_.view();
  const auto old_d2v = old_d2v_rhs_.view();
  const auto ialfa = ialfa_.view();
  const auto ibeta = ibeta_.view();
  const auto k2 = k2_.view();
  const auto der = derivatives_.view();
  const double ni = cfg_.viscosity;
  const double inv_dt_weight = cfg_.implicit_weight / cfg_.dt;
  const double explicit_weight = cfg_.explicit_weight;
  const double history_weight = cfg_.history_weight;
  const Complex mean_pressure = cfg_.mean_pressure;

  Kokkos::parallel_for(
      "dns_nonlinear_initialize_velocity_rhs",
      Kokkos::RangePolicy<ExecutionSpace>(0, count),
      KOKKOS_LAMBDA(const int p) {
        const int y = p / lines;
        const int line = p - y * lines;
        Complex d0v(0.0, 0.0);
        Complex d2v_value(0.0, 0.0);
        Complex d4v_value(0.0, 0.0);
        Complex d0u(0.0, 0.0);
        Complex d2u(0.0, 0.0);
        Complex d0w(0.0, 0.0);
        Complex d2w(0.0, 0.0);

        for (int offset = -2; offset <= 2; ++offset) {
          const int yy = y + offset;
          if (yy < 0 || yy >= ny) continue;
          const int q = yy * lines + line;
          const double d0 = der(derivative_index(y, 0, offset));
          const double d2 = der(derivative_index(y, 2, offset));
          const double d4 = der(derivative_index(y, 3, offset));
          d0v += d0 * v(q);
          d2v_value += d2 * v(q);
          d4v_value += d4 * v(q);
          d0u += d0 * u(q);
          d2u += d2 * u(q);
          d0w += d0 * w(q);
          d2w += d2 * w(q);
        }

        const double k2_line = k2(line);
        const Complex d2v_unknown = d2v_value - k2_line * d0v;
        const Complex d2v_implicit =
            ni * (d4v_value - 2.0 * k2_line * d2v_value + k2_line * k2_line * d0v);
        d2v(p) = inv_dt_weight * d2v_unknown + d2v_implicit - history_weight * old_d2v(p);
        old_d2v(p) = Complex(0.0, 0.0);

        Complex eta_unknown(0.0, 0.0);
        Complex eta_implicit(0.0, 0.0);
        if (k2_line == 0.0) {
          eta_unknown = Complex(d0u.real(), d0w.real());
          eta_implicit = ni * Complex(d2u.real(), d2w.real());
        } else {
          eta_unknown = ibeta(line) * d0u - ialfa(line) * d0w;
          eta_implicit = ni * (ibeta(line) * (d2u - k2_line * d0u) -
                               ialfa(line) * (d2w - k2_line * d0w));
        }
        eta(p) = inv_dt_weight * eta_unknown + eta_implicit - history_weight * old_eta(p);
        if (k2_line == 0.0) eta(p) += explicit_weight * mean_pressure;
        old_eta(p) = Complex(0.0, 0.0);
      });
  fence("dns_nonlinear_initialize_velocity_rhs");
}

void DnsNonlinearVelocityRhsStage::accumulate_product(DnsState& state, VelocityProduct product) {
  check_state(state, "DnsNonlinearVelocityRhsStage::accumulate_product");
  if (!wavenumbers_ready_ || !derivatives_ready_) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage metadata was not uploaded");
  }

  const int product_index = static_cast<int>(product);
  if (product_index < 0 || product_index >= static_cast<int>(cfg_.product_components.size())) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage product index out of range");
  }

  const int ny = grid_.ny;
  const int lines = static_cast<int>(grid_.line_count());
  const int count = static_cast<int>(grid_.values_per_component());
  const auto field = state.component_view(cfg_.product_components[static_cast<std::size_t>(product_index)]);
  const auto eta = state.component_view(cfg_.eta_rhs_component);
  const auto d2v = state.component_view(cfg_.d2v_rhs_component);
  const auto old_eta = old_eta_rhs_.view();
  const auto old_d2v = old_d2v_rhs_.view();
  const auto ialfa = ialfa_.view();
  const auto ibeta = ibeta_.view();
  const auto k2 = k2_.view();
  const auto der = derivatives_.view();
  const double explicit_weight = cfg_.explicit_weight;

  Kokkos::parallel_for(
      "dns_nonlinear_accumulate_velocity_product",
      Kokkos::RangePolicy<ExecutionSpace>(0, count),
      KOKKOS_LAMBDA(const int p) {
        const int y = p / lines;
        const int line = p - y * lines;
        Complex dd0(0.0, 0.0);
        Complex dd1(0.0, 0.0);
        Complex dd2(0.0, 0.0);
        for (int offset = -2; offset <= 2; ++offset) {
          const int yy = y + offset;
          if (yy < 0 || yy >= ny) continue;
          const int q = yy * lines + line;
          dd0 += der(derivative_index(y, 0, offset)) * field(q);
          dd1 += der(derivative_index(y, 1, offset)) * field(q);
          dd2 += der(derivative_index(y, 2, offset)) * field(q);
        }

        const double k2_line = k2(line);
        Complex rhsu(0.0, 0.0);
        Complex rhsw(0.0, 0.0);
        Complex d2v_expl(0.0, 0.0);
        if (product_index == static_cast<int>(VelocityProduct::UU)) {
          rhsu = -ialfa(line) * dd0;
          d2v_expl = ialfa(line) * ialfa(line) * dd1;
        } else if (product_index == static_cast<int>(VelocityProduct::VV)) {
          d2v_expl = k2_line * dd1;
        } else if (product_index == static_cast<int>(VelocityProduct::WW)) {
          rhsw = -ibeta(line) * dd0;
          d2v_expl = ibeta(line) * ibeta(line) * dd1;
        } else if (product_index == static_cast<int>(VelocityProduct::UV)) {
          rhsu = -dd1;
          d2v_expl = ialfa(line) * dd2 + ialfa(line) * k2_line * dd0;
        } else if (product_index == static_cast<int>(VelocityProduct::VW)) {
          rhsw = -dd1;
          d2v_expl = ibeta(line) * dd2 + ibeta(line) * k2_line * dd0;
        } else {
          rhsu = -ibeta(line) * dd0;
          rhsw = -ialfa(line) * dd0;
          d2v_expl = 2.0 * ialfa(line) * ibeta(line) * dd1;
        }

        d2v(p) += explicit_weight * d2v_expl;
        old_d2v(p) += d2v_expl;

        Complex eta_expl = ibeta(line) * rhsu - ialfa(line) * rhsw;
        if (k2_line == 0.0) eta_expl = Complex(rhsu.real(), rhsw.real());
        eta(p) += explicit_weight * eta_expl;
        old_eta(p) += eta_expl;
      });
  fence("dns_nonlinear_accumulate_velocity_product");
}

void DnsNonlinearVelocityRhsStage::apply(DnsState& state) {
  check_state(state, "DnsNonlinearVelocityRhsStage::apply");
  if (cfg_.build_products_from_velocity) {
    build_velocity_products(state);
  }
  initialize_rhs(state);
  for (int product = 0; product < static_cast<int>(cfg_.product_components.size()); ++product) {
    accumulate_product(state, static_cast<VelocityProduct>(product));
  }
}

void DnsNonlinearVelocityRhsStage::check_state(const DnsState& state, const char* caller) const {
  check_prepared(caller);
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DnsNonlinearVelocityRhsStage::check_prepared(const char* caller) const {
  if (!prepared_) {
    throw std::runtime_error(std::string(caller) + " called before prepare");
  }
}

void DnsNonlinearVelocityRhsStage::check_component(int component, const char* name) const {
  if (component < 0 || component >= grid_.components) {
    throw std::runtime_error(std::string("DnsNonlinearVelocityRhsStage ") + name + " out of range");
  }
}

} // namespace channel
