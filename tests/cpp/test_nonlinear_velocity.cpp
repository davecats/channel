#include "channel/dns_solver.hpp"
#include "channel/nonlinear.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <array>
#include <iostream>
#include <string>
#include <vector>

namespace {

int value_index(int y, int line, int lines) {
  return y * lines + line;
}

int derivative_index(int y, int order, int offset) {
  return (y * channel::DnsNonlinearVelocityRhsStage::derivative_orders + order) *
             channel::DnsNonlinearVelocityRhsStage::derivative_stencil +
         (offset + 2);
}

channel::Complex dd(const std::vector<channel::Complex>& values,
                    const std::vector<double>& derivatives,
                    int ny,
                    int lines,
                    int y,
                    int line,
                    int order) {
  channel::Complex sum(0.0, 0.0);
  for (int offset = -2; offset <= 2; ++offset) {
    const int yy = y + offset;
    if (yy < 0 || yy >= ny) continue;
    sum += derivatives[static_cast<std::size_t>(derivative_index(y, order, offset))] *
           values[static_cast<std::size_t>(value_index(yy, line, lines))];
  }
  return sum;
}

struct ReferenceResult {
  std::vector<channel::Complex> eta;
  std::vector<channel::Complex> d2v;
  std::vector<channel::Complex> old_eta;
  std::vector<channel::Complex> old_d2v;
};

ReferenceResult apply_reference(const channel::DnsGrid& grid,
                                const channel::DnsNonlinearVelocityRhsConfig& cfg,
                                const std::vector<channel::Complex>& ialfa,
                                const std::vector<channel::Complex>& ibeta,
                                const std::vector<double>& k2,
                                const std::vector<double>& derivatives,
                                const std::vector<channel::Complex>& u,
                                const std::vector<channel::Complex>& v,
                                const std::vector<channel::Complex>& w,
                                const std::array<std::vector<channel::Complex>, 6>& products,
                                const std::vector<channel::Complex>& old_eta_in,
                                const std::vector<channel::Complex>& old_d2v_in) {
  const int lines = static_cast<int>(grid.line_count());
  ReferenceResult result;
  result.eta.assign(grid.values_per_component(), channel::Complex(0.0, 0.0));
  result.d2v.assign(grid.values_per_component(), channel::Complex(0.0, 0.0));
  result.old_eta.assign(grid.values_per_component(), channel::Complex(0.0, 0.0));
  result.old_d2v.assign(grid.values_per_component(), channel::Complex(0.0, 0.0));

  for (int y = 0; y < grid.ny; ++y) {
    for (int line = 0; line < lines; ++line) {
      const int p = value_index(y, line, lines);
      const channel::Complex d0v = dd(v, derivatives, grid.ny, lines, y, line, 0);
      const channel::Complex d2v_value = dd(v, derivatives, grid.ny, lines, y, line, 2);
      const channel::Complex d4v_value = dd(v, derivatives, grid.ny, lines, y, line, 3);
      const channel::Complex d0u = dd(u, derivatives, grid.ny, lines, y, line, 0);
      const channel::Complex d2u = dd(u, derivatives, grid.ny, lines, y, line, 2);
      const channel::Complex d0w = dd(w, derivatives, grid.ny, lines, y, line, 0);
      const channel::Complex d2w = dd(w, derivatives, grid.ny, lines, y, line, 2);

      const double k2_line = k2[static_cast<std::size_t>(line)];
      const channel::Complex d2v_unknown = d2v_value - k2_line * d0v;
      const channel::Complex d2v_implicit =
          cfg.viscosity * (d4v_value - 2.0 * k2_line * d2v_value + k2_line * k2_line * d0v);
      result.d2v[static_cast<std::size_t>(p)] =
          (cfg.implicit_weight / cfg.dt) * d2v_unknown + d2v_implicit -
          cfg.history_weight * old_d2v_in[static_cast<std::size_t>(p)];

      channel::Complex eta_unknown(0.0, 0.0);
      channel::Complex eta_implicit(0.0, 0.0);
      if (k2_line == 0.0) {
        eta_unknown = channel::Complex(d0u.real(), d0w.real());
        eta_implicit = cfg.viscosity * channel::Complex(d2u.real(), d2w.real());
      } else {
        eta_unknown = ibeta[static_cast<std::size_t>(line)] * d0u -
                      ialfa[static_cast<std::size_t>(line)] * d0w;
        eta_implicit = cfg.viscosity *
                       (ibeta[static_cast<std::size_t>(line)] * (d2u - k2_line * d0u) -
                        ialfa[static_cast<std::size_t>(line)] * (d2w - k2_line * d0w));
      }
      result.eta[static_cast<std::size_t>(p)] =
          (cfg.implicit_weight / cfg.dt) * eta_unknown + eta_implicit -
          cfg.history_weight * old_eta_in[static_cast<std::size_t>(p)];
      if (k2_line == 0.0) result.eta[static_cast<std::size_t>(p)] += cfg.explicit_weight * cfg.mean_pressure;
    }
  }

  for (int product = 0; product < 6; ++product) {
    for (int y = 0; y < grid.ny; ++y) {
      for (int line = 0; line < lines; ++line) {
        const int p = value_index(y, line, lines);
        const channel::Complex dd0 = dd(products[static_cast<std::size_t>(product)], derivatives, grid.ny, lines, y, line, 0);
        const channel::Complex dd1 = dd(products[static_cast<std::size_t>(product)], derivatives, grid.ny, lines, y, line, 1);
        const channel::Complex dd2 = dd(products[static_cast<std::size_t>(product)], derivatives, grid.ny, lines, y, line, 2);
        const double k2_line = k2[static_cast<std::size_t>(line)];
        channel::Complex rhsu(0.0, 0.0);
        channel::Complex rhsw(0.0, 0.0);
        channel::Complex d2v_expl(0.0, 0.0);
        if (product == 0) {
          rhsu = -ialfa[static_cast<std::size_t>(line)] * dd0;
          d2v_expl = ialfa[static_cast<std::size_t>(line)] * ialfa[static_cast<std::size_t>(line)] * dd1;
        } else if (product == 1) {
          d2v_expl = k2_line * dd1;
        } else if (product == 2) {
          rhsw = -ibeta[static_cast<std::size_t>(line)] * dd0;
          d2v_expl = ibeta[static_cast<std::size_t>(line)] * ibeta[static_cast<std::size_t>(line)] * dd1;
        } else if (product == 3) {
          rhsu = -dd1;
          d2v_expl = ialfa[static_cast<std::size_t>(line)] * dd2 +
                     ialfa[static_cast<std::size_t>(line)] * k2_line * dd0;
        } else if (product == 4) {
          rhsw = -dd1;
          d2v_expl = ibeta[static_cast<std::size_t>(line)] * dd2 +
                     ibeta[static_cast<std::size_t>(line)] * k2_line * dd0;
        } else {
          rhsu = -ibeta[static_cast<std::size_t>(line)] * dd0;
          rhsw = -ialfa[static_cast<std::size_t>(line)] * dd0;
          d2v_expl = 2.0 * ialfa[static_cast<std::size_t>(line)] *
                     ibeta[static_cast<std::size_t>(line)] * dd1;
        }

        result.d2v[static_cast<std::size_t>(p)] += cfg.explicit_weight * d2v_expl;
        result.old_d2v[static_cast<std::size_t>(p)] += d2v_expl;
        channel::Complex eta_expl = ibeta[static_cast<std::size_t>(line)] * rhsu -
                                    ialfa[static_cast<std::size_t>(line)] * rhsw;
        if (k2_line == 0.0) eta_expl = channel::Complex(rhsu.real(), rhsw.real());
        result.eta[static_cast<std::size_t>(p)] += cfg.explicit_weight * eta_expl;
        result.old_eta[static_cast<std::size_t>(p)] += eta_expl;
      }
    }
  }
  return result;
}

std::vector<double> make_derivatives(int ny) {
  std::vector<double> derivatives(static_cast<std::size_t>(ny) *
                                  channel::DnsNonlinearVelocityRhsStage::derivative_orders *
                                  channel::DnsNonlinearVelocityRhsStage::derivative_stencil,
                                  0.0);
  for (int y = 0; y < ny; ++y) {
    derivatives[static_cast<std::size_t>(derivative_index(y, 0, 0))] = 1.0;
    derivatives[static_cast<std::size_t>(derivative_index(y, 1, -1))] = -0.5;
    derivatives[static_cast<std::size_t>(derivative_index(y, 1, 1))] = 0.5;
    derivatives[static_cast<std::size_t>(derivative_index(y, 2, -1))] = 1.0;
    derivatives[static_cast<std::size_t>(derivative_index(y, 2, 0))] = -2.0;
    derivatives[static_cast<std::size_t>(derivative_index(y, 2, 1))] = 1.0;
    derivatives[static_cast<std::size_t>(derivative_index(y, 3, -2))] = 0.25;
    derivatives[static_cast<std::size_t>(derivative_index(y, 3, -1))] = -1.0;
    derivatives[static_cast<std::size_t>(derivative_index(y, 3, 0))] = 1.5;
    derivatives[static_cast<std::size_t>(derivative_index(y, 3, 1))] = -1.0;
    derivatives[static_cast<std::size_t>(derivative_index(y, 3, 2))] = 0.25;
  }
  return derivatives;
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1, "nonlinear velocity test runs on one MPI rank");

  const channel::DnsGrid grid{3, 6, 2, 11};
  const int lines = static_cast<int>(grid.line_count());
  channel::DnsNonlinearVelocityRhsConfig cfg;
  cfg.viscosity = 0.03;
  cfg.dt = 0.3;
  cfg.implicit_weight = 1.2;
  cfg.explicit_weight = 0.4;
  cfg.history_weight = 0.25;
  cfg.mean_pressure = channel::Complex(0.8, -0.4);

  std::vector<channel::Complex> ialfa(grid.line_count());
  std::vector<channel::Complex> ibeta(grid.line_count());
  std::vector<double> k2(grid.line_count());
  for (int line = 0; line < lines; ++line) {
    if (line == 0) {
      ialfa[static_cast<std::size_t>(line)] = channel::Complex(0.0, 0.0);
      ibeta[static_cast<std::size_t>(line)] = channel::Complex(0.0, 0.0);
      k2[static_cast<std::size_t>(line)] = 0.0;
    } else {
      ialfa[static_cast<std::size_t>(line)] = channel::Complex(0.0, 0.11 * (line + 1));
      ibeta[static_cast<std::size_t>(line)] = channel::Complex(0.0, -0.07 * (line + 2));
      k2[static_cast<std::size_t>(line)] = 0.2 + 0.03 * line;
    }
  }
  const auto derivatives = make_derivatives(grid.ny);

  std::vector<channel::Complex> u(grid.values_per_component());
  std::vector<channel::Complex> v(grid.values_per_component());
  std::vector<channel::Complex> w(grid.values_per_component());
  std::array<std::vector<channel::Complex>, 6> products;
  for (auto& product : products) product.resize(grid.values_per_component());

  for (int y = 0; y < grid.ny; ++y) {
    for (int line = 0; line < lines; ++line) {
      const int p = value_index(y, line, lines);
      u[static_cast<std::size_t>(p)] = channel::Complex(0.2 + 0.03 * y + 0.02 * line, -0.01 * (y + line));
      v[static_cast<std::size_t>(p)] = channel::Complex(-0.1 + 0.04 * y - 0.015 * line, 0.02 * y);
      w[static_cast<std::size_t>(p)] = channel::Complex(0.05 - 0.02 * y + 0.025 * line, -0.03 * line);
      for (int product = 0; product < 6; ++product) {
        products[static_cast<std::size_t>(product)][static_cast<std::size_t>(p)] =
            channel::Complex(0.07 * (product + 1) + 0.01 * y - 0.004 * line,
                             -0.03 * (product + 1) + 0.002 * y + 0.006 * line);
      }
    }
  }

  channel::DnsState state;
  state.resize(grid);
  state.copy_component_from_host(cfg.u_component, u);
  state.copy_component_from_host(cfg.v_component, v);
  state.copy_component_from_host(cfg.w_component, w);
  for (int product = 0; product < 6; ++product) {
    state.copy_component_from_host(cfg.product_components[static_cast<std::size_t>(product)],
                                   products[static_cast<std::size_t>(product)]);
  }

  channel::DnsNonlinearVelocityRhsStage stage;
  stage.prepare(state, cfg);
  stage.copy_line_wavenumbers_from_host(ialfa, ibeta, k2);
  stage.copy_y_derivatives_from_host(derivatives);

  std::vector<channel::Complex> old_eta(grid.values_per_component(), channel::Complex(0.0, 0.0));
  std::vector<channel::Complex> old_d2v(grid.values_per_component(), channel::Complex(0.0, 0.0));
  auto expected = apply_reference(grid, cfg, ialfa, ibeta, k2, derivatives, u, v, w, products, old_eta, old_d2v);
  stage.apply(state);
  auto got_eta = state.component_host(cfg.eta_rhs_component);
  auto got_d2v = state.component_host(cfg.d2v_rhs_component);
  for (std::size_t i = 0; i < got_eta.size(); ++i) {
    channel::test::require_near(got_eta[i], expected.eta[i], 1.0e-12, "first nonlinear eta RHS");
    channel::test::require_near(got_d2v[i], expected.d2v[i], 1.0e-12, "first nonlinear D2v RHS");
  }
  auto got_old_eta = stage.eta_history_host();
  auto got_old_d2v = stage.d2v_history_host();
  for (std::size_t i = 0; i < got_old_eta.size(); ++i) {
    channel::test::require_near(got_old_eta[i], expected.old_eta[i], 1.0e-12, "stored eta nonlinear history");
    channel::test::require_near(got_old_d2v[i], expected.old_d2v[i], 1.0e-12, "stored D2v nonlinear history");
  }

  auto second_expected =
      apply_reference(grid, cfg, ialfa, ibeta, k2, derivatives, u, v, w, products, expected.old_eta, expected.old_d2v);
  stage.apply(state);
  got_eta = state.component_host(cfg.eta_rhs_component);
  got_d2v = state.component_host(cfg.d2v_rhs_component);
  for (std::size_t i = 0; i < got_eta.size(); ++i) {
    channel::test::require_near(got_eta[i], second_expected.eta[i], 1.0e-12, "second nonlinear eta RHS");
    channel::test::require_near(got_d2v[i], second_expected.d2v[i], 1.0e-12, "second nonlinear D2v RHS");
  }

  channel::DnsState product_state;
  product_state.resize(grid);
  product_state.copy_component_from_host(cfg.u_component, u);
  product_state.copy_component_from_host(cfg.v_component, v);
  product_state.copy_component_from_host(cfg.w_component, w);
  channel::DnsNonlinearVelocityRhsConfig product_cfg = cfg;
  product_cfg.product_factor = 2.5;
  channel::DnsNonlinearVelocityRhsStage product_stage;
  product_stage.prepare(product_state, product_cfg);
  product_stage.build_velocity_products(product_state);
  const auto built_uv = product_state.component_host(product_cfg.product_components[static_cast<int>(channel::VelocityProduct::UV)]);
  for (std::size_t i = 0; i < built_uv.size(); ++i) {
    const channel::Complex expected_uv(2.5 * u[i].real() * v[i].real(), 0.0);
    channel::test::require_near(built_uv[i], expected_uv, 0.0, "velocity product builder uses physical real fields");
  }

  channel::DnsState wrapped_state;
  wrapped_state.resize(grid);
  wrapped_state.copy_component_from_host(cfg.u_component, u);
  wrapped_state.copy_component_from_host(cfg.v_component, v);
  wrapped_state.copy_component_from_host(cfg.w_component, w);
  for (int product = 0; product < 6; ++product) {
    wrapped_state.copy_component_from_host(cfg.product_components[static_cast<std::size_t>(product)],
                                           products[static_cast<std::size_t>(product)]);
  }
  channel::DnsVelocityStep velocity_step;
  channel::DnsVelocityStepConfig velocity_cfg;
  velocity_cfg.enable_nonlinear_velocity_rhs = true;
  velocity_cfg.nonlinear_velocity_rhs = cfg;
  velocity_step.prepare(wrapped_state, velocity_cfg);
  velocity_step.copy_nonlinear_line_wavenumbers_from_host(ialfa, ibeta, k2);
  velocity_step.copy_nonlinear_y_derivatives_from_host(derivatives);
  velocity_step.advance(wrapped_state);
  const auto wrapped_eta = wrapped_state.component_host(cfg.eta_rhs_component);
  const auto wrapped_d2v = wrapped_state.component_host(cfg.d2v_rhs_component);
  for (std::size_t i = 0; i < wrapped_eta.size(); ++i) {
    channel::test::require_near(wrapped_eta[i], expected.eta[i], 1.0e-12, "velocity step nonlinear eta RHS");
    channel::test::require_near(wrapped_d2v[i], expected.d2v[i], 1.0e-12, "velocity step nonlinear D2v RHS");
  }
  channel::test::require(velocity_step.has_nonlinear_velocity_rhs(), "velocity step owns nonlinear RHS stage");

  if (runtime.rank() == 0) std::cout << "nonlinear velocity RHS test PASSED\n";
  return 0;
}
