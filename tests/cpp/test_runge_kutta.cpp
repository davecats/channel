#include "channel/runge_kutta.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

int derivative_index(int y, int order, int offset) {
  return (y * channel::DnsNonlinearVelocityRhsStage::derivative_orders + order) *
             channel::DnsNonlinearVelocityRhsStage::derivative_stencil +
         (offset + 2);
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1, "Runge-Kutta timestepper test runs on one MPI rank");

  const auto weights = channel::channel_rk3_weights();
  channel::test::require(std::abs(weights[0][0] - 120.0 / 32.0) < 1.0e-15, "RK stage 1 implicit weight");
  channel::test::require(std::abs(weights[1][1] - 50.0 / 8.0) < 1.0e-15, "RK stage 2 explicit weight");
  channel::test::require(std::abs(weights[2][2] - 50.0 / 20.0) < 1.0e-15, "RK stage 3 history weight");

  const channel::DnsGrid grid{1, 4, 1, 11};
  channel::DnsState state;
  state.resize(grid);
  state.fill(channel::Complex(0.0, 0.0));

  channel::DnsNonlinearVelocityRhsConfig nonlinear_cfg;
  nonlinear_cfg.u_component = 0;
  nonlinear_cfg.v_component = 1;
  nonlinear_cfg.w_component = 2;
  nonlinear_cfg.eta_rhs_component = 3;
  nonlinear_cfg.d2v_rhs_component = 4;
  nonlinear_cfg.product_components = {5, 6, 7, 8, 9, 10};
  nonlinear_cfg.viscosity = 0.0;

  std::vector<channel::Complex> uu(grid.values_per_component());
  for (int y = 0; y < grid.ny; ++y) {
    uu[static_cast<std::size_t>(y)] = channel::Complex(0.3 + 0.02 * y, -0.1 + 0.01 * y);
  }
  state.copy_component_from_host(nonlinear_cfg.product_components[0], uu);

  std::vector<channel::Complex> ialfa = {channel::Complex(0.0, 0.4)};
  std::vector<channel::Complex> ibeta = {channel::Complex(0.0, -0.25)};
  std::vector<double> k2 = {0.41};
  std::vector<double> derivatives(static_cast<std::size_t>(grid.ny) *
                                  channel::DnsNonlinearVelocityRhsStage::derivative_orders *
                                  channel::DnsNonlinearVelocityRhsStage::derivative_stencil,
                                  0.0);
  for (int y = 0; y < grid.ny; ++y) {
    derivatives[static_cast<std::size_t>(derivative_index(y, 0, 0))] = 1.0;
    derivatives[static_cast<std::size_t>(derivative_index(y, 1, 0))] = 1.0;
  }

  const double dt = 0.125;
  channel::DnsRungeKuttaTimestepperConfig cfg;
  cfg.dt = dt;
  cfg.time = 2.0;
  cfg.enable_nonlinear_velocity_rhs = true;
  cfg.nonlinear_velocity_rhs = nonlinear_cfg;
  cfg.stages.resize(weights.size());
  for (std::size_t stage = 0; stage < weights.size(); ++stage) {
    cfg.stages[stage].weights = weights[stage];
  }

  channel::DnsRungeKuttaTimestepper stepper;
  stepper.prepare(state, cfg);
  stepper.copy_nonlinear_line_wavenumbers_from_host(ialfa, ibeta, k2);
  stepper.copy_nonlinear_y_derivatives_from_host(derivatives);
  stepper.advance_one_step(state);

  const double expected_time =
      2.0 + 2.0 / weights[0][0] * dt + 2.0 / weights[1][0] * dt + 2.0 / weights[2][0] * dt;
  channel::test::require(std::abs(stepper.time() - expected_time) < 1.0e-15, "RK timestepper advances substep time");
  channel::test::require(stepper.stage_count() == 3, "RK timestepper owns three stages");
  channel::test::require(stepper.has_nonlinear_velocity_rhs(), "RK timestepper owns nonlinear RHS");

  const auto got_eta = state.component_host(nonlinear_cfg.eta_rhs_component);
  const auto got_d2v = state.component_host(nonlinear_cfg.d2v_rhs_component);
  for (int y = 0; y < grid.ny; ++y) {
    const auto product = uu[static_cast<std::size_t>(y)];
    const auto rhsu = -ialfa[0] * product;
    const auto eta_expl = ibeta[0] * rhsu;
    const auto d2v_expl = ialfa[0] * ialfa[0] * product;
    const double final_factor = weights[2][1] - weights[2][2];
    channel::test::require_near(got_eta[static_cast<std::size_t>(y)], final_factor * eta_expl, 1.0e-13,
                                "RK final eta RHS uses shared oldrhs history");
    channel::test::require_near(got_d2v[static_cast<std::size_t>(y)], final_factor * d2v_expl, 1.0e-13,
                                "RK final D2v RHS uses shared oldrhs history");
  }

  channel::DnsState closed_state;
  const channel::DnsGrid closed_grid{1, 3, 2, 11};
  closed_state.resize(closed_grid);
  closed_state.fill(channel::Complex(0.0, 0.0));
  std::vector<channel::Complex> closed_u(closed_grid.values_per_component());
  std::vector<channel::Complex> closed_v(closed_grid.values_per_component());
  std::vector<channel::Complex> closed_w(closed_grid.values_per_component());
  for (int y = 0; y < closed_grid.ny; ++y) {
    for (int line = 0; line < static_cast<int>(closed_grid.line_count()); ++line) {
      const int p = y * static_cast<int>(closed_grid.line_count()) + line;
      closed_u[static_cast<std::size_t>(p)] = channel::Complex(0.2 + 0.03 * y + 0.01 * line, 0.0);
      closed_v[static_cast<std::size_t>(p)] = channel::Complex(-0.1 + 0.02 * y - 0.01 * line, 0.0);
      closed_w[static_cast<std::size_t>(p)] = channel::Complex(0.05 - 0.01 * y + 0.02 * line, 0.0);
    }
  }
  closed_state.copy_component_from_host(0, closed_u);
  closed_state.copy_component_from_host(1, closed_v);
  closed_state.copy_component_from_host(2, closed_w);

  std::vector<channel::Complex> closed_ialfa = {channel::Complex(0.0, 0.0), channel::Complex(0.0, 0.3)};
  std::vector<channel::Complex> closed_ibeta = {channel::Complex(0.0, 0.0), channel::Complex(0.0, -0.2)};
  std::vector<double> closed_k2 = {0.0, 0.13};
  std::vector<double> closed_derivatives(static_cast<std::size_t>(closed_grid.ny) *
                                         channel::DnsNonlinearVelocityRhsStage::derivative_orders *
                                         channel::DnsNonlinearVelocityRhsStage::derivative_stencil,
                                         0.0);
  for (int y = 0; y < closed_grid.ny; ++y) {
    closed_derivatives[static_cast<std::size_t>(derivative_index(y, 0, 0))] = 1.0;
    closed_derivatives[static_cast<std::size_t>(derivative_index(y, 1, 0))] = 1.0;
  }

  channel::DnsRungeKuttaTimestepperConfig closed_cfg;
  closed_cfg.dt = 0.5;
  closed_cfg.time = 0.0;
  closed_cfg.enable_nonlinear_product_transform = true;
  closed_cfg.nonlinear_product_transform.u_component = 0;
  closed_cfg.nonlinear_product_transform.v_component = 1;
  closed_cfg.nonlinear_product_transform.w_component = 2;
  closed_cfg.nonlinear_product_transform.product_components = {5, 6, 7, 8, 9, 10};
  closed_cfg.nonlinear_product_transform.enable_velocity_inverse_ffts = false;
  closed_cfg.nonlinear_product_transform.enable_product_forward_ffts = false;
  closed_cfg.enable_nonlinear_velocity_rhs = true;
  closed_cfg.nonlinear_velocity_rhs.u_component = 0;
  closed_cfg.nonlinear_velocity_rhs.v_component = 1;
  closed_cfg.nonlinear_velocity_rhs.w_component = 2;
  closed_cfg.nonlinear_velocity_rhs.eta_rhs_component = 0;
  closed_cfg.nonlinear_velocity_rhs.d2v_rhs_component = 1;
  closed_cfg.nonlinear_velocity_rhs.product_components = {5, 6, 7, 8, 9, 10};
  closed_cfg.enable_velocity_recovery = true;
  closed_cfg.velocity_recovery.u_component = 0;
  closed_cfg.velocity_recovery.v_component = 1;
  closed_cfg.velocity_recovery.w_component = 2;
  closed_cfg.velocity_recovery.eta_component = 0;
  closed_cfg.velocity_recovery.dvdy_component = 2;
  closed_cfg.stages.resize(1);
  closed_cfg.stages[0].weights = {1.0, 0.0, 0.0};
  closed_cfg.stages[0].linear.enable_implicit_yline = true;
  closed_cfg.stages[0].linear.implicit_yline.components = {0, 1};
  closed_cfg.stages[0].linear.implicit_yline.npy = 1;
  closed_cfg.stages[0].linear.implicit_yline.ipy = 0;
  closed_cfg.stages[0].linear.implicit_yline.comm_y = channel::world_comm();

  channel::DnsRungeKuttaTimestepper closed_stepper;
  closed_stepper.prepare(closed_state, closed_cfg);
  closed_stepper.copy_nonlinear_line_wavenumbers_from_host(closed_ialfa, closed_ibeta, closed_k2);
  closed_stepper.copy_nonlinear_y_derivatives_from_host(closed_derivatives);
  closed_stepper.copy_recovery_line_wavenumbers_from_host(closed_ialfa, closed_ibeta, closed_k2);
  closed_stepper.copy_recovery_y_derivatives_from_host(closed_derivatives);
  std::vector<channel::Complex> zero(closed_state.values_per_component(), channel::Complex(0.0, 0.0));
  std::vector<channel::Complex> one(closed_state.values_per_component(), channel::Complex(1.0, 0.0));
  closed_stepper.copy_yline_coefficients_from_host(0, 0, zero, zero, one, zero, zero);
  closed_stepper.copy_yline_coefficients_from_host(0, 1, zero, zero, one, zero, zero);
  closed_stepper.advance_one_step(closed_state);
  channel::test::require(closed_stepper.has_nonlinear_product_transform(), "RK timestepper owns product transform");
  channel::test::require(closed_stepper.has_velocity_recovery(), "RK timestepper owns velocity recovery");

  const auto closed_got_u = closed_state.component_host(0);
  const auto closed_got_w = closed_state.component_host(2);
  for (int y = 0; y < closed_grid.ny; ++y) {
    const int p = y * static_cast<int>(closed_grid.line_count()) + 1;
    const auto eta_rhs = (closed_ibeta[1] * closed_u[static_cast<std::size_t>(p)] -
                          closed_ialfa[1] * closed_w[static_cast<std::size_t>(p)]) /
                         closed_cfg.dt;
    const auto v_rhs = (-closed_k2[1] * closed_v[static_cast<std::size_t>(p)]) / closed_cfg.dt;
    const auto want_u = (closed_ialfa[1] * v_rhs - closed_ibeta[1] * eta_rhs) / closed_k2[1];
    const auto want_w = (closed_ibeta[1] * v_rhs + closed_ialfa[1] * eta_rhs) / closed_k2[1];
    channel::test::require_near(closed_got_u[static_cast<std::size_t>(p)], want_u, 1.0e-13,
                                "closed RK velocity pipeline recovers u");
    channel::test::require_near(closed_got_w[static_cast<std::size_t>(p)], want_w, 1.0e-13,
                                "closed RK velocity pipeline recovers w");
  }

  if (runtime.rank() == 0) std::cout << "Runge-Kutta timestepper test PASSED\n";
  return 0;
}
