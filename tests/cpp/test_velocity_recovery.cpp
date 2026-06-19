#include "channel/runtime.hpp"
#include "channel/velocity_recovery.hpp"

#include "test_common.hpp"

#include <iostream>
#include <vector>

namespace {

int value_index(int y, int line, int lines) {
  return y * lines + line;
}

int derivative_index(int y, int order, int offset) {
  return (y * channel::DnsVelocityRecoveryStage::derivative_orders + order) *
             channel::DnsVelocityRecoveryStage::derivative_stencil +
         (offset + 2);
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1, "velocity recovery test runs on one MPI rank");

  const channel::DnsGrid grid{3, 5, 2, 5};
  const int lines = static_cast<int>(grid.line_count());
  channel::DnsState state;
  state.resize(grid);

  std::vector<channel::Complex> u(grid.values_per_component());
  std::vector<channel::Complex> v(grid.values_per_component());
  std::vector<channel::Complex> w(grid.values_per_component());
  std::vector<channel::Complex> eta(grid.values_per_component());
  for (int y = 0; y < grid.ny; ++y) {
    for (int line = 0; line < lines; ++line) {
      const int p = value_index(y, line, lines);
      u[static_cast<std::size_t>(p)] = channel::Complex(10.0 + y, -line);
      v[static_cast<std::size_t>(p)] = channel::Complex(0.2 + 0.03 * y + 0.01 * line,
                                                        -0.4 + 0.02 * y - 0.005 * line);
      w[static_cast<std::size_t>(p)] = channel::Complex(-7.0 - y, line);
      eta[static_cast<std::size_t>(p)] = channel::Complex(0.1 + 0.04 * y - 0.02 * line,
                                                          -0.3 + 0.01 * y + 0.03 * line);
    }
  }
  state.copy_component_from_host(0, u);
  state.copy_component_from_host(1, v);
  state.copy_component_from_host(2, w);
  state.copy_component_from_host(3, eta);

  std::vector<channel::Complex> ialfa(grid.line_count());
  std::vector<channel::Complex> ibeta(grid.line_count());
  std::vector<double> k2(grid.line_count());
  for (int line = 0; line < lines; ++line) {
    if (line == 0) {
      ialfa[static_cast<std::size_t>(line)] = channel::Complex(0.0, 0.0);
      ibeta[static_cast<std::size_t>(line)] = channel::Complex(0.0, 0.0);
      k2[static_cast<std::size_t>(line)] = 0.0;
    } else {
      ialfa[static_cast<std::size_t>(line)] = channel::Complex(0.0, 0.12 * (line + 1));
      ibeta[static_cast<std::size_t>(line)] = channel::Complex(0.0, -0.08 * (line + 2));
      k2[static_cast<std::size_t>(line)] = 0.3 + 0.02 * line;
    }
  }

  std::vector<double> derivatives(static_cast<std::size_t>(grid.ny) *
                                  channel::DnsVelocityRecoveryStage::derivative_orders *
                                  channel::DnsVelocityRecoveryStage::derivative_stencil,
                                  0.0);
  for (int y = 0; y < grid.ny; ++y) {
    derivatives[static_cast<std::size_t>(derivative_index(y, 1, -1))] = -0.5;
    derivatives[static_cast<std::size_t>(derivative_index(y, 1, 1))] = 0.5;
  }

  channel::DnsVelocityRecoveryStage recovery;
  channel::DnsVelocityRecoveryConfig cfg;
  cfg.u_component = 0;
  cfg.v_component = 1;
  cfg.w_component = 2;
  cfg.eta_component = 3;
  cfg.dvdy_component = 4;
  recovery.prepare(state, cfg);
  recovery.copy_line_wavenumbers_from_host(ialfa, ibeta, k2);
  recovery.copy_y_derivatives_from_host(derivatives);
  recovery.apply(state);

  const auto got_u = state.component_host(0);
  const auto got_w = state.component_host(2);
  const auto got_dvdy = state.component_host(4);
  for (int y = 0; y < grid.ny; ++y) {
    for (int line = 0; line < lines; ++line) {
      const int p = value_index(y, line, lines);
      channel::Complex dvdy(0.0, 0.0);
      for (int offset = -2; offset <= 2; ++offset) {
        const int yy = y + offset;
        if (yy < 0 || yy >= grid.ny) continue;
        dvdy += derivatives[static_cast<std::size_t>(derivative_index(y, 1, offset))] *
                v[static_cast<std::size_t>(value_index(yy, line, lines))];
      }
      channel::test::require_near(got_dvdy[static_cast<std::size_t>(p)], dvdy, 1.0e-13,
                                  "recovered dvdy");
      if (line == 0) {
        channel::test::require_near(got_u[static_cast<std::size_t>(p)], u[static_cast<std::size_t>(p)], 0.0,
                                    "zero-mode u preserved");
        channel::test::require_near(got_w[static_cast<std::size_t>(p)], w[static_cast<std::size_t>(p)], 0.0,
                                    "zero-mode w preserved");
      } else {
        const auto want_u = (ialfa[static_cast<std::size_t>(line)] * dvdy -
                             ibeta[static_cast<std::size_t>(line)] * eta[static_cast<std::size_t>(p)]) /
                            k2[static_cast<std::size_t>(line)];
        const auto want_w = (ibeta[static_cast<std::size_t>(line)] * dvdy +
                             ialfa[static_cast<std::size_t>(line)] * eta[static_cast<std::size_t>(p)]) /
                            k2[static_cast<std::size_t>(line)];
        channel::test::require_near(got_u[static_cast<std::size_t>(p)], want_u, 1.0e-13, "recovered u");
        channel::test::require_near(got_w[static_cast<std::size_t>(p)], want_w, 1.0e-13, "recovered w");
      }
    }
  }

  if (runtime.rank() == 0) std::cout << "velocity recovery test PASSED\n";
  return 0;
}
