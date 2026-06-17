#include "channel/nonlinear.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <iostream>
#include <vector>

namespace {

int value_index(int y, int line, int lines) {
  return y * lines + line;
}

void fill_velocity(channel::DnsState& state, double u_value, double v_value, double w_value) {
  const auto& grid = state.grid();
  const int lines = static_cast<int>(grid.line_count());
  std::vector<channel::Complex> u(grid.values_per_component());
  std::vector<channel::Complex> v(grid.values_per_component());
  std::vector<channel::Complex> w(grid.values_per_component());
  for (int y = 0; y < grid.ny; ++y) {
    for (int line = 0; line < lines; ++line) {
      const int p = value_index(y, line, lines);
      u[static_cast<std::size_t>(p)] = channel::Complex(u_value, 0.1 * y);
      v[static_cast<std::size_t>(p)] = channel::Complex(v_value, -0.2 * line);
      w[static_cast<std::size_t>(p)] = channel::Complex(w_value, 0.3 * (y + line));
    }
  }
  state.copy_component_from_host(0, u);
  state.copy_component_from_host(1, v);
  state.copy_component_from_host(2, w);
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1, "nonlinear product transform test runs on one MPI rank");

  const channel::DnsGrid grid{2, 4, 1, 11};
  const int lines = static_cast<int>(grid.line_count());

  channel::DnsState physical_state;
  physical_state.resize(grid);
  fill_velocity(physical_state, 2.0, 3.0, 5.0);

  channel::DnsNonlinearProductTransformConfig physical_cfg;
  physical_cfg.enable_velocity_inverse_ffts = false;
  physical_cfg.enable_product_forward_ffts = false;
  physical_cfg.product_factor = 0.25;
  channel::DnsNonlinearProductTransformStage physical_stage;
  physical_stage.prepare(physical_state, physical_cfg);
  physical_stage.apply(physical_state);
  channel::test::require(physical_stage.velocity_fft_count() == 0, "disabled velocity FFT count");
  channel::test::require(physical_stage.product_fft_count() == 0, "disabled product FFT count");

  const auto uu = physical_state.component_host(5);
  const auto vv = physical_state.component_host(6);
  const auto ww = physical_state.component_host(7);
  const auto uv = physical_state.component_host(8);
  const auto vw = physical_state.component_host(9);
  const auto uw = physical_state.component_host(10);
  for (std::size_t i = 0; i < uu.size(); ++i) {
    channel::test::require_near(uu[i], channel::Complex(1.0, 0.0), 0.0, "uu physical product");
    channel::test::require_near(vv[i], channel::Complex(2.25, 0.0), 0.0, "vv physical product");
    channel::test::require_near(ww[i], channel::Complex(6.25, 0.0), 0.0, "ww physical product");
    channel::test::require_near(uv[i], channel::Complex(1.5, 0.0), 0.0, "uv physical product");
    channel::test::require_near(vw[i], channel::Complex(3.75, 0.0), 0.0, "vw physical product");
    channel::test::require_near(uw[i], channel::Complex(2.5, 0.0), 0.0, "uw physical product");
  }

  channel::DnsState spectral_state;
  spectral_state.resize(grid);
  fill_velocity(spectral_state, 2.0, 3.0, 5.0);
  channel::DnsNonlinearProductTransformConfig spectral_cfg;
  spectral_cfg.enable_velocity_inverse_ffts = false;
  spectral_cfg.enable_product_forward_ffts = true;
  spectral_cfg.product_forward_normalization = channel::FftNormalization::None;
  spectral_cfg.product_factor = 0.25;
  channel::DnsNonlinearProductTransformStage spectral_stage;
  spectral_stage.prepare(spectral_state, spectral_cfg);
  spectral_stage.apply(spectral_state);
  channel::test::require(spectral_stage.product_fft_count() == 6, "product FFT count");

  const auto uu_spectrum = spectral_state.component_host(5);
  for (int line = 0; line < lines; ++line) {
    channel::test::require_near(
        uu_spectrum[static_cast<std::size_t>(value_index(0, line, lines))],
        channel::Complex(static_cast<double>(grid.ny), 0.0),
        1.0e-12,
        "constant product zero y-mode");
    for (int y = 1; y < grid.ny; ++y) {
      channel::test::require_near(
          uu_spectrum[static_cast<std::size_t>(value_index(y, line, lines))],
          channel::Complex(0.0, 0.0),
          1.0e-12,
          "constant product nonzero y-mode");
    }
  }

  if (runtime.rank() == 0) std::cout << "nonlinear product transform test PASSED\n";
  return 0;
}
