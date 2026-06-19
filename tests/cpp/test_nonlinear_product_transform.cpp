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

channel::Complex spectral_value(int component, int y, int global_x, int z_storage, int storage_z) {
  const int retained_z = (storage_z - 1) / 2;
  const int signed_z = z_storage <= retained_z ? z_storage : z_storage - storage_z;
  const double real = 0.25 * component + 0.125 * y + 0.5 * global_x - 0.2 * signed_z;
  const double imag = -0.1 * component + 0.0625 * y - 0.3 * global_x + 0.15 * signed_z;
  return {real, imag};
}

void fill_spectral_velocity(channel::DnsState& state, int global_x_first, int storage_z) {
  const auto& grid = state.grid();
  const int lines = static_cast<int>(grid.line_count());
  for (int component = 0; component < 3; ++component) {
    std::vector<channel::Complex> values(grid.values_per_component());
    for (int y = 0; y < grid.ny; ++y) {
      for (int z = 0; z < grid.nz; ++z) {
        for (int x = 0; x < grid.nx; ++x) {
          const int p = value_index(y, z * grid.nx + x, lines);
          values[static_cast<std::size_t>(p)] =
              spectral_value(component, y, global_x_first + x, z, storage_z);
        }
      }
    }
    state.copy_component_from_host(component, values);
  }
}

void exercise_single_rank_transform() {
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
  for (int y = 0; y < grid.ny; ++y) {
    for (int line = 0; line < lines; ++line) {
      const auto expected = line == 0 ? channel::Complex(static_cast<double>(grid.nx * grid.nz), 0.0)
                                      : channel::Complex(0.0, 0.0);
      channel::test::require_near(
          uu_spectrum[static_cast<std::size_t>(value_index(y, line, lines))],
          expected,
          1.0e-12,
          "constant product xz slab spectrum");
    }
  }
}

void exercise_distributed_transform(const channel::Runtime& runtime) {
  constexpr int global_x = 4;
  constexpr int local_x = global_x / 2;
  constexpr int storage_z = 5;
  constexpr int padded_z = 6;
  constexpr int physical_x = 8;
  constexpr int ny = 3;
  constexpr int components = 11;
  const int global_x_first = runtime.rank() * local_x;

  channel::DnsState local_state;
  local_state.resize({local_x, ny, storage_z, components});
  local_state.fill(channel::Complex(0.0, 0.0));
  fill_spectral_velocity(local_state, global_x_first, storage_z);

  channel::DnsNonlinearProductTransformConfig distributed_cfg;
  distributed_cfg.dealiased_physical_x = physical_x;
  distributed_cfg.dealiased_physical_z = padded_z;
  distributed_cfg.product_factor = 0.125;
  distributed_cfg.distributed_dealiased_fft = true;
  distributed_cfg.distributed_spectral_x_total = global_x;
  distributed_cfg.distributed_spectral_x_first = global_x_first;
  distributed_cfg.distributed_npxz = runtime.size();
  distributed_cfg.distributed_ipxz = runtime.rank();
  distributed_cfg.distributed_comm_x = channel::world_comm();
  channel::DnsNonlinearProductTransformStage distributed_stage;
  distributed_stage.prepare(local_state, distributed_cfg);
  distributed_stage.apply(local_state);

  channel::DnsState expected_state;
  expected_state.resize({global_x, ny, storage_z, components});
  expected_state.fill(channel::Complex(0.0, 0.0));
  fill_spectral_velocity(expected_state, 0, storage_z);

  channel::DnsNonlinearProductTransformConfig expected_cfg;
  expected_cfg.dealiased_physical_x = physical_x;
  expected_cfg.dealiased_physical_z = padded_z;
  expected_cfg.product_factor = distributed_cfg.product_factor;
  channel::DnsNonlinearProductTransformStage expected_stage;
  expected_stage.prepare(expected_state, expected_cfg);
  expected_stage.apply(expected_state);

  const int local_lines = static_cast<int>(local_state.line_count());
  const int expected_lines = static_cast<int>(expected_state.line_count());
  for (int component = 5; component <= 10; ++component) {
    const auto got = local_state.component_host(component);
    const auto want = expected_state.component_host(component);
    for (int y = 0; y < ny; ++y) {
      for (int z = 0; z < storage_z; ++z) {
        for (int x = 0; x < local_x; ++x) {
          const int local_p = value_index(y, z * local_x + x, local_lines);
          const int expected_p = value_index(y, z * global_x + global_x_first + x, expected_lines);
          channel::test::require_near(got[static_cast<std::size_t>(local_p)],
                                      want[static_cast<std::size_t>(expected_p)],
                                      1.0e-10,
                                      "distributed dealiased product spectrum");
        }
      }
    }
  }
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  if (runtime.size() == 1) {
    exercise_single_rank_transform();
  } else {
    channel::test::require(runtime.size() == 2, "distributed nonlinear product transform test expects two ranks");
    exercise_distributed_transform(runtime);
  }

  if (runtime.rank() == 0) std::cout << "nonlinear product transform test PASSED\n";
  return 0;
}
