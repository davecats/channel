#include "channel/io.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

namespace {

std::filesystem::path fixture_path(const char* name) {
  const std::vector<std::filesystem::path> candidates = {
      std::filesystem::path("tests/data") / name,
      std::filesystem::path("../tests/data") / name,
      std::filesystem::path("../../tests/data") / name,
  };
  for (const auto& candidate : candidates) {
    if (std::filesystem::exists(candidate)) return candidate;
  }
  throw std::runtime_error(std::string("could not find fixture ") + name);
}

std::size_t fortran_index(const channel::LegacyRestartMetadata& header,
                          int y_file,
                          int z_file,
                          int x_file,
                          int component) {
  const auto y_extent = static_cast<std::size_t>(header.ny + 3);
  const auto z_extent = static_cast<std::size_t>(2 * header.nz + 1);
  const auto x_extent = static_cast<std::size_t>(header.nx + 1);
  return static_cast<std::size_t>(y_file) +
         y_extent * (static_cast<std::size_t>(z_file) +
                     z_extent * (static_cast<std::size_t>(x_file) +
                                 x_extent * static_cast<std::size_t>(component)));
}

std::vector<channel::Complex> raw_payload(const std::filesystem::path& path) {
  constexpr std::streamsize header_bytes = 3 * static_cast<std::streamsize>(sizeof(std::int32_t)) +
                                           7 * static_cast<std::streamsize>(sizeof(double));
  std::ifstream in(path, std::ios::binary);
  channel::test::require(static_cast<bool>(in), "open restart fixture");
  in.seekg(0, std::ios::end);
  const auto file_size = static_cast<std::streamsize>(in.tellg());
  channel::test::require(file_size > header_bytes, "restart fixture has payload");
  std::vector<channel::Complex> payload(static_cast<std::size_t>(file_size - header_bytes) /
                                        sizeof(channel::Complex));
  in.seekg(header_bytes, std::ios::beg);
  in.read(reinterpret_cast<char*>(payload.data()),
          static_cast<std::streamsize>(payload.size() * sizeof(channel::Complex)));
  channel::test::require(static_cast<bool>(in), "read restart fixture payload");
  return payload;
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 1, "restart IO test runs on one MPI rank");

  const auto path = fixture_path("start_field.out");
  channel::DnsState state;
  const auto header = channel::read_legacy_fortran_state_for_tests(path.string(), state);
  channel::test::require(header.nx == 5 && header.ny == 16 && header.nz == 8, "legacy restart metadata dimensions");
  channel::test::require(header.components == 3, "legacy restart component count");
  channel::test::require(std::abs(header.alfa0 - 0.5) < 1.0e-15, "legacy restart alfa0");
  channel::test::require(std::abs(header.beta0 - 1.0) < 1.0e-15, "legacy restart beta0");
  channel::test::require(std::abs(header.ni - 1.0 / 12431.0) < 1.0e-15, "legacy restart viscosity");
  channel::test::require(state.grid().nx == 5 && state.grid().ny == 16 && state.grid().nz == 8,
                         "legacy restart state dimensions");
  channel::test::require(state.grid().components == 3, "legacy restart state components");

  const auto payload = raw_payload(path);
  const auto u = state.component_host(0);
  const auto v = state.component_host(1);
  const int nx = header.nx;
  const int nz = header.nz;
  const int lines = nx * nz;
  const auto check_value = [&](int y, int iz, int ix, int component, const std::vector<channel::Complex>& values) {
    const int signed_iz = iz <= nz / 2 ? iz : iz - nz;
    const int file_z = signed_iz + nz;
    const auto want = payload[fortran_index(header, y + 1, file_z, ix, component)];
    const auto got = values[static_cast<std::size_t>(y * lines + iz * nx + ix)];
    channel::test::require_near(got, want, 0.0, "legacy restart Fortran-order payload mapping");
  };
  check_value(0, 0, 0, 0, u);
  check_value(5, 3, 2, 0, u);
  check_value(15, 7, 4, 1, v);

  channel::DnsState full_state;
  const auto full_header = channel::read_legacy_fortran_full_spectrum_with_ghosts_for_tests(path.string(), full_state);
  channel::test::require(full_header.nx == header.nx && full_header.ny == header.ny && full_header.nz == header.nz,
                         "full-spectrum restart metadata dimensions");
  channel::test::require(full_state.grid().nx == header.nx + 1 && full_state.grid().ny == header.ny + 3 &&
                             full_state.grid().nz == 2 * header.nz + 1,
                         "full-spectrum restart state dimensions");
  const auto full_u = full_state.component_host(0);
  const int full_lines = full_state.grid().nx * full_state.grid().nz;
  const auto check_full_value = [&](int fortran_y, int signed_iz, int ix, int component,
                                    const std::vector<channel::Complex>& values) {
    const int canonical_z = signed_iz >= 0 ? signed_iz : signed_iz + 2 * header.nz + 1;
    const auto want = payload[fortran_index(header, fortran_y + 1, signed_iz + header.nz, ix, component)];
    const auto got = values[static_cast<std::size_t>((fortran_y + 1) * full_lines +
                                                     canonical_z * full_state.grid().nx + ix)];
    channel::test::require_near(got, want, 0.0, "full-spectrum legacy restart payload mapping");
  };
  check_full_value(-1, -header.nz, header.nx, 0, full_u);
  check_full_value(0, 0, 0, 0, full_u);
  check_full_value(header.ny + 1, header.nz, 3, 0, full_u);

  if (runtime.rank() == 0) std::cout << "Legacy restart IO test PASSED\n";
  return 0;
}
