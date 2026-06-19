#include "channel/io.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct NativeHeader {
  std::uint32_t magic = 0x43484E4C; // "CHNL"
  std::uint32_t version = 1;
  std::uint32_t nx = 0;
  std::uint32_t ny = 0;
  std::uint32_t nz = 0;
  std::uint32_t components = 0;
};

constexpr std::streamsize legacy_header_bytes = 3 * static_cast<std::streamsize>(sizeof(std::int32_t)) +
                                              7 * static_cast<std::streamsize>(sizeof(double));

[[nodiscard]] std::size_t total_payload_bytes(const channel::LegacyRestartMetadata& h) {
  return static_cast<std::size_t>(h.ny + 3) * static_cast<std::size_t>(2 * h.nz + 1) *
         static_cast<std::size_t>(h.nx + 1) * sizeof(channel::Complex);
}

[[nodiscard]] std::size_t legacy_fortran_index(const channel::LegacyRestartMetadata& h,
                                               int y_file,
                                               int z_file,
                                               int x_file,
                                               int component) {
  const auto y_extent = static_cast<std::size_t>(h.ny + 3);
  const auto z_extent = static_cast<std::size_t>(2 * h.nz + 1);
  const auto x_extent = static_cast<std::size_t>(h.nx + 1);
  return static_cast<std::size_t>(y_file) +
         y_extent * (static_cast<std::size_t>(z_file) +
                     z_extent * (static_cast<std::size_t>(x_file) +
                                 x_extent * static_cast<std::size_t>(component)));
}

[[nodiscard]] int canonical_signed_z(int z, int nz) {
  return z <= nz ? z : z - (2 * nz + 1);
}

} // namespace

namespace channel {

void write_state_to_file(const DnsState& state, std::string_view path) {
  if (state.grid().components < 1) {
    throw std::runtime_error("write_state_to_file requires a resized DnsState");
  }

  const auto& grid = state.grid();
  NativeHeader header;
  header.nx = static_cast<std::uint32_t>(grid.nx);
  header.ny = static_cast<std::uint32_t>(grid.ny);
  header.nz = static_cast<std::uint32_t>(grid.nz);
  header.components = static_cast<std::uint32_t>(grid.components);

  std::ofstream out(std::string(path), std::ios::binary);
  if (!out) {
    throw std::runtime_error("Failed to open native checkpoint for writing: " + std::string(path));
  }
  out.write(reinterpret_cast<const char*>(&header), static_cast<std::streamsize>(sizeof(header)));
  if (!out) {
    throw std::runtime_error("Failed writing native checkpoint header: " + std::string(path));
  }

  for (int component = 0; component < grid.components; ++component) {
    const auto values = state.component_host(component);
    out.write(reinterpret_cast<const char*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(Complex)));
    if (!out) {
      throw std::runtime_error("Failed writing native checkpoint body: " + std::string(path));
    }
  }
}

struct LegacyRestartPayload {
  LegacyRestartMetadata header;
  std::vector<Complex> values;
};

LegacyRestartPayload read_legacy_payload(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("Failed to open legacy restart fixture: " + path);
  }

  LegacyRestartPayload payload;
  auto& header = payload.header;
  in.read(reinterpret_cast<char*>(&header.nx), sizeof(header.nx));
  in.read(reinterpret_cast<char*>(&header.ny), sizeof(header.ny));
  in.read(reinterpret_cast<char*>(&header.nz), sizeof(header.nz));
  in.read(reinterpret_cast<char*>(&header.alfa0), sizeof(header.alfa0));
  in.read(reinterpret_cast<char*>(&header.beta0), sizeof(header.beta0));
  in.read(reinterpret_cast<char*>(&header.ni), sizeof(header.ni));
  in.read(reinterpret_cast<char*>(&header.stretching), sizeof(header.stretching));
  in.read(reinterpret_cast<char*>(&header.ymin), sizeof(header.ymin));
  in.read(reinterpret_cast<char*>(&header.ymax), sizeof(header.ymax));
  in.read(reinterpret_cast<char*>(&header.time), sizeof(header.time));
  if (!in) {
    throw std::runtime_error("Failed reading legacy restart metadata: " + path);
  }

  const std::size_t payload_bytes_per_component = total_payload_bytes(header);
  in.seekg(0, std::ios::end);
  const std::streampos end = in.tellg();
  if (end == std::streampos(-1)) {
    throw std::runtime_error("Unable to determine legacy restart fixture size: " + path);
  }
  const std::streamsize file_size = static_cast<std::streamsize>(end);
  if (file_size < legacy_header_bytes) {
    throw std::runtime_error("Restart fixture is smaller than expected header: " + path);
  }
  const std::size_t data_bytes =
      static_cast<std::size_t>(file_size) - static_cast<std::size_t>(legacy_header_bytes);
  if ((data_bytes % payload_bytes_per_component) != 0) {
    throw std::runtime_error("Legacy fixture payload size mismatch for: " + path);
  }
  const auto legacy_components = data_bytes / payload_bytes_per_component;
  header.components = static_cast<int>(legacy_components);
  if ((payload_bytes_per_component % sizeof(Complex)) != 0) {
    throw std::runtime_error("Restart fixture payload is not an integer Complex multiple: " + path);
  }
  const std::size_t payload_values_per_component =
      payload_bytes_per_component / sizeof(Complex);
  if (legacy_components == 0 || legacy_components > 1000000) {
    throw std::runtime_error("Legacy fixture component count appears invalid: " + path);
  }

  const int expected_nx = static_cast<int>(header.nx);
  const int expected_ny = static_cast<int>(header.ny);
  const int expected_nz = static_cast<int>(header.nz);

  if (expected_nx < 1 || expected_ny < 1 || expected_nz < 1) {
    throw std::runtime_error("Legacy fixture has invalid domain metadata: " + path);
  }

  payload.values.resize(data_bytes / sizeof(Complex));
  in.seekg(legacy_header_bytes, std::ios::beg);
  in.read(reinterpret_cast<char*>(payload.values.data()), static_cast<std::streamsize>(data_bytes));
  if (!in) {
    throw std::runtime_error("Failed reading legacy restart payload: " + path);
  }
  return payload;
}

LegacyRestartMetadata read_legacy_fortran_state_for_tests(const std::string& path, DnsState& state) {
  const auto payload = read_legacy_payload(path);
  const auto& header = payload.header;
  const auto& file_values = payload.values;

  const int expected_nx = static_cast<int>(header.nx);
  const int expected_ny = static_cast<int>(header.ny);
  const int expected_nz = static_cast<int>(header.nz);

  const auto current_grid = state.grid();
  if (current_grid.components > 0 || current_grid.nx > 0 || current_grid.ny > 0 || current_grid.nz > 0) {
    if (current_grid.nx != expected_nx || current_grid.ny != expected_ny || current_grid.nz != expected_nz) {
      throw std::runtime_error("Legacy fixture grid mismatch with target DnsState: " + path);
    }
    if (current_grid.components != header.components) {
      throw std::runtime_error("Legacy fixture component mismatch with target DnsState: " + path);
    }
  }

  state.resize({expected_nx, expected_ny, expected_nz, header.components});

  const std::size_t native_line_count =
      static_cast<std::size_t>(expected_nx) * static_cast<std::size_t>(expected_nz);
  std::vector<Complex> host_values(state.values_per_component(), Complex(0.0, 0.0));
  for (int component = 0; component < header.components; ++component) {
    std::fill(host_values.begin(), host_values.end(), Complex(0.0, 0.0));
    for (int y = 0; y < expected_ny; ++y) {
      const int legacy_y = y + 1; // File y index 0 is Fortran iy=-1 for one-rank fixtures.
      const auto native_row = static_cast<std::size_t>(y) * native_line_count;
      for (int iz = 0; iz < expected_nz; ++iz) {
        const int signed_iz = iz <= expected_nz / 2 ? iz : iz - expected_nz;
        const int legacy_iz = signed_iz + expected_nz;
        for (int ix = 0; ix < expected_nx; ++ix) {
          const std::size_t legacy_index = legacy_fortran_index(header, legacy_y, legacy_iz, ix, component);
          const std::size_t line = static_cast<std::size_t>(iz) * expected_nx + static_cast<std::size_t>(ix);
          host_values[native_row + line] = file_values[legacy_index];
        }
      }
    }
    state.copy_component_from_host(component, host_values);
  }
  return header;
}

LegacyRestartMetadata read_legacy_fortran_state_with_ghosts_for_tests(const std::string& path, DnsState& state) {
  const auto payload = read_legacy_payload(path);
  const auto& header = payload.header;
  const auto& file_values = payload.values;

  const int expected_nx = static_cast<int>(header.nx);
  const int expected_ny = static_cast<int>(header.ny);
  const int expected_nz = static_cast<int>(header.nz);

  if (expected_nx < 1 || expected_ny < 1 || expected_nz < 1) {
    throw std::runtime_error("Legacy fixture has invalid domain metadata: " + path);
  }

  const auto current_grid = state.grid();
  const int storage_ny = expected_ny + 3;
  if (current_grid.components > 0 || current_grid.nx > 0 || current_grid.ny > 0 || current_grid.nz > 0) {
    if (current_grid.nx != expected_nx || current_grid.ny != storage_ny || current_grid.nz != expected_nz) {
      throw std::runtime_error("Legacy ghost fixture grid mismatch with target DnsState: " + path);
    }
    if (current_grid.components != header.components) {
      throw std::runtime_error("Legacy ghost fixture component mismatch with target DnsState: " + path);
    }
  }

  state.resize({expected_nx,
                storage_ny,
                expected_nz,
                header.components,
                -1,
                1,
                expected_ny - 1});

  const std::size_t native_line_count =
      static_cast<std::size_t>(expected_nx) * static_cast<std::size_t>(expected_nz);
  std::vector<Complex> host_values(state.values_per_component(), Complex(0.0, 0.0));
  for (int component = 0; component < header.components; ++component) {
    std::fill(host_values.begin(), host_values.end(), Complex(0.0, 0.0));
    for (int y_storage = 0; y_storage < storage_ny; ++y_storage) {
      const int legacy_y = y_storage;
      const auto native_row = static_cast<std::size_t>(y_storage) * native_line_count;
      for (int iz = 0; iz < expected_nz; ++iz) {
        const int signed_iz = iz <= expected_nz / 2 ? iz : iz - expected_nz;
        const int legacy_iz = signed_iz + expected_nz;
        for (int ix = 0; ix < expected_nx; ++ix) {
          const std::size_t legacy_index = legacy_fortran_index(header, legacy_y, legacy_iz, ix, component);
          const std::size_t line = static_cast<std::size_t>(iz) * expected_nx + static_cast<std::size_t>(ix);
          host_values[native_row + line] = file_values[legacy_index];
        }
      }
    }
    state.copy_component_from_host(component, host_values);
  }
  return header;
}

LegacyRestartMetadata read_legacy_fortran_full_spectrum_with_ghosts_for_tests(const std::string& path, DnsState& state) {
  const auto payload = read_legacy_payload(path);
  const auto& header = payload.header;
  const auto& file_values = payload.values;

  const int expected_nx = static_cast<int>(header.nx);
  const int expected_ny = static_cast<int>(header.ny);
  const int expected_nz = static_cast<int>(header.nz);

  if (expected_nx < 1 || expected_ny < 1 || expected_nz < 1) {
    throw std::runtime_error("Legacy fixture has invalid domain metadata: " + path);
  }

  const int storage_nx = expected_nx + 1;
  const int storage_ny = expected_ny + 3;
  const int storage_nz = 2 * expected_nz + 1;
  const auto current_grid = state.grid();
  if (current_grid.components > 0 || current_grid.nx > 0 || current_grid.ny > 0 || current_grid.nz > 0) {
    if (current_grid.nx != storage_nx || current_grid.ny != storage_ny || current_grid.nz != storage_nz) {
      throw std::runtime_error("Legacy full-spectrum fixture grid mismatch with target DnsState: " + path);
    }
    if (current_grid.components != header.components) {
      throw std::runtime_error("Legacy full-spectrum fixture component mismatch with target DnsState: " + path);
    }
  }

  state.resize({storage_nx,
                storage_ny,
                storage_nz,
                header.components,
                -1,
                1,
                expected_ny - 1});

  const std::size_t native_line_count =
      static_cast<std::size_t>(storage_nx) * static_cast<std::size_t>(storage_nz);
  std::vector<Complex> host_values(state.values_per_component(), Complex(0.0, 0.0));
  for (int component = 0; component < header.components; ++component) {
    std::fill(host_values.begin(), host_values.end(), Complex(0.0, 0.0));
    for (int y_storage = 0; y_storage < storage_ny; ++y_storage) {
      const int legacy_y = y_storage;
      const auto native_row = static_cast<std::size_t>(y_storage) * native_line_count;
      for (int iz = 0; iz < storage_nz; ++iz) {
        const int signed_iz = canonical_signed_z(iz, expected_nz);
        const int legacy_iz = signed_iz + expected_nz;
        for (int ix = 0; ix < storage_nx; ++ix) {
          const std::size_t legacy_index = legacy_fortran_index(header, legacy_y, legacy_iz, ix, component);
          const std::size_t line = static_cast<std::size_t>(iz) * storage_nx + static_cast<std::size_t>(ix);
          host_values[native_row + line] = file_values[legacy_index];
        }
      }
    }
    state.copy_component_from_host(component, host_values);
  }
  return header;
}

} // namespace channel
