#pragma once

#include "channel/dnsdata.hpp"

#include <cstdint>
#include <string>
#include <string_view>

namespace channel {

struct LegacyRestartMetadata {
  std::int32_t nx = 0;
  std::int32_t ny = 0;
  std::int32_t nz = 0;
  double alfa0 = 0.0;
  double beta0 = 0.0;
  double ni = 0.0;
  double stretching = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  double time = 0.0;
  int components = 0;
};

void write_state_to_file(const DnsState& state, std::string_view path);

LegacyRestartMetadata read_legacy_fortran_state_for_tests(const std::string& path, DnsState& state);
LegacyRestartMetadata read_legacy_fortran_state_with_ghosts_for_tests(const std::string& path, DnsState& state);
LegacyRestartMetadata read_legacy_fortran_full_spectrum_with_ghosts_for_tests(const std::string& path, DnsState& state);

} // namespace channel
