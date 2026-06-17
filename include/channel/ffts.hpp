#pragma once

#include "channel/device_vector.hpp"
#include "channel/dnsdata.hpp"
#include "channel/types.hpp"

#include <string>

namespace channel {

enum class FftDirection {
  Forward,
  Inverse,
};

enum class FftNormalization {
  None,
  InverseLength,
  Ortho,
};

struct LocalFftConfig {
  int length = 0;
  int batch_count = 0;
  FftNormalization normalization = FftNormalization::InverseLength;
};

class LocalFftPlan {
public:
  LocalFftPlan() = default;
  explicit LocalFftPlan(LocalFftConfig cfg);

  void configure(LocalFftConfig cfg);
  void execute(DeviceVector<Complex>& values,
               FftDirection direction,
               const std::string& label = "local_fft") const;

  [[nodiscard]] int length() const { return cfg_.length; }
  [[nodiscard]] int batch_count() const { return cfg_.batch_count; }
  [[nodiscard]] FftNormalization normalization() const { return cfg_.normalization; }

private:
  LocalFftConfig cfg_;
};

void local_fft(DeviceVector<Complex>& values,
               const LocalFftConfig& cfg,
               FftDirection direction,
               const std::string& label = "local_fft");

struct DnsComponentFftConfig {
  int component = 0;
  FftNormalization normalization = FftNormalization::InverseLength;
};

class DnsComponentFftPlan {
public:
  void configure(const DnsState& state, DnsComponentFftConfig cfg);
  void execute(DnsState& state,
               FftDirection direction,
               const std::string& label = "dns_component_fft");

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsComponentFftConfig& config() const { return cfg_; }
  [[nodiscard]] DeviceVector<Complex>& scratch() { return scratch_; }
  [[nodiscard]] const DeviceVector<Complex>& scratch() const { return scratch_; }

private:
  DnsGrid grid_;
  DnsComponentFftConfig cfg_;
  LocalFftPlan fft_;
  DeviceVector<Complex> scratch_;
};

} // namespace channel
