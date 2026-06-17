#include "channel/ffts.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  constexpr int length = 8;
  constexpr int batch = 3;
  channel::DeviceVector<channel::Complex> values("fft_values", length * batch);

  std::vector<channel::Complex> host(length * batch);
  for (int line = 0; line < batch; ++line) {
    for (int y = 0; y < length; ++y) {
      host[static_cast<std::size_t>(line + y * batch)] =
          channel::Complex(0.25 * y + line, 0.5 * line - 0.125 * y);
    }
  }

  values.copy_from_host(host);
  channel::LocalFftPlan plan({length, batch, channel::FftNormalization::InverseLength});
  plan.execute(values, channel::FftDirection::Forward);
  plan.execute(values, channel::FftDirection::Inverse);

  const auto roundtrip = values.copy_to_host();
  for (std::size_t i = 0; i < host.size(); ++i) {
    channel::test::require_near(roundtrip[i], host[i], 1.0e-11, "fft forward/inverse roundtrip");
  }

  std::vector<channel::Complex> impulse(length * batch, channel::Complex(0.0, 0.0));
  for (int line = 0; line < batch; ++line) {
    impulse[static_cast<std::size_t>(line)] = channel::Complex(1.0 + line, 0.0);
  }
  values.copy_from_host(impulse);
  plan.execute(values, channel::FftDirection::Forward);
  const auto spectrum = values.copy_to_host();
  for (int line = 0; line < batch; ++line) {
    for (int k = 0; k < length; ++k) {
      channel::test::require_near(
          spectrum[static_cast<std::size_t>(line + k * batch)],
          channel::Complex(1.0 + line, 0.0),
          1.0e-12,
          "fft impulse spectrum");
    }
  }
  if (runtime.rank() == 0) std::cout << "FFT test PASSED\n";
  return 0;
}
