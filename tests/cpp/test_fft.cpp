#include "channel/ffts.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <KokkosFFT.hpp>

#include <iostream>
#include <vector>

namespace {

void test_fft_axis0(channel::ComplexView3D values,
                    channel::ComplexView3D scratch,
                    channel::FftDirection direction,
                    const char* label) {
  channel::ExecutionSpace exec;
  if (direction == channel::FftDirection::Forward) {
    KokkosFFT::fft(exec, values, scratch, KokkosFFT::Normalization::backward, 0);
  } else {
    KokkosFFT::ifft(exec, values, scratch, KokkosFFT::Normalization::backward, 0);
  }
  Kokkos::deep_copy(exec, values, scratch);
  exec.fence(label);
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  constexpr int length = 8;
  constexpr int nz = 3;
  constexpr int nx = 2;
  channel::ComplexView3D values("fft_values", length, nz, nx);
  channel::ComplexView3D scratch("fft_scratch", length, nz, nx);

  std::vector<channel::Complex> host(length * nz * nx);
  for (int z = 0; z < nz; ++z) {
    for (int x = 0; x < nx; ++x) {
      const int line = z * nx + x;
      for (int y = 0; y < length; ++y) {
        host[static_cast<std::size_t>(line + y * nz * nx)] =
            channel::Complex(0.25 * y + line, 0.5 * line - 0.125 * y);
      }
    }
  }

  auto values_host = Kokkos::create_mirror_view(values);
  for (int y = 0; y < length; ++y) {
    for (int z = 0; z < nz; ++z) {
      for (int x = 0; x < nx; ++x) {
        const int line = z * nx + x;
        values_host(y, z, x) = host[static_cast<std::size_t>(line + y * nz * nx)];
      }
    }
  }
  Kokkos::deep_copy(values, values_host);
  test_fft_axis0(values, scratch, channel::FftDirection::Forward, "test_fft_forward");
  test_fft_axis0(values, scratch, channel::FftDirection::Inverse, "test_fft_inverse");

  Kokkos::deep_copy(values_host, values);
  for (int y = 0; y < length; ++y) {
    for (int z = 0; z < nz; ++z) {
      for (int x = 0; x < nx; ++x) {
        const int line = z * nx + x;
        const auto i = static_cast<std::size_t>(line + y * nz * nx);
        channel::test::require_near(values_host(y, z, x), host[i], 1.0e-11, "fft forward/inverse roundtrip");
      }
    }
  }

  for (int y = 0; y < length; ++y) {
    for (int z = 0; z < nz; ++z) {
      for (int x = 0; x < nx; ++x) {
        const int line = z * nx + x;
        values_host(y, z, x) = y == 0 ? channel::Complex(1.0 + line, 0.0) : channel::Complex(0.0, 0.0);
      }
    }
  }
  Kokkos::deep_copy(values, values_host);
  test_fft_axis0(values, scratch, channel::FftDirection::Forward, "test_fft_impulse_forward");
  Kokkos::deep_copy(values_host, values);
  for (int z = 0; z < nz; ++z) {
    for (int x = 0; x < nx; ++x) {
      const int line = z * nx + x;
      for (int k = 0; k < length; ++k) {
        channel::test::require_near(values_host(k, z, x),
                                    channel::Complex(1.0 + line, 0.0),
                                    1.0e-12,
                                    "fft impulse spectrum");
      }
    }
  }
  if (runtime.rank() == 0) std::cout << "FFT test PASSED\n";
  return 0;
}
