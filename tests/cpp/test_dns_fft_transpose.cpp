#include "channel/ffts.hpp"
#include "channel/mpi_transpose.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace {

int dns_index(int y, int line, int batch) {
  return line + y * batch;
}

channel::Complex fft_value(int y, int line) {
  return channel::Complex(0.35 * y + 0.2 * line, -0.12 * y + 0.5 * line);
}

channel::Complex payload(int source, int target, int value) {
  return {1000.0 * source + 100.0 * target + value, -10.0 * source - target - 0.25 * value};
}

void exercise_dns_fft(const channel::Runtime& runtime) {
  constexpr int nx = 3;
  constexpr int ny = 8;
  constexpr int nz = 2;
  constexpr int batch = nx * nz;
  channel::DnsState state;
  state.resize({nx, ny, nz, 2});
  std::vector<channel::Complex> values(state.values_per_component());
  std::vector<channel::Complex> sentinel(state.values_per_component(), channel::Complex(-7.0, 2.0));

  for (int y = 0; y < ny; ++y) {
    for (int line = 0; line < batch; ++line) {
      values[static_cast<std::size_t>(dns_index(y, line, batch))] = fft_value(y, line);
    }
  }
  state.copy_component_from_host(0, values);
  state.copy_component_from_host(1, sentinel);

  channel::DnsComponentFftPlan plan;
  plan.configure(state, {0, channel::FftNormalization::InverseLength});
  plan.execute(state, channel::FftDirection::Forward, "dns_fft_forward");
  plan.execute(state, channel::FftDirection::Inverse, "dns_fft_inverse");

  const auto roundtrip = state.component_host(0);
  for (std::size_t i = 0; i < values.size(); ++i) {
    channel::test::require_near(roundtrip[i], values[i], 1.0e-11, "DNS FFT component roundtrip");
  }
  const auto untouched = state.component_host(1);
  for (std::size_t i = 0; i < untouched.size(); ++i) {
    channel::test::require_near(untouched[i], sentinel[i], 0.0, "DNS FFT inactive component untouched");
  }

  std::vector<channel::Complex> impulse(state.values_per_component(), channel::Complex(0.0, 0.0));
  for (int y = 0; y < ny; ++y) {
    impulse[static_cast<std::size_t>(dns_index(y, 0, batch))] =
        channel::Complex(1.0 + runtime.rank() + 0.1 * y, 0.0);
  }
  state.copy_component_from_host(0, impulse);
  plan.execute(state, channel::FftDirection::Forward, "dns_fft_impulse");
  const auto spectrum = state.component_host(0);
  for (int y = 0; y < ny; ++y) {
    for (int line = 0; line < batch; ++line) {
      const auto expected = impulse[static_cast<std::size_t>(dns_index(y, 0, batch))];
      channel::test::require_near(spectrum[static_cast<std::size_t>(dns_index(y, line, batch))],
                                  expected, 1.0e-12, "DNS FFT slab impulse spectrum");
    }
  }
}

void exercise_dns_transpose(const channel::Runtime& runtime,
                            channel::ExchangeMode mode,
                            const std::string& label,
                            int values_per_peer) {
  const int total = values_per_peer * runtime.size();
  channel::DnsState state;
  state.resize({total, 1, 1, 2});
  std::vector<channel::Complex> component(static_cast<std::size_t>(total), channel::Complex(-1.0, -1.0));
  std::vector<channel::Complex> sentinel(static_cast<std::size_t>(total), channel::Complex(4.0, -3.0));

  if (mode == channel::ExchangeMode::AllGather) {
    for (int i = 0; i < values_per_peer; ++i) {
      component[static_cast<std::size_t>(i)] = payload(runtime.rank(), 0, i);
    }
  } else {
    for (int target = 0; target < runtime.size(); ++target) {
      for (int i = 0; i < values_per_peer; ++i) {
        component[static_cast<std::size_t>(target * values_per_peer + i)] =
            payload(runtime.rank(), target, i);
      }
    }
  }

  state.copy_component_from_host(0, component);
  state.copy_component_from_host(1, sentinel);

  channel::DnsComponentTransposePlan plan;
  plan.configure(state, {0, values_per_peer, mode, channel::world_comm(), label});
  plan.execute(state);

  const auto got = state.component_host(0);
  for (int source = 0; source < runtime.size(); ++source) {
    for (int i = 0; i < values_per_peer; ++i) {
      const auto want = mode == channel::ExchangeMode::AllGather ? payload(source, 0, i)
                                                                 : payload(source, runtime.rank(), i);
      channel::test::require_near(got[static_cast<std::size_t>(source * values_per_peer + i)], want, 0.0,
                                  label + " DNS component transpose payload");
    }
  }

  const auto untouched = state.component_host(1);
  for (std::size_t i = 0; i < untouched.size(); ++i) {
    channel::test::require_near(untouched[i], sentinel[i], 0.0, label + " inactive component untouched");
  }
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();

  exercise_dns_fft(runtime);
  exercise_dns_transpose(runtime, channel::ExchangeMode::AllToAll, "dns_transpose_alltoall", 5);
  exercise_dns_transpose(runtime, channel::ExchangeMode::AllGather, "dns_transpose_allgather", 4);

  if (runtime.rank() == 0) std::cout << "DNS FFT/transpose test PASSED\n";
  return 0;
}
