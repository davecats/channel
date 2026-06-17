#include "channel/ffts.hpp"

#include "channel/memory.hpp"

#include <limits>
#include <stdexcept>

namespace channel {

namespace {

KokkosFFT::Normalization kokkosfft_normalization(FftNormalization normalization) {
  switch (normalization) {
    case FftNormalization::None:
      return KokkosFFT::Normalization::none;
    case FftNormalization::InverseLength:
      return KokkosFFT::Normalization::backward;
    case FftNormalization::Ortho:
      return KokkosFFT::Normalization::ortho;
  }
  return KokkosFFT::Normalization::backward;
}

} // namespace

LocalFftPlan::LocalFftPlan(LocalFftConfig cfg) {
  configure(cfg);
}

void LocalFftPlan::configure(LocalFftConfig cfg) {
  if (cfg.length < 1 || cfg.batch_count < 1) {
    throw std::runtime_error("LocalFftPlan requires positive length and batch_count");
  }
  cfg_ = cfg;
}

void LocalFftPlan::execute(DeviceVector<Complex>& values,
                           FftDirection direction,
                           const std::string& label) const {
  local_fft(values, cfg_, direction, label);
}

void local_fft(DeviceVector<Complex>& values,
               const LocalFftConfig& cfg,
               FftDirection direction,
               const std::string& label) {
  (void)label;
  if (cfg.length < 1 || cfg.batch_count < 1) {
    throw std::runtime_error("local_fft requires positive length and batch_count");
  }
  const auto expected = static_cast<std::size_t>(cfg.length) * static_cast<std::size_t>(cfg.batch_count);
  if (values.size() != expected) {
    throw std::runtime_error("local_fft input size mismatch");
  }

  using UnmanagedMatrix =
      Kokkos::View<Complex**, Kokkos::LayoutRight, DefaultMemorySpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  UnmanagedMatrix input(values.data(), cfg.length, cfg.batch_count);
  initialize_memory_pool();
  Kokkos::View<Complex**, Kokkos::LayoutRight, DefaultMemorySpace> output(label + "_tmp", cfg.length, cfg.batch_count);
  const auto norm = kokkosfft_normalization(cfg.normalization);

  ExecutionSpace exec;
  if (direction == FftDirection::Forward) {
    KokkosFFT::fft(exec, input, output, norm, 0);
  } else {
    KokkosFFT::ifft(exec, input, output, norm, 0);
  }
  Kokkos::deep_copy(exec, input, output);
  exec.fence(label);
}

void DnsComponentFftPlan::configure(const DnsState& state, DnsComponentFftConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  [[maybe_unused]] const auto component_origin = state.index(cfg.component, 0, 0);
  if (grid.line_count() > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      grid.values_per_component() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsComponentFftPlan grid exceeds int range");
  }
  grid_ = grid;
  cfg_ = cfg;
  scratch_.resize("dns_component_fft_scratch", grid_.values_per_component());
  fft_.configure({grid_.ny, static_cast<int>(grid_.line_count()), cfg_.normalization});
}

void DnsComponentFftPlan::execute(DnsState& state, FftDirection direction, const std::string& label) {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error("DnsComponentFftPlan::execute grid mismatch");
  }
  auto component = state.component_view(cfg_.component);
  auto scratch = scratch_.view();
  Kokkos::parallel_for(
      label + "_pack", Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(grid_.values_per_component())),
      KOKKOS_LAMBDA(const int i) { scratch(i) = component(i); });
  Kokkos::fence(label + "_pack");

  fft_.execute(scratch_, direction, label);

  Kokkos::parallel_for(
      label + "_unpack", Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(grid_.values_per_component())),
      KOKKOS_LAMBDA(const int i) { component(i) = scratch(i); });
  Kokkos::fence(label + "_unpack");
}

} // namespace channel
