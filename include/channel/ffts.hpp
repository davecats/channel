#pragma once

#include "channel/dnsdata.hpp"
#include "channel/mpi.hpp"
#include "channel/types.hpp"

#include <memory>
#include <string>

namespace channel {

using RealView3D = Kokkos::View<double***, DefaultMemorySpace>;
using ComplexView3D = Kokkos::View<Complex***, DefaultMemorySpace>;

using ComplexFft1DPlan = KokkosFFT::Plan<ExecutionSpace, ComplexView3D, ComplexView3D, 1>;
using ComplexFft2DPlan = KokkosFFT::Plan<ExecutionSpace, ComplexView3D, ComplexView3D, 2>;
using RealToComplexFft1DPlan = KokkosFFT::Plan<ExecutionSpace, RealView3D, ComplexView3D, 1>;
using ComplexToRealFft1DPlan = KokkosFFT::Plan<ExecutionSpace, ComplexView3D, RealView3D, 1>;
using RealToComplexFft2DPlan = KokkosFFT::Plan<ExecutionSpace, RealView3D, ComplexView3D, 2>;
using ComplexToRealFft2DPlan = KokkosFFT::Plan<ExecutionSpace, ComplexView3D, RealView3D, 2>;

enum class FftDirection {
  Forward,
  Inverse,
};

enum class FftNormalization {
  None,
  InverseLength,
  Ortho,
};

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
  [[nodiscard]] ComplexView3D scratch() const { return scratch_; }

private:
  DnsGrid grid_;
  DnsComponentFftConfig cfg_;
  ComplexView3D scratch_;
  std::unique_ptr<ComplexFft2DPlan> forward_plan_;
  std::unique_ptr<ComplexFft2DPlan> inverse_plan_;
};

struct DistributedDealiasedFft2DConfig {
  int spectral_x_total = 0;
  int spectral_x_first = 0;
  int physical_x = 0;
  int padded_z = 0;
  int npxz = 1;
  int ipxz = 0;
  MpiComm comm_x = world_comm();
};

class DistributedDealiasedFft2DPlan {
public:
  void configure(const DnsGrid& compact_grid, DistributedDealiasedFft2DConfig cfg);
  void inverse_component_to_physical(const DnsState& state,
                                     int component,
                                     RealView3D physical,
                                     const std::string& label);
  void forward_physical_to_component(DnsState& state,
                                     RealView3D physical,
                                     int component,
                                     const std::string& label);

  [[nodiscard]] const DnsGrid& compact_grid() const { return grid_; }
  [[nodiscard]] const DistributedDealiasedFft2DConfig& config() const { return cfg_; }
  [[nodiscard]] int z_count() const { return z_count_; }
  [[nodiscard]] std::size_t physical_values() const;

  using ComplexView1D = Kokkos::View<Complex*, DefaultMemorySpace>;

  // Public because CUDA extended lambdas cannot be enclosed by private member functions.
  void pack_component_to_vz(const DnsState& state, int component, const std::string& label);
  void unpack_vz_to_component(DnsState& state, int component, const std::string& label);
  void transpose_z_to_x(const std::string& label);
  void transpose_x_to_z(const std::string& label);

private:
  void check_state(const DnsState& state, const char* caller) const;
  void check_physical(RealView3D physical, const char* caller) const;

  DnsGrid grid_;
  DistributedDealiasedFft2DConfig cfg_;
  int x_half_ = 0;
  int z_count_ = 0;
  int sendcount_ = 0;
  ComplexView3D vz_;
  ComplexView3D vz_scratch_;
  ComplexView3D vx_;
  ComplexView3D vx_scratch_;
  RealView3D physical_plan_view_;
  ComplexView1D send_;
  ComplexView1D recv_;
  std::unique_ptr<ComplexFft1DPlan> z_inverse_plan_;
  std::unique_ptr<RealToComplexFft1DPlan> x_forward_plan_;
  std::unique_ptr<ComplexToRealFft1DPlan> x_inverse_plan_;
  std::unique_ptr<ComplexFft1DPlan> z_forward_plan_;
  bool configured_ = false;
};

} // namespace channel
