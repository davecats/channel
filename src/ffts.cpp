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
  initialize_memory_pool();
  scratch_ = ComplexView3D("dns_component_fft_scratch", grid_.ny, grid_.nz, grid_.nx);
  auto component = const_cast<DnsState&>(state).component_view_3d(cfg_.component);
  ExecutionSpace exec;
  forward_plan_ = std::make_unique<ComplexFft2DPlan>(
      exec, scratch_, component, KokkosFFT::Direction::forward, KokkosFFT::axis_type<2>({1, 2}));
  inverse_plan_ = std::make_unique<ComplexFft2DPlan>(
      exec, scratch_, component, KokkosFFT::Direction::backward, KokkosFFT::axis_type<2>({1, 2}));
}

void DnsComponentFftPlan::execute(DnsState& state, FftDirection direction, const std::string& label) {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error("DnsComponentFftPlan::execute grid mismatch");
  }
  auto component = state.component_view_3d(cfg_.component);
  auto scratch = scratch_;
  Kokkos::parallel_for(
      label + "_pack",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {grid_.ny, grid_.nz, grid_.nx}),
      KOKKOS_LAMBDA(const int y, const int z, const int x) { scratch(y, z, x) = component(y, z, x); });
  Kokkos::fence(label + "_pack");

  ExecutionSpace exec;
  const auto norm = kokkosfft_normalization(cfg_.normalization);
  if (direction == FftDirection::Forward) {
    KokkosFFT::execute(*forward_plan_, scratch, component, norm);
  } else {
    KokkosFFT::execute(*inverse_plan_, scratch, component, norm);
  }
  exec.fence(label);

}

void DistributedDealiasedFft2DPlan::configure(const DnsGrid& compact_grid, DistributedDealiasedFft2DConfig cfg) {
  compact_grid.validate();
  if (cfg.npxz < 1 || cfg.ipxz < 0 || cfg.ipxz >= cfg.npxz) {
    throw std::runtime_error("DistributedDealiasedFft2DPlan received invalid rank metadata");
  }
  const int comm_size = mpi_size(cfg.comm_x);
  const int comm_rank = mpi_rank(cfg.comm_x);
  if (comm_size != cfg.npxz || comm_rank != cfg.ipxz) {
    throw std::runtime_error("DistributedDealiasedFft2DPlan communicator metadata mismatch");
  }
  if (cfg.physical_x < 1 || cfg.padded_z < 1 || (cfg.physical_x % 2) != 0) {
    throw std::runtime_error("DistributedDealiasedFft2DPlan requires positive padded FFT dimensions");
  }
  x_half_ = cfg.physical_x / 2 + 1;
  if (cfg.spectral_x_total < compact_grid.nx ||
      cfg.spectral_x_first < 0 ||
      cfg.spectral_x_first + compact_grid.nx > cfg.spectral_x_total) {
    throw std::runtime_error("DistributedDealiasedFft2DPlan received invalid compact x ownership");
  }
  if (cfg.spectral_x_total % cfg.npxz != 0 || cfg.padded_z % cfg.npxz != 0) {
    throw std::runtime_error("DistributedDealiasedFft2DPlan requires npxz to divide nx+1 and nzd");
  }
  if (compact_grid.nx != cfg.spectral_x_total / cfg.npxz ||
      cfg.spectral_x_first != cfg.ipxz * compact_grid.nx) {
    throw std::runtime_error("DistributedDealiasedFft2DPlan currently requires equal Fortran-style x blocks");
  }

  grid_ = compact_grid;
  cfg_ = cfg;
  z_count_ = cfg_.padded_z / cfg_.npxz;
  sendcount_ = grid_.ny * grid_.nx * z_count_;
  const auto buffer_size = static_cast<std::size_t>(sendcount_) * static_cast<std::size_t>(cfg_.npxz);
  initialize_memory_pool();
  vz_ = ComplexView3D("distributed_fft2d_vz", grid_.ny, cfg_.padded_z, grid_.nx);
  vz_scratch_ = ComplexView3D("distributed_fft2d_vz_scratch", grid_.ny, cfg_.padded_z, grid_.nx);
  vx_ = ComplexView3D("distributed_fft2d_vx", grid_.ny, z_count_, x_half_);
  vx_scratch_ = ComplexView3D("distributed_fft2d_vx_scratch", grid_.ny, z_count_, x_half_);
  physical_plan_view_ = RealView3D("distributed_fft2d_physical_plan", grid_.ny, z_count_, cfg_.physical_x);
  send_ = ComplexView1D("distributed_fft2d_send", buffer_size);
  recv_ = ComplexView1D("distributed_fft2d_recv", buffer_size);
  ExecutionSpace exec;
  z_inverse_plan_ =
      std::make_unique<ComplexFft1DPlan>(exec, vz_, vz_scratch_, KokkosFFT::Direction::backward, 1);
  x_forward_plan_ =
      std::make_unique<RealToComplexFft1DPlan>(exec, physical_plan_view_, vx_scratch_, KokkosFFT::Direction::forward, 2);
  x_inverse_plan_ =
      std::make_unique<ComplexToRealFft1DPlan>(exec, vx_, physical_plan_view_, KokkosFFT::Direction::backward, 2);
  z_forward_plan_ =
      std::make_unique<ComplexFft1DPlan>(exec, vz_, vz_scratch_, KokkosFFT::Direction::forward, 1);
  configured_ = true;
}

std::size_t DistributedDealiasedFft2DPlan::physical_values() const {
  if (!configured_) return 0;
  return static_cast<std::size_t>(grid_.ny) *
         static_cast<std::size_t>(z_count_) *
         static_cast<std::size_t>(cfg_.physical_x);
}

void DistributedDealiasedFft2DPlan::check_state(const DnsState& state, const char* caller) const {
  if (!configured_) {
    throw std::runtime_error(std::string(caller) + " called before configure");
  }
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DistributedDealiasedFft2DPlan::check_physical(RealView3D physical, const char* caller) const {
  const auto values = static_cast<std::size_t>(physical.extent(0)) *
                      static_cast<std::size_t>(physical.extent(1)) *
                      static_cast<std::size_t>(physical.extent(2));
  if (values != physical_values()) {
    throw std::runtime_error(std::string(caller) + " physical buffer size mismatch");
  }
}

void DistributedDealiasedFft2DPlan::pack_component_to_vz(const DnsState& state,
                                                         int component_id,
                                                         const std::string& label) {
  const auto component = state.component_view_3d(component_id);
  auto vz = vz_;
  Kokkos::deep_copy(vz_, Complex(0.0, 0.0));
  const int retained_z = (grid_.nz - 1) / 2;
  const int storage_z = grid_.nz;
  const int padded_z = cfg_.padded_z;
  Kokkos::parallel_for(
      label + "_pack_component_to_vz",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {grid_.ny, grid_.nz, grid_.nx}),
      KOKKOS_LAMBDA(const int y, const int z, const int x) {
        const int signed_z = z <= retained_z ? z : z - storage_z;
        const int padded_index_z = signed_z >= 0 ? signed_z : signed_z + padded_z;
        vz(y, padded_index_z, x) = component(y, z, x);
      });
  Kokkos::fence((label + "_pack_component_to_vz").c_str());
}

void DistributedDealiasedFft2DPlan::unpack_vz_to_component(DnsState& state, int component_id, const std::string& label) {
  auto component = state.component_view_3d(component_id);
  auto vz_scratch = vz_scratch_;
  const int retained_z = (grid_.nz - 1) / 2;
  const int storage_z = grid_.nz;
  const int padded_z = cfg_.padded_z;
  Kokkos::parallel_for(
      label + "_unpack_vz_to_component",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {grid_.ny, grid_.nz, grid_.nx}),
      KOKKOS_LAMBDA(const int y, const int z, const int x) {
        const int signed_z = z <= retained_z ? z : z - storage_z;
        const int padded_index_z = signed_z >= 0 ? signed_z : signed_z + padded_z;
        component(y, z, x) = vz_scratch(y, padded_index_z, x);
      });
  Kokkos::fence((label + "_unpack_vz_to_component").c_str());
}

void DistributedDealiasedFft2DPlan::transpose_z_to_x(const std::string& label) {
  auto send = send_;
  auto recv = recv_;
  auto vz_scratch = vz_scratch_;
  auto vx = vx_;
  const int nx_local = grid_.nx;
  const int z_local = z_count_;
  const int sendcount = sendcount_;
  Kokkos::parallel_for(
      label + "_z_to_x_pack",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>, ExecutionSpace>({0, 0, 0, 0}, {cfg_.npxz, grid_.ny, nx_local, z_local}),
      KOKKOS_LAMBDA(const int dest, const int y, const int x, const int z) {
        const int p = dest * sendcount + z + z_local * x + z_local * nx_local * y;
        send(p) = vz_scratch(y, dest * z_local + z, x);
      });
  Kokkos::fence((label + "_z_to_x_pack").c_str());
  alltoall_complex_device(send_, recv_, sendcount_, cfg_.comm_x, label + "_z_to_x");
  Kokkos::deep_copy(vx_, Complex(0.0, 0.0));
  Kokkos::parallel_for(
      label + "_z_to_x_unpack",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>, ExecutionSpace>({0, 0, 0, 0}, {cfg_.npxz, grid_.ny, nx_local, z_local}),
      KOKKOS_LAMBDA(const int src, const int y, const int x, const int z) {
        const int p = src * sendcount + z + z_local * x + z_local * nx_local * y;
        vx(y, z, src * nx_local + x) = recv(p);
      });
  Kokkos::fence((label + "_z_to_x_unpack").c_str());
}

void DistributedDealiasedFft2DPlan::transpose_x_to_z(const std::string& label) {
  auto send = send_;
  auto recv = recv_;
  auto vx_scratch = vx_scratch_;
  auto vz = vz_;
  const int nx_local = grid_.nx;
  const int z_local = z_count_;
  const int sendcount = sendcount_;
  Kokkos::parallel_for(
      label + "_x_to_z_pack",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>, ExecutionSpace>({0, 0, 0, 0}, {cfg_.npxz, grid_.ny, z_local, nx_local}),
      KOKKOS_LAMBDA(const int dest, const int y, const int z, const int x) {
        const int p = dest * sendcount + x + nx_local * z + nx_local * z_local * y;
        send(p) = vx_scratch(y, z, dest * nx_local + x);
      });
  Kokkos::fence((label + "_x_to_z_pack").c_str());
  alltoall_complex_device(send_, recv_, sendcount_, cfg_.comm_x, label + "_x_to_z");
  Kokkos::deep_copy(vz_, Complex(0.0, 0.0));
  Kokkos::parallel_for(
      label + "_x_to_z_unpack",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>, ExecutionSpace>({0, 0, 0, 0}, {cfg_.npxz, grid_.ny, z_local, nx_local}),
      KOKKOS_LAMBDA(const int src, const int y, const int z, const int x) {
        const int p = src * sendcount + x + nx_local * z + nx_local * z_local * y;
        vz(y, src * z_local + z, x) = recv(p);
      });
  Kokkos::fence((label + "_x_to_z_unpack").c_str());
}

void DistributedDealiasedFft2DPlan::inverse_component_to_physical(const DnsState& state,
                                                                  int component,
                                                                  RealView3D physical,
                                                                  const std::string& label) {
  check_state(state, "DistributedDealiasedFft2DPlan::inverse_component_to_physical");
  check_physical(physical, "DistributedDealiasedFft2DPlan::inverse_component_to_physical");

  pack_component_to_vz(state, component, label);
  ExecutionSpace exec;
  KokkosFFT::execute(*z_inverse_plan_, vz_, vz_scratch_, KokkosFFT::Normalization::none);
  exec.fence(label + "_z_ifft");
  transpose_z_to_x(label);
  KokkosFFT::execute(*x_inverse_plan_, vx_, physical, KokkosFFT::Normalization::none);
  exec.fence(label + "_x_irfft");
}

void DistributedDealiasedFft2DPlan::forward_physical_to_component(DnsState& state,
                                                                  RealView3D physical,
                                                                  int component,
                                                                  const std::string& label) {
  check_state(state, "DistributedDealiasedFft2DPlan::forward_physical_to_component");
  check_physical(physical, "DistributedDealiasedFft2DPlan::forward_physical_to_component");

  ExecutionSpace exec;
  KokkosFFT::execute(*x_forward_plan_, physical, vx_scratch_, KokkosFFT::Normalization::none);
  exec.fence(label + "_x_rfft");
  transpose_x_to_z(label);
  KokkosFFT::execute(*z_forward_plan_, vz_, vz_scratch_, KokkosFFT::Normalization::none);
  exec.fence(label + "_z_fft");
  unpack_vz_to_component(state, component, label);
}

} // namespace channel
