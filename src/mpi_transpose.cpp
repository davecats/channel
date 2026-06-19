#include "channel/mpi_transpose.hpp"

#include <limits>
#include <stdexcept>
#include <utility>

namespace channel {

MpiTransposePlan::MpiTransposePlan(MpiTransposeConfig cfg) {
  configure(std::move(cfg));
}

void MpiTransposePlan::configure(MpiTransposeConfig cfg) {
  if (cfg.values_per_peer < 1) {
    throw std::runtime_error("MpiTransposePlan requires values_per_peer >= 1");
  }
  cfg_ = std::move(cfg);
  peer_count_ = mpi_size(cfg_.comm);
  if (peer_count_ < 1) throw std::runtime_error("MpiTransposePlan saw an empty communicator");
  send_.resize(cfg_.label + "_send", buffer_size());
  recv_.resize(cfg_.label + "_recv", buffer_size());
}

void MpiTransposePlan::copy_send_from_host(const std::vector<Complex>& values) {
  if (values.size() != buffer_size()) {
    throw std::runtime_error("MpiTransposePlan::copy_send_from_host size mismatch");
  }
  send_.copy_from_host(values);
}

void MpiTransposePlan::execute() {
  if (cfg_.exchange_mode == ExchangeMode::AllGather) {
    allgather_complex_device(send_.view(), cfg_.values_per_peer, recv_.view(), cfg_.comm, cfg_.label);
    return;
  }
  alltoall_complex_device(send_.view(), recv_.view(), cfg_.values_per_peer, cfg_.comm, cfg_.label);
}

std::size_t MpiTransposePlan::buffer_size() const {
  return static_cast<std::size_t>(cfg_.values_per_peer) * static_cast<std::size_t>(peer_count_);
}

void DnsComponentTransposePlan::configure(const DnsState& state, DnsComponentTransposeConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  [[maybe_unused]] const auto component_origin = state.index(cfg.component, 0, 0);
  if (grid.values_per_component() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsComponentTransposePlan component size exceeds int range");
  }
  MpiTransposeConfig transpose_cfg;
  transpose_cfg.values_per_peer = cfg.values_per_peer;
  transpose_cfg.exchange_mode = cfg.exchange_mode;
  transpose_cfg.comm = cfg.comm;
  transpose_cfg.label = cfg.label;
  transpose_.configure(std::move(transpose_cfg));
  if (transpose_.buffer_size() != grid.values_per_component()) {
    throw std::runtime_error("DnsComponentTransposePlan component size must match MPI transpose buffer size");
  }
  grid_ = grid;
  cfg_ = std::move(cfg);
}

void DnsComponentTransposePlan::execute(DnsState& state) {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error("DnsComponentTransposePlan::execute grid mismatch");
  }
  auto component = state.component_view_3d(cfg_.component);
  auto send = transpose_.send_buffer().view();
  auto recv = transpose_.recv_buffer().view();
  const int line_count = static_cast<int>(grid_.line_count());
  const int nx = grid_.nx;
  Kokkos::parallel_for(
      cfg_.label + "_pack", Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(grid_.values_per_component())),
      KOKKOS_LAMBDA(const int i) {
        const int y = i / line_count;
        const int line = i - y * line_count;
        const int z = line / nx;
        const int x = line - z * nx;
        send(i) = component(y, z, x);
      });
  Kokkos::fence(cfg_.label + "_pack");

  transpose_.execute();

  Kokkos::parallel_for(
      cfg_.label + "_unpack", Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(grid_.values_per_component())),
      KOKKOS_LAMBDA(const int i) {
        const int y = i / line_count;
        const int line = i - y * line_count;
        const int z = line / nx;
        const int x = line - z * nx;
        component(y, z, x) = recv(i);
      });
  Kokkos::fence(cfg_.label + "_unpack");
}

} // namespace channel
