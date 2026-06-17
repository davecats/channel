#include "channel/mpi_transpose.hpp"
#include "channel/memory.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <Kokkos_Core.hpp>

#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

namespace {

channel::Complex payload(int source, int target, int value) {
  return {1000.0 * source + 100.0 * target + value, -10.0 * source - target - 0.25 * value};
}

void require_buffers_use_default_memory(channel::MpiTransposePlan& plan) {
  using SendView = decltype(plan.send_buffer().view());
  using RecvView = decltype(plan.recv_buffer().view());
  static_assert(std::is_same_v<typename SendView::memory_space, channel::DefaultMemorySpace>);
  static_assert(std::is_same_v<typename RecvView::memory_space, channel::DefaultMemorySpace>);

  channel::test::require(plan.send_buffer().data() != nullptr, "transpose send buffer is allocated");
  channel::test::require(plan.recv_buffer().data() != nullptr, "transpose recv buffer is allocated");

  if constexpr (!Kokkos::SpaceAccessibility<Kokkos::HostSpace, channel::MemorySpace>::accessible) {
    channel::test::require(channel::memory_pool_resource_name() == "DEVICE",
                           "device execution defaults MPI transpose buffers to the Umpire DEVICE resource");
  }
}

void exercise_alltoall(const channel::Runtime& runtime, int values_per_peer) {
  channel::MpiTransposePlan plan(
      {values_per_peer, channel::ExchangeMode::AllToAll, channel::world_comm(), "transpose_alltoall"});
  require_buffers_use_default_memory(plan);

  std::vector<channel::Complex> send(plan.buffer_size());
  for (int target = 0; target < runtime.size(); ++target) {
    for (int i = 0; i < values_per_peer; ++i) {
      send[static_cast<std::size_t>(target * values_per_peer + i)] = payload(runtime.rank(), target, i);
    }
  }

  plan.copy_send_from_host(send);
  plan.execute();

  const auto recv = plan.recv_host();
  channel::test::require(recv.size() == send.size(), "AllToAll transpose recv size");
  for (int source = 0; source < runtime.size(); ++source) {
    for (int i = 0; i < values_per_peer; ++i) {
      const auto got = recv[static_cast<std::size_t>(source * values_per_peer + i)];
      const auto want = payload(source, runtime.rank(), i);
      channel::test::require_near(got, want, 0.0, "AllToAll transpose rank-routed payload");
    }
  }
}

void exercise_allgather(const channel::Runtime& runtime, int values_per_peer) {
  channel::MpiTransposePlan plan(
      {values_per_peer, channel::ExchangeMode::AllGather, channel::world_comm(), "transpose_allgather"});
  require_buffers_use_default_memory(plan);

  std::vector<channel::Complex> send(plan.buffer_size(), {-1.0, -1.0});
  for (int i = 0; i < values_per_peer; ++i) {
    send[static_cast<std::size_t>(i)] = payload(runtime.rank(), 0, i);
  }

  plan.copy_send_from_host(send);
  plan.execute();

  const auto recv = plan.recv_host();
  channel::test::require(recv.size() == send.size(), "AllGather transpose recv size");
  for (int source = 0; source < runtime.size(); ++source) {
    for (int i = 0; i < values_per_peer; ++i) {
      const auto got = recv[static_cast<std::size_t>(source * values_per_peer + i)];
      const auto want = payload(source, 0, i);
      channel::test::require_near(got, want, 0.0, "AllGather transpose rank-routed payload");
    }
  }
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(sizeof(channel::Complex) == 2 * sizeof(double),
                         "Complex payload size matches MPI_C_DOUBLE_COMPLEX");
  channel::test::require(!std::string(channel::default_memory_space_name()).empty(),
                         "DefaultMemorySpace reports a backend name");

  exercise_alltoall(runtime, 7);
  exercise_allgather(runtime, 5);

  if (runtime.rank() == 0) std::cout << "MPI transpose test PASSED\n";
  return 0;
}
