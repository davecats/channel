#include "channel/runtime.hpp"
#include "channel/workspace.hpp"

#include "test_common.hpp"

#include <iostream>

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::WorkspaceArena arena(channel::WorkspaceBackend::UmpireSpace);
  channel::test::require(arena.backend() == channel::WorkspaceBackend::UmpireSpace,
                         "UmpireSpace workspace backend is active");

  {
    auto lease = arena.lease(4096, "test");
    auto* doubles = lease.slice<double>(8);
    auto* complex_values = lease.slice<channel::Complex>(4);
    channel::test::require(doubles != nullptr, "double slice is non-null");
    channel::test::require(complex_values != nullptr, "complex slice is non-null");
    channel::test::require(lease.used() <= lease.bytes(), "lease usage stays inside capacity");
    channel::test::require(arena.owned(), "arena is owned while lease lives");
  }

  channel::test::require(!arena.owned(), "arena releases after lease destruction");
  channel::test::require(arena.high_water_bytes() >= 4096, "arena high-water records lease size");
  if (runtime.rank() == 0) std::cout << "WorkspaceArena test PASSED\n";
  return 0;
}
