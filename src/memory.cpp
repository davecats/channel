#include "channel/memory.hpp"

#include <cstdlib>
#include <mutex>

#if defined(CHANNEL_HAS_UMPIRE_SPACE)
#include "umpire/ResourceManager.hpp"
#include "umpire/strategy/QuickPool.hpp"
#endif

namespace channel {

namespace {

std::string configured_resource_name() {
  if (const char* resource = std::getenv("CHANNEL_UMPIRE_RESOURCE")) {
    if (*resource != '\0') return resource;
  }
#if defined(CHANNEL_HAS_UMPIRE_SPACE)
  if constexpr (Kokkos::SpaceAccessibility<Kokkos::HostSpace, MemorySpace>::accessible) {
    return "HOST";
  } else {
    return "DEVICE";
  }
#else
  return MemorySpace::name();
#endif
}

} // namespace

void initialize_memory_pool() {
#if defined(CHANNEL_HAS_UMPIRE_SPACE)
  static std::once_flag once;
  std::call_once(once, [] {
    auto& rm = umpire::ResourceManager::getInstance();
    const std::string resource_name = configured_resource_name();
    const std::string pool_name = "channel_default_pool";
    if (!rm.isAllocator(pool_name)) {
      rm.makeAllocator<umpire::strategy::QuickPool>(pool_name, rm.getAllocator(resource_name), 128 * 1024 * 1024);
    }
    DefaultMemorySpace::set_allocator(pool_name);
  });
#endif
}

std::string memory_pool_resource_name() {
  return configured_resource_name();
}

const char* default_memory_space_name() {
  return DefaultMemorySpace::name();
}

} // namespace channel

#if defined(CHANNEL_HAS_UMPIRE_SPACE)
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include <impl/Kokkos_SharedAlloc_timpl.hpp>
#undef KOKKOS_IMPL_PUBLIC_INCLUDE

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(channel::DefaultMemorySpace);
#else
KOKKOS_IMPL_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(channel::DefaultMemorySpace);
#endif
#endif
