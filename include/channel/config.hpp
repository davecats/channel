#pragma once

#include <Kokkos_Core.hpp>

#if defined(CHANNEL_HAS_KOKKOSFFT)
#include <KokkosFFT.hpp>
#endif

#if defined(CHANNEL_HAS_UMPIRE_SPACE)
#include <UmpireSpace.hpp>
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include <impl/Kokkos_SharedAlloc_timpl.hpp>
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

namespace channel {

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = ExecutionSpace::memory_space;

#if defined(CHANNEL_HAS_UMPIRE_SPACE)
struct DefaultUmpireMemoryTag {};
using DefaultMemorySpace = ::UmpireSpace<MemorySpace, DefaultUmpireMemoryTag>;
#else
using DefaultMemorySpace = MemorySpace;
#endif

} // namespace channel

#if defined(CHANNEL_HAS_UMPIRE_SPACE)
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_SYCL)
KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_SPECIALIZATION(channel::DefaultMemorySpace);
#else
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(channel::DefaultMemorySpace);
#endif
#endif
