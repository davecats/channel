#include "channel/disabled_hooks.hpp"

namespace channel {

DisabledHookResult convvelo_hook(DnsState& state) {
  (void)state;
  return {false, "convvelo is intentionally disabled for the first C++/Kokkos milestone"};
}

DisabledHookResult pressure_hook(DnsState& state) {
  (void)state;
  return {false, "pressure is intentionally disabled for the first C++/Kokkos milestone"};
}

} // namespace channel
