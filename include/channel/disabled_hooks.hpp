#pragma once

#include "channel/dnsdata.hpp"

#include <string>

namespace channel {

struct DisabledHookResult {
  bool executed = false;
  std::string reason;
};

DisabledHookResult convvelo_hook(DnsState& state);
DisabledHookResult pressure_hook(DnsState& state);

} // namespace channel
