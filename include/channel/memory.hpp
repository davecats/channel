#pragma once

#include "channel/config.hpp"

#include <string>

namespace channel {

void initialize_memory_pool();
[[nodiscard]] std::string memory_pool_resource_name();
[[nodiscard]] const char* default_memory_space_name();

} // namespace channel
