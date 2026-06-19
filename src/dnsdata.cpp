#include "channel/dnsdata.hpp"

#include <string>

namespace channel {

void DnsGrid::validate() const {
  if (nx < 1 || ny < 1 || nz < 1 || components < 1) {
    throw std::runtime_error("DnsGrid requires positive nx, ny, nz, and components");
  }
  if (active_y_count < 0) {
    throw std::runtime_error("DnsGrid active_y_count must be nonnegative");
  }
  const int active_count = active_y_count == 0 ? ny : active_y_count;
  const int active_first = active_y_count == 0 ? y_first : active_y_first;
  const int active_storage_first = active_first - y_first;
  if (active_storage_first < 0 || active_storage_first + active_count > ny) {
    throw std::runtime_error("DnsGrid active y range is outside storage");
  }
}

std::size_t DnsGrid::line_count() const {
  validate();
  return static_cast<std::size_t>(nx) * static_cast<std::size_t>(nz);
}

std::size_t DnsGrid::values_per_component() const {
  return line_count() * static_cast<std::size_t>(ny);
}

std::size_t DnsGrid::total_values() const {
  return values_per_component() * static_cast<std::size_t>(components);
}

int DnsGrid::active_y_storage_first() const {
  validate();
  return active_y_count == 0 ? 0 : active_y_first - y_first;
}

int DnsGrid::active_y_storage_last_exclusive() const {
  validate();
  return active_y_storage_first() + (active_y_count == 0 ? ny : active_y_count);
}

std::size_t DnsGrid::active_values_per_component() const {
  validate();
  const int count = active_y_count == 0 ? ny : active_y_count;
  return line_count() * static_cast<std::size_t>(count);
}

void DnsState::resize(const DnsGrid& grid) {
  grid.validate();
  grid_ = grid;
  components_.clear();
  components_.reserve(static_cast<std::size_t>(grid_.components));
  for (int component = 0; component < grid_.components; ++component) {
    components_.emplace_back("dns_state_component_" + std::to_string(component), grid_.ny, grid_.nz, grid_.nx);
  }
}

void DnsState::fill(const Complex& value) {
  for (auto& component : components_) {
    Kokkos::deep_copy(component, value);
  }
}

std::size_t DnsState::index(int component, int y, int line) const {
  check_component(component);
  if (y < 0 || y >= grid_.ny) throw std::runtime_error("DnsState y index out of range");
  if (line < 0 || static_cast<std::size_t>(line) >= grid_.line_count()) {
    throw std::runtime_error("DnsState line index out of range");
  }
  const auto component_offset = static_cast<std::size_t>(component) * grid_.values_per_component();
  return component_offset + static_cast<std::size_t>(y) * grid_.line_count() +
         static_cast<std::size_t>(line);
}

void DnsState::copy_component_from_host(int component, std::span<const Complex> values) {
  check_component(component);
  if (values.size() != grid_.values_per_component()) {
    throw std::runtime_error("DnsState::copy_component_from_host size mismatch");
  }
  auto component_values = component_view_3d(component);
  auto mirror = Kokkos::create_mirror_view(component_values);
  for (int y = 0; y < grid_.ny; ++y) {
    for (int z = 0; z < grid_.nz; ++z) {
      for (int x = 0; x < grid_.nx; ++x) {
        const auto line = static_cast<std::size_t>(z) * static_cast<std::size_t>(grid_.nx) +
                          static_cast<std::size_t>(x);
        mirror(y, z, x) = values[static_cast<std::size_t>(y) * grid_.line_count() + line];
      }
    }
  }
  Kokkos::deep_copy(component_values, mirror);
}

std::vector<Complex> DnsState::component_host(int component) const {
  check_component(component);
  auto component_values = component_view_3d(component);
  auto mirror = Kokkos::create_mirror_view(component_values);
  Kokkos::deep_copy(mirror, component_values);
  std::vector<Complex> values(grid_.values_per_component());
  for (int y = 0; y < grid_.ny; ++y) {
    for (int z = 0; z < grid_.nz; ++z) {
      for (int x = 0; x < grid_.nx; ++x) {
        const auto line = static_cast<std::size_t>(z) * static_cast<std::size_t>(grid_.nx) +
                          static_cast<std::size_t>(x);
        values[static_cast<std::size_t>(y) * grid_.line_count() + line] = mirror(y, z, x);
      }
    }
  }
  return values;
}

void DnsState::check_component(int component) const {
  if (grid_.components < 1 || components_.empty()) throw std::runtime_error("DnsState is not initialized");
  if (component < 0 || component >= grid_.components) {
    throw std::runtime_error("DnsState component index out of range");
  }
}

} // namespace channel
