#include "channel/dnsdata.hpp"

#include <algorithm>

namespace channel {

void DnsGrid::validate() const {
  if (nx < 1 || ny < 1 || nz < 1 || components < 1) {
    throw std::runtime_error("DnsGrid requires positive nx, ny, nz, and components");
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

void DnsState::resize(const DnsGrid& grid) {
  grid.validate();
  grid_ = grid;
  values_.resize("dns_state", grid_.total_values());
}

void DnsState::fill(const Complex& value) {
  values_.fill(value);
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

Complex* DnsState::component_data(int component) {
  check_component(component);
  return values_.data() + static_cast<std::size_t>(component) * grid_.values_per_component();
}

const Complex* DnsState::component_data(int component) const {
  check_component(component);
  return values_.data() + static_cast<std::size_t>(component) * grid_.values_per_component();
}

void DnsState::copy_component_from_host(int component, std::span<const Complex> values) {
  check_component(component);
  if (values.size() != grid_.values_per_component()) {
    throw std::runtime_error("DnsState::copy_component_from_host size mismatch");
  }
  auto component_values = component_view(component);
  auto mirror = Kokkos::create_mirror_view(component_values);
  std::copy(values.begin(), values.end(), mirror.data());
  Kokkos::deep_copy(component_values, mirror);
}

std::vector<Complex> DnsState::component_host(int component) const {
  check_component(component);
  auto component_values = component_view(component);
  auto mirror = Kokkos::create_mirror_view(component_values);
  Kokkos::deep_copy(mirror, component_values);
  std::vector<Complex> values(grid_.values_per_component());
  std::copy(mirror.data(), mirror.data() + static_cast<std::ptrdiff_t>(values.size()), values.begin());
  return values;
}

void DnsState::check_component(int component) const {
  if (grid_.components < 1) throw std::runtime_error("DnsState is not initialized");
  if (component < 0 || component >= grid_.components) {
    throw std::runtime_error("DnsState component index out of range");
  }
}

} // namespace channel
