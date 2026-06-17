#pragma once

#include "channel/device_vector.hpp"
#include "channel/types.hpp"

#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace channel {

struct DnsGrid {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  int components = 0;

  void validate() const;
  [[nodiscard]] std::size_t line_count() const;
  [[nodiscard]] std::size_t values_per_component() const;
  [[nodiscard]] std::size_t total_values() const;
};

class DnsState {
public:
  void resize(const DnsGrid& grid);
  void fill(const Complex& value);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] std::size_t line_count() const { return grid_.line_count(); }
  [[nodiscard]] std::size_t values_per_component() const { return grid_.values_per_component(); }
  [[nodiscard]] std::size_t total_values() const { return grid_.total_values(); }

  [[nodiscard]] std::size_t index(int component, int y, int line) const;
  [[nodiscard]] Complex* component_data(int component);
  [[nodiscard]] const Complex* component_data(int component) const;
  [[nodiscard]] auto component_view(int component);
  [[nodiscard]] auto component_view(int component) const;
  [[nodiscard]] DeviceVector<Complex>& values() { return values_; }
  [[nodiscard]] const DeviceVector<Complex>& values() const { return values_; }

  void copy_component_from_host(int component, std::span<const Complex> values);
  [[nodiscard]] std::vector<Complex> component_host(int component) const;

private:
  void check_component(int component) const;

  DnsGrid grid_;
  DeviceVector<Complex> values_;
};

inline auto DnsState::component_view(int component) {
  check_component(component);
  const auto offset = static_cast<std::size_t>(component) * grid_.values_per_component();
  return Kokkos::subview(values_.view(), std::make_pair(offset, offset + grid_.values_per_component()));
}

inline auto DnsState::component_view(int component) const {
  check_component(component);
  const auto offset = static_cast<std::size_t>(component) * grid_.values_per_component();
  return Kokkos::subview(values_.view(), std::make_pair(offset, offset + grid_.values_per_component()));
}

} // namespace channel
