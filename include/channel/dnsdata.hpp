#pragma once

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
  int y_first = 0;
  int active_y_first = 0;
  int active_y_count = 0;

  void validate() const;
  [[nodiscard]] std::size_t line_count() const;
  [[nodiscard]] std::size_t values_per_component() const;
  [[nodiscard]] std::size_t total_values() const;
  [[nodiscard]] int active_y_storage_first() const;
  [[nodiscard]] int active_y_storage_last_exclusive() const;
  [[nodiscard]] std::size_t active_values_per_component() const;
};

class DnsState {
public:
  using ComponentView3D = Kokkos::View<Complex***, DefaultMemorySpace>;
  using ConstComponentView3D = Kokkos::View<const Complex***, DefaultMemorySpace>;

  void resize(const DnsGrid& grid);
  void fill(const Complex& value);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] std::size_t line_count() const { return grid_.line_count(); }
  [[nodiscard]] std::size_t values_per_component() const { return grid_.values_per_component(); }
  [[nodiscard]] std::size_t total_values() const { return grid_.total_values(); }

  [[nodiscard]] std::size_t index(int component, int y, int line) const;
  [[nodiscard]] ComponentView3D component_view_3d(int component);
  [[nodiscard]] ConstComponentView3D component_view_3d(int component) const;

  void copy_component_from_host(int component, std::span<const Complex> values);
  [[nodiscard]] std::vector<Complex> component_host(int component) const;

private:
  void check_component(int component) const;

  DnsGrid grid_;
  std::vector<ComponentView3D> components_;
};

inline DnsState::ComponentView3D DnsState::component_view_3d(int component) {
  check_component(component);
  return components_[static_cast<std::size_t>(component)];
}

inline DnsState::ConstComponentView3D DnsState::component_view_3d(int component) const {
  check_component(component);
  return components_[static_cast<std::size_t>(component)];
}

} // namespace channel
