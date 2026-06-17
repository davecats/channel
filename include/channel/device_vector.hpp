#pragma once

#include "channel/config.hpp"
#include "channel/memory.hpp"

#include <algorithm>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace channel {

template <class T>
class DeviceVector {
public:
  DeviceVector() = default;

  DeviceVector(std::string label, std::size_t size) {
    resize(std::move(label), size);
  }

  void resize(std::string label, std::size_t size) {
    label_ = std::move(label);
    initialize_memory_pool();
    view_ = Kokkos::View<T*, DefaultMemorySpace>(label_, size);
  }

  [[nodiscard]] std::size_t size() const {
    return view_.extent(0);
  }

  [[nodiscard]] bool empty() const { return size() == 0; }

  T* data() {
    return view_.data();
  }

  const T* data() const {
    return view_.data();
  }

  auto view() { return view_; }
  auto view() const { return view_; }

  void copy_from_host(std::span<const T> values) {
    if (values.size() != size()) {
      throw std::runtime_error("DeviceVector::copy_from_host size mismatch");
    }
    auto mirror = Kokkos::create_mirror_view(view_);
    for (std::size_t i = 0; i < values.size(); ++i) mirror(i) = values[i];
    Kokkos::deep_copy(view_, mirror);
  }

  [[nodiscard]] std::vector<T> copy_to_host() const {
    std::vector<T> values(size());
    auto mirror = Kokkos::create_mirror_view(view_);
    Kokkos::deep_copy(mirror, view_);
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = mirror(i);
    return values;
  }

  void copy_from(const DeviceVector<T>& other) {
    if (other.size() != size()) {
      throw std::runtime_error("DeviceVector::copy_from size mismatch");
    }
    Kokkos::deep_copy(view_, other.view_);
  }

  void fill(const T& value) {
    Kokkos::deep_copy(view_, value);
  }

private:
  std::string label_;
  Kokkos::View<T*, DefaultMemorySpace> view_;
};

} // namespace channel
