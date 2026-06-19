#pragma once

#include "channel/config.hpp"
#include "channel/types.hpp"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iostream>
#include <string>

namespace channel::test {

inline double magnitude(const Complex& value) {
  return std::hypot(static_cast<double>(value.real()), static_cast<double>(value.imag()));
}

inline void require(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAILED: " << message << '\n';
    std::exit(1);
  }
}

inline void require_kokkos_runtime() {
  require(Kokkos::is_initialized(), "Kokkos runtime is initialized");
}

inline void require_throws(const std::function<void()>& fn, const std::string& message) {
  try {
    fn();
  } catch (const std::exception&) {
    return;
  }
  std::cerr << "FAILED: " << message << " did not throw\n";
  std::exit(1);
}

inline void require_near(const Complex& got, const Complex& expected, double tol, const std::string& message) {
  if (magnitude(got - expected) > tol) {
    std::cerr << "FAILED: " << message << " got=(" << got.real() << "," << got.imag()
              << ") expected=(" << expected.real() << "," << expected.imag() << ")\n";
    std::exit(1);
  }
}

} // namespace channel::test
