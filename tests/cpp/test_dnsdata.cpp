#include "channel/disabled_hooks.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::DnsState state;
  state.resize({4, 5, 3, 3});
  state.fill(channel::Complex(-1.0, 0.0));

  std::vector<channel::Complex> values(state.values_per_component());
  for (std::size_t i = 0; i < values.size(); ++i) {
    values[i] = channel::Complex(static_cast<double>(i), -static_cast<double>(i));
  }

  state.copy_component_from_host(1, values);
  const auto got = state.component_host(1);
  channel::test::require(got.size() == values.size(), "component host size");
  for (std::size_t i = 0; i < got.size(); ++i) {
    channel::test::require_near(got[i], values[i], 0.0, "component value roundtrip");
  }
  for (const int component : {0, 2}) {
    const auto untouched = state.component_host(component);
    for (std::size_t i = 0; i < untouched.size(); ++i) {
      channel::test::require_near(untouched[i], channel::Complex(-1.0, 0.0), 0.0,
                                  "component copy leaves neighbors untouched");
    }
  }

  auto component_two = state.component_view(2);
  Kokkos::parallel_for(
      "dns_state_component_subview_fill", Kokkos::RangePolicy<channel::ExecutionSpace>(0, component_two.extent_int(0)),
      KOKKOS_LAMBDA(const int i) { component_two(i) = channel::Complex(3.0 + 0.1 * i, -0.5); });
  channel::fence("dns_state_component_subview_fill");
  const auto changed = state.component_host(2);
  for (std::size_t i = 0; i < changed.size(); ++i) {
    channel::test::require_near(changed[i], channel::Complex(3.0 + 0.1 * static_cast<double>(i), -0.5), 1.0e-12,
                                "component subview kernel writes selected component");
  }

  const auto idx = state.index(1, 2, 7);
  const auto expected = state.values_per_component() + 2 * state.line_count() + 7;
  channel::test::require(idx == expected, "dns state line-major y index");

  const auto conv = channel::convvelo_hook(state);
  const auto pressure = channel::pressure_hook(state);
  channel::test::require(!conv.executed && !pressure.executed, "disabled hooks stay no-op");
  if (runtime.rank() == 0) std::cout << "DnsState test PASSED\n";
  return 0;
}
