#include "channel/runtime.hpp"
#include "channel/yline.hpp"

#include "test_common.hpp"

#include <iostream>
#include <vector>

namespace {

std::vector<channel::Complex> apply_penta(const std::vector<channel::Complex>& ds,
                                          const std::vector<channel::Complex>& dl,
                                          const std::vector<channel::Complex>& d,
                                          const std::vector<channel::Complex>& du,
                                          const std::vector<channel::Complex>& dw,
                                          const std::vector<channel::Complex>& exact,
                                          int n,
                                          int batch) {
  std::vector<channel::Complex> rhs(exact.size(), channel::Complex(0.0, 0.0));
  for (int line = 0; line < batch; ++line) {
    for (int row = 0; row < n; ++row) {
      const int p = line + row * batch;
      rhs[p] += d[p] * exact[p];
      if (row >= 1) rhs[p] += dl[p] * exact[line + (row - 1) * batch];
      if (row >= 2) rhs[p] += ds[p] * exact[line + (row - 2) * batch];
      if (row + 1 < n) rhs[p] += du[p] * exact[line + (row + 1) * batch];
      if (row + 2 < n) rhs[p] += dw[p] * exact[line + (row + 2) * batch];
    }
  }
  return rhs;
}

void run_case(int n, int batch) {
  const auto count = static_cast<std::size_t>(n) * static_cast<std::size_t>(batch);
  std::vector<channel::Complex> ds(count), dl(count), d(count), du(count), dw(count), exact(count);

  for (int line = 0; line < batch; ++line) {
    for (int row = 0; row < n; ++row) {
      const int p = line + row * batch;
      exact[p] = channel::Complex(1.0 + 0.25 * row + 0.1 * line, -0.2 + 0.05 * row);
      ds[p] = (row >= 2) ? channel::Complex(-0.03, 0.01) : channel::Complex(0.0, 0.0);
      dl[p] = (row >= 1) ? channel::Complex(-0.15, -0.02) : channel::Complex(0.0, 0.0);
      d[p] = channel::Complex(2.4 + 0.01 * row, 0.05);
      du[p] = (row + 1 < n) ? channel::Complex(-0.11, 0.03) : channel::Complex(0.0, 0.0);
      dw[p] = (row + 2 < n) ? channel::Complex(-0.02, -0.01) : channel::Complex(0.0, 0.0);
    }
  }

  auto rhs = apply_penta(ds, dl, d, du, dw, exact, n, batch);
  channel::YLineWorkspace ws;
  ws.resize(n, batch);
  ws.copy_from_host(ds, dl, d, du, dw, rhs);
  ws.solve("test_pentadiagonal");
  auto solved = ws.x_host();
  for (std::size_t i = 0; i < solved.size(); ++i) {
    channel::test::require_near(solved[i], exact[i], 1.0e-11, "pentadiagonal solution");
  }
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  run_case(1, 5);
  run_case(2, 4);
  run_case(7, 6);
  if (runtime.rank() == 0) std::cout << "Pentadiagonal batch test PASSED\n";
  return 0;
}
