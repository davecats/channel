#include "channel/mean_correction.hpp"
#include "channel/dns_solver.hpp"
#include "channel/runtime.hpp"

#include "test_common.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <vector>

namespace {

int line_index(int y) {
  return y + 1;
}

int coeff_index(int row, int offset) {
  return row * 5 + (offset + 2);
}

double magnitude(const channel::Complex& value) {
  return std::hypot(static_cast<double>(value.real()), static_cast<double>(value.imag()));
}

std::vector<channel::Complex> dense_solve(std::vector<channel::Complex> a,
                                          std::vector<channel::Complex> b,
                                          int n) {
  for (int col = 0; col < n; ++col) {
    int pivot = col;
    double pivot_norm = magnitude(a[static_cast<std::size_t>(col * n + col)]);
    for (int row = col + 1; row < n; ++row) {
      const double norm = magnitude(a[static_cast<std::size_t>(row * n + col)]);
      if (norm > pivot_norm) {
        pivot = row;
        pivot_norm = norm;
      }
    }
    channel::test::require(pivot_norm > 0.0, "dense compact reference is nonsingular");
    if (pivot != col) {
      for (int j = col; j < n; ++j) {
        std::swap(a[static_cast<std::size_t>(col * n + j)], a[static_cast<std::size_t>(pivot * n + j)]);
      }
      std::swap(b[static_cast<std::size_t>(col)], b[static_cast<std::size_t>(pivot)]);
    }

    const auto diag = a[static_cast<std::size_t>(col * n + col)];
    for (int row = col + 1; row < n; ++row) {
      const auto factor = a[static_cast<std::size_t>(row * n + col)] / diag;
      a[static_cast<std::size_t>(row * n + col)] = channel::Complex(0.0, 0.0);
      for (int j = col + 1; j < n; ++j) {
        a[static_cast<std::size_t>(row * n + j)] -= factor * a[static_cast<std::size_t>(col * n + j)];
      }
      b[static_cast<std::size_t>(row)] -= factor * b[static_cast<std::size_t>(col)];
    }
  }

  std::vector<channel::Complex> x(static_cast<std::size_t>(n));
  for (int row = n - 1; row >= 0; --row) {
    auto value = b[static_cast<std::size_t>(row)];
    for (int col = row + 1; col < n; ++col) {
      value -= a[static_cast<std::size_t>(row * n + col)] * x[static_cast<std::size_t>(col)];
    }
    x[static_cast<std::size_t>(row)] = value / a[static_cast<std::size_t>(row * n + row)];
  }
  return x;
}

std::vector<channel::Complex> compact_reference(int ny,
                                                const std::vector<double>& matrix,
                                                const std::vector<channel::Complex>& rhs_line,
                                                const channel::CompactBoundaryRows& boundaries) {
  const int n = ny + 3;
  std::vector<channel::Complex> a(static_cast<std::size_t>(n * n), channel::Complex(0.0, 0.0));
  std::vector<channel::Complex> b(static_cast<std::size_t>(n), channel::Complex(0.0, 0.0));

  auto add_centered_row = [&](int eq, int center, const std::array<double, 5>& row) {
    for (int offset = -2; offset <= 2; ++offset) {
      const int y = center + offset;
      a[static_cast<std::size_t>(eq * n + line_index(y))] += row[static_cast<std::size_t>(offset + 2)];
    }
  };

  add_centered_row(0, 1, boundaries.lower_ghost);
  b[0] = boundaries.rhs_lower_ghost;
  add_centered_row(1, 1, boundaries.lower);
  b[1] = boundaries.rhs_lower;

  for (int iy = 1; iy <= ny - 1; ++iy) {
    const int eq = iy + 1;
    for (int offset = -2; offset <= 2; ++offset) {
      a[static_cast<std::size_t>(eq * n + line_index(iy + offset))] +=
          matrix[static_cast<std::size_t>(coeff_index(iy, offset))];
    }
    b[static_cast<std::size_t>(eq)] = rhs_line[static_cast<std::size_t>(line_index(iy))];
  }

  add_centered_row(ny + 1, ny - 1, boundaries.upper);
  b[static_cast<std::size_t>(ny + 1)] = boundaries.rhs_upper;
  add_centered_row(ny + 2, ny - 1, boundaries.upper_ghost);
  b[static_cast<std::size_t>(ny + 2)] = boundaries.rhs_upper_ghost;

  return dense_solve(std::move(a), std::move(b), n);
}

void require_solver_storage_uses_default_memory(channel::FullLineCompactSolver& solver) {
  using LineView = decltype(solver.line().view());
  static_assert(std::is_same_v<typename LineView::memory_space, channel::DefaultMemorySpace>);
  channel::test::require(solver.line().data() != nullptr, "mean-correction line buffer is allocated");
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();

  constexpr int ny = 9;
  std::vector<double> matrix(static_cast<std::size_t>(ny + 1) * 5, 0.0);
  std::vector<channel::Complex> rhs_line(static_cast<std::size_t>(ny + 3), channel::Complex(0.0, 0.0));

  for (int iy = 1; iy <= ny - 1; ++iy) {
    matrix[static_cast<std::size_t>(coeff_index(iy, -2))] = -0.020 - 0.001 * iy;
    matrix[static_cast<std::size_t>(coeff_index(iy, -1))] = -0.110 + 0.002 * iy;
    matrix[static_cast<std::size_t>(coeff_index(iy, 0))] = 2.300 + 0.030 * iy;
    matrix[static_cast<std::size_t>(coeff_index(iy, 1))] = -0.085 - 0.001 * iy;
    matrix[static_cast<std::size_t>(coeff_index(iy, 2))] = -0.017 + 0.0005 * iy;
    rhs_line[static_cast<std::size_t>(line_index(iy))] =
        channel::Complex(1.0 + 0.07 * iy, -0.20 + 0.03 * iy);
  }

  channel::CompactBoundaryRows boundaries;
  boundaries.lower_ghost = {1.0, 0.0, 0.0, 0.0, 0.0};
  boundaries.lower = {0.0, 1.0, 0.0, 0.0, 0.0};
  boundaries.upper = {0.0, 0.0, 0.0, 1.0, 0.0};
  boundaries.upper_ghost = {0.0, 0.0, 0.0, 0.0, 1.0};
  boundaries.rhs_lower_ghost = channel::Complex(0.3, -0.1);
  boundaries.rhs_lower = channel::Complex(-0.2, 0.05);
  boundaries.rhs_upper = channel::Complex(0.4, 0.2);
  boundaries.rhs_upper_ghost = channel::Complex(-0.1, 0.3);

  const auto expected = compact_reference(ny, matrix, rhs_line, boundaries);

  channel::FullLineCompactSolver solver;
  solver.resize(ny);
  require_solver_storage_uses_default_memory(solver);
  if (runtime.rank() == 0) {
    solver.copy_system_from_host(matrix, rhs_line);
  } else {
    std::vector<double> unused_matrix(matrix.size(), 0.0);
    std::vector<channel::Complex> sentinel_line(rhs_line.size(), channel::Complex(99.0, -99.0));
    solver.copy_system_from_host(unused_matrix, sentinel_line);
  }
  solver.solve_on_root_and_broadcast(boundaries, 0, channel::world_comm(), "test_mean_correction");

  const auto got = solver.line_host();
  channel::test::require(got.size() == expected.size(), "mean-correction full line size");
  for (std::size_t i = 0; i < got.size(); ++i) {
    channel::test::require_near(got[i], expected[i], 2.0e-12, "mean-correction full-line solve");
  }

  const auto local_range = channel::split_range(runtime.rank(), ny - 1, runtime.size());
  channel::test::require(local_range.count > 0, "mean-correction DNS test has local active rows");
  channel::DnsState state;
  state.resize({2, local_range.count, 1, 2});
  std::vector<channel::Complex> component0(state.values_per_component(), channel::Complex(-7.0, 7.0));
  std::vector<channel::Complex> component1(state.values_per_component(), channel::Complex(3.0, -4.0));
  const int line_count = static_cast<int>(state.line_count());
  for (int local_y = 0; local_y < local_range.count; ++local_y) {
    const int global_y = local_range.first + local_y + 1;
    component0[static_cast<std::size_t>(local_y * line_count)] =
        rhs_line[static_cast<std::size_t>(line_index(global_y))];
  }
  state.copy_component_from_host(0, component0);
  state.copy_component_from_host(1, component1);

  channel::DnsMeanCorrectionConfig cfg;
  cfg.component = 0;
  cfg.line = 0;
  cfg.root = 0;
  cfg.comm_y = channel::world_comm();
  cfg.boundaries = boundaries;
  cfg.label = "test_dns_mean_correction";
  channel::DnsMeanCorrectionStep dns_step;
  dns_step.prepare(state, cfg);
  dns_step.copy_matrix_from_host(matrix);
  dns_step.apply(state);

  const auto dns_component0 = state.component_host(0);
  const auto dns_component1 = state.component_host(1);
  for (int local_y = 0; local_y < local_range.count; ++local_y) {
    const int global_y = local_range.first + local_y + 1;
    const auto solved = dns_component0[static_cast<std::size_t>(local_y * line_count)];
    channel::test::require_near(solved, expected[static_cast<std::size_t>(line_index(global_y))], 2.0e-12,
                                "DNS mean-correction selected line");
    const auto untouched_line = dns_component0[static_cast<std::size_t>(local_y * line_count + 1)];
    channel::test::require_near(untouched_line, component0[static_cast<std::size_t>(local_y * line_count + 1)], 0.0,
                                "DNS mean-correction preserves neighboring line");
  }
  for (std::size_t i = 0; i < dns_component1.size(); ++i) {
    channel::test::require_near(dns_component1[i], component1[i], 0.0,
                                "DNS mean-correction preserves inactive component");
  }
  channel::test::require(dns_step.global_ny() == ny, "DNS mean-correction infers global ny");
  channel::test::require(dns_step.first_global_y() == local_range.first + 1,
                         "DNS mean-correction infers first global y");

  channel::DnsState linear_state;
  linear_state.resize({2, local_range.count, 1, 2});
  linear_state.copy_component_from_host(0, component0);
  linear_state.copy_component_from_host(1, component1);
  channel::DnsLinearStepConfig linear_cfg;
  linear_cfg.mean_corrections = {cfg};
  channel::DnsLinearStep linear_step;
  linear_step.prepare(linear_state, linear_cfg);
  linear_step.copy_mean_correction_matrix_from_host(0, matrix);
  linear_step.advance(linear_state);
  channel::test::require(linear_step.mean_correction_count() == 1,
                         "DNS linear step owns one mean-correction stage");
  const auto linear_component0 = linear_state.component_host(0);
  for (int local_y = 0; local_y < local_range.count; ++local_y) {
    const int global_y = local_range.first + local_y + 1;
    const auto solved = linear_component0[static_cast<std::size_t>(local_y * line_count)];
    channel::test::require_near(solved, expected[static_cast<std::size_t>(line_index(global_y))], 2.0e-12,
                                "DNS linear-step mean-correction selected line");
  }

  if (runtime.rank() == 0) std::cout << "Mean correction compact line test PASSED\n";
  return 0;
}
