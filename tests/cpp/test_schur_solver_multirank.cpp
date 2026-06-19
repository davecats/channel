#include "channel/runtime.hpp"
#include "channel/schur.hpp"

#include "test_common.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

namespace {

constexpr int row_width = channel::SchurLevel::row_width;
constexpr int value_width = channel::SchurLevel::value_width;
constexpr int max_rows = channel::SchurLevel::max_rows;
constexpr int bandwidth = channel::SchurLevel::bandwidth;

int band_index(int row, int col) {
  return row * (2 * bandwidth + 1) + col;
}

int rhs_index(int row, int rhs) {
  return row * 5 + rhs;
}

int recovery_index(int line, int rhs, int row, int nrows) {
  return line * 5 * nrows + rhs * nrows + row;
}

int line_row(int line, int row) {
  return line * row_width + row;
}

int line_value(int line, int value) {
  return line * value_width + value;
}

channel::Complex leaf_row_value(int child, int line, int k) {
  const double base = 0.08 * static_cast<double>(child + 1) + 0.01 * static_cast<double>(line + 1);
  if (k < 4) {
    return channel::Complex(1.0 + base + 0.02 * static_cast<double>(k + 1),
                            0.003 * static_cast<double>((child + 1) * (k + 1)));
  }
  return channel::Complex(0.015 + 0.002 * static_cast<double>(k) + 0.001 * static_cast<double>(line),
                          -0.0005 * static_cast<double>((child + 1) * (k + 1)));
}

void factor_banded(std::array<channel::Complex, max_rows * (2 * bandwidth + 1)>& a, int n) {
  for (int i = 0; i < n; ++i) {
    const auto piv = a[static_cast<std::size_t>(band_index(i, bandwidth))];
    const int last = std::min(bandwidth, n - i - 1);
    for (int j = 1; j <= last; ++j) {
      const auto factor = a[static_cast<std::size_t>(band_index(i + j, bandwidth - j))] / piv;
      a[static_cast<std::size_t>(band_index(i + j, bandwidth - j))] = factor;
      for (int t = 1; t <= last; ++t) {
        a[static_cast<std::size_t>(band_index(i + j, bandwidth + t - j))] -=
            factor * a[static_cast<std::size_t>(band_index(i, bandwidth + t))];
      }
    }
  }
}

void solve_banded(std::array<channel::Complex, max_rows * 5>& rhs,
                  const std::array<channel::Complex, max_rows * (2 * bandwidth + 1)>& a,
                  int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = std::max(0, i - bandwidth); j < i; ++j) {
      const auto factor = a[static_cast<std::size_t>(band_index(i, bandwidth + j - i))];
      for (int irhs = 0; irhs < 5; ++irhs) {
        rhs[static_cast<std::size_t>(rhs_index(i, irhs))] -=
            factor * rhs[static_cast<std::size_t>(rhs_index(j, irhs))];
      }
    }
  }

  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j <= std::min(n - 1, i + bandwidth); ++j) {
      const auto factor = a[static_cast<std::size_t>(band_index(i, bandwidth + j - i))];
      for (int irhs = 0; irhs < 5; ++irhs) {
        rhs[static_cast<std::size_t>(rhs_index(i, irhs))] -=
            factor * rhs[static_cast<std::size_t>(rhs_index(j, irhs))];
      }
    }
    const auto piv = a[static_cast<std::size_t>(band_index(i, bandwidth))];
    for (int irhs = 0; irhs < 5; ++irhs) {
      rhs[static_cast<std::size_t>(rhs_index(i, irhs))] /= piv;
    }
  }
}

std::vector<channel::Complex> expected_leaf_values(int child_id, int arity, int nlines) {
  const int nrows = 4 * arity;
  std::vector<channel::Complex> expected(static_cast<std::size_t>(nlines * value_width),
                                         channel::Complex(0.0, 0.0));

  for (int owner = 0; owner < arity; ++owner) {
    const auto owned = channel::split_range(owner, nlines, arity);
    for (int local_line = 0; local_line < owned.count; ++local_line) {
      const int line = owned.first + local_line;
      std::array<channel::Complex, max_rows * (2 * bandwidth + 1)> a{};
      std::array<channel::Complex, max_rows * 5> rhs{};

      for (int child = 0; child < arity; ++child) {
        const int row0 = 4 * child;
        for (int k = 0; k < 4; ++k) {
          const int row = row0 + k;
          a[static_cast<std::size_t>(band_index(row, bandwidth))] = channel::Complex(1.0, 0.0);
          rhs[static_cast<std::size_t>(rhs_index(row, 0))] = leaf_row_value(child, line, k);
          if (child > 0) {
            int col = row0 - 2;
            a[static_cast<std::size_t>(band_index(row, bandwidth + col - row))] = -leaf_row_value(child, line, 4 + k);
            col = row0 - 1;
            a[static_cast<std::size_t>(band_index(row, bandwidth + col - row))] = -leaf_row_value(child, line, 8 + k);
          } else {
            rhs[static_cast<std::size_t>(rhs_index(row, 1))] = leaf_row_value(child, line, 4 + k);
            rhs[static_cast<std::size_t>(rhs_index(row, 2))] = leaf_row_value(child, line, 8 + k);
          }
          if (child < arity - 1) {
            int col = row0 + 4;
            a[static_cast<std::size_t>(band_index(row, bandwidth + col - row))] = -leaf_row_value(child, line, 12 + k);
            col = row0 + 5;
            a[static_cast<std::size_t>(band_index(row, bandwidth + col - row))] = -leaf_row_value(child, line, 16 + k);
          } else {
            rhs[static_cast<std::size_t>(rhs_index(row, 3))] = leaf_row_value(child, line, 12 + k);
            rhs[static_cast<std::size_t>(rhs_index(row, 4))] = leaf_row_value(child, line, 16 + k);
          }
        }
      }

      factor_banded(a, nrows);
      solve_banded(rhs, a, nrows);

      const int row0 = 4 * child_id;
      expected[static_cast<std::size_t>(line_value(line, 0))] = rhs[static_cast<std::size_t>(rhs_index(row0 + 0, 0))];
      expected[static_cast<std::size_t>(line_value(line, 1))] = rhs[static_cast<std::size_t>(rhs_index(row0 + 1, 0))];
      expected[static_cast<std::size_t>(line_value(line, 2))] = rhs[static_cast<std::size_t>(rhs_index(row0 + 2, 0))];
      expected[static_cast<std::size_t>(line_value(line, 3))] = rhs[static_cast<std::size_t>(rhs_index(row0 + 3, 0))];
      expected[static_cast<std::size_t>(line_value(line, 4))] =
          child_id > 0 ? rhs[static_cast<std::size_t>(rhs_index(row0 - 2, 0))]
                       : channel::Complex(0.0, 0.0);
      expected[static_cast<std::size_t>(line_value(line, 5))] =
          child_id > 0 ? rhs[static_cast<std::size_t>(rhs_index(row0 - 1, 0))]
                       : channel::Complex(0.0, 0.0);
      expected[static_cast<std::size_t>(line_value(line, 6))] =
          child_id < arity - 1 ? rhs[static_cast<std::size_t>(rhs_index(row0 + 4, 0))]
                               : channel::Complex(0.0, 0.0);
      expected[static_cast<std::size_t>(line_value(line, 7))] =
          child_id < arity - 1 ? rhs[static_cast<std::size_t>(rhs_index(row0 + 5, 0))]
                               : channel::Complex(0.0, 0.0);
    }
  }

  return expected;
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();
  channel::test::require(runtime.size() == 2, "Schur solver multirank test must run with two MPI ranks");

  constexpr int nlines = 5;
  std::vector<channel::Complex> leaf_rows(static_cast<std::size_t>(nlines * row_width));
  for (int line = 0; line < nlines; ++line) {
    for (int k = 0; k < row_width; ++k) {
      leaf_rows[static_cast<std::size_t>(line_row(line, k))] = leaf_row_value(runtime.rank(), line, k);
    }
  }

  channel::DeviceVector<channel::Complex> leaf_rows_device("schur_leaf_rows", leaf_rows.size());
  channel::DeviceVector<channel::Complex> leaf_values_device("schur_leaf_values", nlines * value_width);
  leaf_rows_device.copy_from_host(leaf_rows);

  channel::SchurSolver solver;
  channel::SchurSolverConfig cfg;
  cfg.npy = runtime.size();
  cfg.ipy = runtime.rank();
  cfg.nlines = nlines;
  cfg.pass_counts = {runtime.size()};
  cfg.exchange_mode = channel::ExchangeMode::AllToAll;
  cfg.comm_y = channel::world_comm();
  solver.prepare(cfg);
  solver.solve_from_leaf_rows(leaf_rows_device, leaf_values_device);

  const auto actual = leaf_values_device.copy_to_host();
  const auto expected = expected_leaf_values(runtime.rank(), runtime.size(), nlines);
  for (std::size_t i = 0; i < actual.size(); ++i) {
    channel::test::require_near(actual[i], expected[i], 5.0e-12,
                                "Schur solver multirank leaf value rank " + std::to_string(runtime.rank()) +
                                    " index " + std::to_string(i));
  }

  if (runtime.rank() == 0) std::cout << "Schur solver multirank test PASSED\n";
  return 0;
}
