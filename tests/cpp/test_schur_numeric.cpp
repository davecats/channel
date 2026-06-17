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

void factor_banded(std::array<channel::Complex, max_rows * (2 * bandwidth + 1)>& a, int n) {
  for (int i = 0; i < n; ++i) {
    const auto piv = a[band_index(i, bandwidth)];
    const int last = std::min(bandwidth, n - i - 1);
    for (int j = 1; j <= last; ++j) {
      const auto factor = a[band_index(i + j, bandwidth - j)] / piv;
      a[band_index(i + j, bandwidth - j)] = factor;
      for (int t = 1; t <= last; ++t) {
        a[band_index(i + j, bandwidth + t - j)] -= factor * a[band_index(i, bandwidth + t)];
      }
    }
  }
}

void solve_banded(std::array<channel::Complex, max_rows * 5>& rhs,
                  const std::array<channel::Complex, max_rows * (2 * bandwidth + 1)>& a,
                  int n) {
  for (int i = 0; i < n; ++i) {
    for (int j = std::max(0, i - bandwidth); j < i; ++j) {
      const auto factor = a[band_index(i, bandwidth + j - i)];
      for (int irhs = 0; irhs < 5; ++irhs) {
        rhs[rhs_index(i, irhs)] -= factor * rhs[rhs_index(j, irhs)];
      }
    }
  }

  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j <= std::min(n - 1, i + bandwidth); ++j) {
      const auto factor = a[band_index(i, bandwidth + j - i)];
      for (int irhs = 0; irhs < 5; ++irhs) {
        rhs[rhs_index(i, irhs)] -= factor * rhs[rhs_index(j, irhs)];
      }
    }
    const auto piv = a[band_index(i, bandwidth)];
    for (int irhs = 0; irhs < 5; ++irhs) {
      rhs[rhs_index(i, irhs)] /= piv;
    }
  }
}

void reference_compose(const channel::SchurLevelConfig& cfg,
                       const std::vector<channel::Complex>& exchanged,
                       std::vector<channel::Complex>& rows,
                       std::vector<channel::Complex>& basis) {
  const int nrows = 4 * cfg.arity;
  const int nsolve = cfg.owned.count;
  rows.assign(static_cast<std::size_t>(nsolve * row_width), channel::Complex(0.0, 0.0));
  basis.assign(static_cast<std::size_t>(nsolve * 5 * nrows), channel::Complex(0.0, 0.0));

  for (int line = 0; line < nsolve; ++line) {
    std::array<channel::Complex, max_rows * (2 * bandwidth + 1)> a{};
    std::array<channel::Complex, max_rows * 5> rhs{};

    for (int child = 0; child < cfg.arity; ++child) {
      const int offset = child * row_width * cfg.owned.count + line * row_width;
      const int row0 = 4 * child;
      for (int k = 0; k < 4; ++k) {
        const int row = row0 + k;
        a[band_index(row, bandwidth)] = channel::Complex(1.0, 0.0);
        rhs[rhs_index(row, 0)] = exchanged[static_cast<std::size_t>(offset + k)];
        if (child > 0) {
          int col = row0 - 2;
          a[band_index(row, bandwidth + col - row)] = -exchanged[static_cast<std::size_t>(offset + 4 + k)];
          col = row0 - 1;
          a[band_index(row, bandwidth + col - row)] = -exchanged[static_cast<std::size_t>(offset + 8 + k)];
        } else {
          rhs[rhs_index(row, 1)] = exchanged[static_cast<std::size_t>(offset + 4 + k)];
          rhs[rhs_index(row, 2)] = exchanged[static_cast<std::size_t>(offset + 8 + k)];
        }
        if (child < cfg.arity - 1) {
          int col = row0 + 4;
          a[band_index(row, bandwidth + col - row)] = -exchanged[static_cast<std::size_t>(offset + 12 + k)];
          col = row0 + 5;
          a[band_index(row, bandwidth + col - row)] = -exchanged[static_cast<std::size_t>(offset + 16 + k)];
        } else {
          rhs[rhs_index(row, 3)] = exchanged[static_cast<std::size_t>(offset + 12 + k)];
          rhs[rhs_index(row, 4)] = exchanged[static_cast<std::size_t>(offset + 16 + k)];
        }
      }
    }

    factor_banded(a, nrows);
    solve_banded(rhs, a, nrows);

    for (int irhs = 0; irhs < 5; ++irhs) {
      for (int row = 0; row < nrows; ++row) {
        basis[static_cast<std::size_t>(recovery_index(line, irhs, row, nrows))] = rhs[rhs_index(row, irhs)];
      }
    }
    for (int k = 0; k < 4; ++k) {
      const int exposed = k < 2 ? k : 4 * (cfg.arity - 1) + k;
      rows[static_cast<std::size_t>(line_row(line, k))] = rhs[rhs_index(exposed, 0)];
      for (int ext = 0; ext < 4; ++ext) {
        rows[static_cast<std::size_t>(line_row(line, 4 * (ext + 1) + k))] = rhs[rhs_index(exposed, ext + 1)];
      }
    }
  }
}

std::vector<channel::Complex> reference_pack(const channel::SchurLevelConfig& cfg,
                                             const std::vector<channel::Complex>& values,
                                             const std::vector<channel::Complex>& basis) {
  const int nrows = 4 * cfg.arity;
  const int nsolve = cfg.owned.count;
  std::vector<channel::Complex> packed(static_cast<std::size_t>(cfg.arity * value_width * nsolve),
                                       channel::Complex(0.0, 0.0));

  for (int line = 0; line < nsolve; ++line) {
    std::array<channel::Complex, max_rows> recovered{};
    const auto prev1 = values[static_cast<std::size_t>(line_value(line, 4))];
    const auto prev2 = values[static_cast<std::size_t>(line_value(line, 5))];
    const auto next1 = values[static_cast<std::size_t>(line_value(line, 6))];
    const auto next2 = values[static_cast<std::size_t>(line_value(line, 7))];
    for (int row = 0; row < nrows; ++row) {
      recovered[static_cast<std::size_t>(row)] =
          basis[static_cast<std::size_t>(recovery_index(line, 0, row, nrows))] +
          basis[static_cast<std::size_t>(recovery_index(line, 1, row, nrows))] * prev1 +
          basis[static_cast<std::size_t>(recovery_index(line, 2, row, nrows))] * prev2 +
          basis[static_cast<std::size_t>(recovery_index(line, 3, row, nrows))] * next1 +
          basis[static_cast<std::size_t>(recovery_index(line, 4, row, nrows))] * next2;
    }

    for (int child = 0; child < cfg.arity; ++child) {
      const int row0 = 4 * child;
      const int offset = child * value_width * nsolve + line * value_width;
      packed[static_cast<std::size_t>(offset + 0)] = recovered[static_cast<std::size_t>(row0 + 0)];
      packed[static_cast<std::size_t>(offset + 1)] = recovered[static_cast<std::size_t>(row0 + 1)];
      packed[static_cast<std::size_t>(offset + 2)] = recovered[static_cast<std::size_t>(row0 + 2)];
      packed[static_cast<std::size_t>(offset + 3)] = recovered[static_cast<std::size_t>(row0 + 3)];
      packed[static_cast<std::size_t>(offset + 4)] = child > 0 ? recovered[static_cast<std::size_t>(row0 - 2)] : prev1;
      packed[static_cast<std::size_t>(offset + 5)] = child > 0 ? recovered[static_cast<std::size_t>(row0 - 1)] : prev2;
      packed[static_cast<std::size_t>(offset + 6)] =
          child < cfg.arity - 1 ? recovered[static_cast<std::size_t>(row0 + 4)] : next1;
      packed[static_cast<std::size_t>(offset + 7)] =
          child < cfg.arity - 1 ? recovered[static_cast<std::size_t>(row0 + 5)] : next2;
    }
  }

  return packed;
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  channel::test::require_kokkos_runtime();

  channel::SchurLevelConfig cfg;
  cfg.level_index = 0;
  cfg.arity = 3;
  cfg.child_id = 1;
  cfg.prev_first = 0;
  cfg.prev_count = 6;
  cfg.owned = {2, 2};
  cfg.exchange_mode = channel::ExchangeMode::AllToAll;
  cfg.comm = channel::world_comm();

  std::vector<channel::Complex> exchanged(static_cast<std::size_t>(cfg.arity * row_width * cfg.owned.count));
  for (int child = 0; child < cfg.arity; ++child) {
    for (int line = 0; line < cfg.owned.count; ++line) {
      const int offset = child * row_width * cfg.owned.count + line * row_width;
      for (int k = 0; k < 4; ++k) {
        exchanged[static_cast<std::size_t>(offset + k)] =
            channel::Complex(0.7 + 0.1 * child + 0.03 * line + 0.01 * k, -0.02 * k);
        exchanged[static_cast<std::size_t>(offset + 4 + k)] =
            channel::Complex(0.015 + 0.002 * child + 0.001 * k, 0.001);
        exchanged[static_cast<std::size_t>(offset + 8 + k)] =
            channel::Complex(-0.011 - 0.001 * child + 0.001 * k, -0.002);
        exchanged[static_cast<std::size_t>(offset + 12 + k)] =
            channel::Complex(0.012 + 0.001 * child - 0.001 * k, 0.003);
        exchanged[static_cast<std::size_t>(offset + 16 + k)] =
            channel::Complex(-0.009 + 0.001 * child + 0.0005 * k, -0.001);
      }
    }
  }

  std::vector<channel::Complex> ref_rows;
  std::vector<channel::Complex> ref_basis;
  reference_compose(cfg, exchanged, ref_rows, ref_basis);

  channel::DeviceVector<channel::Complex> exchanged_device("schur_numeric_exchanged", exchanged.size());
  exchanged_device.copy_from_host(exchanged);
  channel::SchurLevel level(cfg);
  level.compose_from_exchanged_rows(exchanged_device);

  const auto rows = level.rows_host();
  const auto basis = level.recovery_basis_host();
  for (std::size_t i = 0; i < ref_rows.size(); ++i) {
    channel::test::require_near(rows[i], ref_rows[i], 1.0e-12, "composed Schur row");
  }
  for (std::size_t i = 0; i < ref_basis.size(); ++i) {
    channel::test::require_near(basis[i], ref_basis[i], 1.0e-12, "Schur recovery basis");
  }

  level.seed_root_values();
  channel::DeviceVector<channel::Complex> packed("schur_numeric_packed", cfg.arity * value_width * cfg.owned.count);
  level.pack_recovered_values(packed);
  const auto values = level.values_host();
  const auto ref_packed = reference_pack(cfg, values, ref_basis);
  const auto got_packed = packed.copy_to_host();
  for (std::size_t i = 0; i < ref_packed.size(); ++i) {
    channel::test::require_near(got_packed[i], ref_packed[i], 1.0e-12, "packed recovered Schur value");
  }

  if (runtime.rank() == 0) std::cout << "Schur numeric compose/recover test PASSED\n";
  return 0;
}
