#include "channel/schur.hpp"
#include "channel/runtime.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace channel {

namespace {

static_assert(SchurLevel::bandwidth == 5, "update schur_bandwidth() if SchurLevel::bandwidth changes");
static_assert(SchurLevel::row_width == 20, "update schur_row_width() if SchurLevel::row_width changes");
static_assert(SchurLevel::value_width == 8, "update schur_value_width() if SchurLevel::value_width changes");

KOKKOS_INLINE_FUNCTION int schur_bandwidth() {
  return 5;
}

KOKKOS_INLINE_FUNCTION int schur_row_width() {
  return 20;
}

KOKKOS_INLINE_FUNCTION int schur_value_width() {
  return 8;
}

KOKKOS_INLINE_FUNCTION int band_index(int row, int col) {
  return row * (2 * schur_bandwidth() + 1) + col;
}

KOKKOS_INLINE_FUNCTION int rhs_index(int row, int rhs) {
  return row * 5 + rhs;
}

KOKKOS_INLINE_FUNCTION int row_index(int line, int row) {
  return line * schur_row_width() + row;
}

KOKKOS_INLINE_FUNCTION int value_index(int line, int value) {
  return line * schur_value_width() + value;
}

KOKKOS_INLINE_FUNCTION int recovery_index(int line, int rhs, int row, int nrows) {
  return line * 5 * nrows + rhs * nrows + row;
}

KOKKOS_INLINE_FUNCTION void factor_banded_fixed(Complex* a, int n) {
  for (int i = 0; i < n; ++i) {
    const Complex piv = a[band_index(i, schur_bandwidth())];
    const int last = Kokkos::min(schur_bandwidth(), n - i - 1);
    for (int j = 1; j <= last; ++j) {
      const Complex factor = a[band_index(i + j, schur_bandwidth() - j)] / piv;
      a[band_index(i + j, schur_bandwidth() - j)] = factor;
      for (int t = 1; t <= last; ++t) {
        a[band_index(i + j, schur_bandwidth() + t - j)] -=
            factor * a[band_index(i, schur_bandwidth() + t)];
      }
    }
  }
}

KOKKOS_INLINE_FUNCTION void solve_factored_banded_multi_fixed(Complex* rhs, const Complex* a, int n, int nrhs) {
  for (int i = 0; i < n; ++i) {
    for (int j = Kokkos::max(0, i - schur_bandwidth()); j < i; ++j) {
      const Complex factor = a[band_index(i, schur_bandwidth() + j - i)];
      for (int irhs = 0; irhs < nrhs; ++irhs) {
        rhs[rhs_index(i, irhs)] -= factor * rhs[rhs_index(j, irhs)];
      }
    }
  }

  for (int i = n - 1; i >= 0; --i) {
    for (int j = i + 1; j <= Kokkos::min(n - 1, i + schur_bandwidth()); ++j) {
      const Complex factor = a[band_index(i, schur_bandwidth() + j - i)];
      for (int irhs = 0; irhs < nrhs; ++irhs) {
        rhs[rhs_index(i, irhs)] -= factor * rhs[rhs_index(j, irhs)];
      }
    }
    const Complex piv = a[band_index(i, schur_bandwidth())];
    for (int irhs = 0; irhs < nrhs; ++irhs) {
      rhs[rhs_index(i, irhs)] /= piv;
    }
  }
}

ExchangeMode resolved_exchange_mode(ExchangeMode requested, int level_index, int pass_count, int arity) {
  if (requested != ExchangeMode::Auto) {
    if (requested == ExchangeMode::AllGather && level_index + 1 != pass_count) {
      return ExchangeMode::AllToAll;
    }
    return requested;
  }
  if (level_index + 1 == pass_count && arity == 2) {
    return ExchangeMode::AllGather;
  }
  return ExchangeMode::AllToAll;
}

} // namespace

SchurLevel::SchurLevel(SchurLevelConfig cfg) : cfg_(std::move(cfg)) {
  resize_buffers();
}

SchurLevel::SchurLevel(SchurLevel&& other) noexcept
    : cfg_(other.cfg_),
      rows_(std::move(other.rows_)),
      values_(std::move(other.values_)),
      recovery_basis_(std::move(other.recovery_basis_)),
      row_send_(std::move(other.row_send_)),
      row_recv_(std::move(other.row_recv_)),
      value_send_(std::move(other.value_send_)),
      value_recv_(std::move(other.value_recv_)),
      child_ranges_(std::move(other.child_ranges_)),
      row_send_counts_(std::move(other.row_send_counts_)),
      row_send_displs_(std::move(other.row_send_displs_)),
      row_recv_counts_(std::move(other.row_recv_counts_)),
      row_recv_displs_(std::move(other.row_recv_displs_)),
      value_send_counts_(std::move(other.value_send_counts_)),
      value_send_displs_(std::move(other.value_send_displs_)),
      value_recv_counts_(std::move(other.value_recv_counts_)),
      value_recv_displs_(std::move(other.value_recv_displs_)) {
  other.cfg_.comm = world_comm();
}

SchurLevel& SchurLevel::operator=(SchurLevel&& other) noexcept {
  if (this != &other) {
    release_comm();
    cfg_ = other.cfg_;
    rows_ = std::move(other.rows_);
    values_ = std::move(other.values_);
    recovery_basis_ = std::move(other.recovery_basis_);
    row_send_ = std::move(other.row_send_);
    row_recv_ = std::move(other.row_recv_);
    value_send_ = std::move(other.value_send_);
    value_recv_ = std::move(other.value_recv_);
    child_ranges_ = std::move(other.child_ranges_);
    row_send_counts_ = std::move(other.row_send_counts_);
    row_send_displs_ = std::move(other.row_send_displs_);
    row_recv_counts_ = std::move(other.row_recv_counts_);
    row_recv_displs_ = std::move(other.row_recv_displs_);
    value_send_counts_ = std::move(other.value_send_counts_);
    value_send_displs_ = std::move(other.value_send_displs_);
    value_recv_counts_ = std::move(other.value_recv_counts_);
    value_recv_displs_ = std::move(other.value_recv_displs_);
    other.cfg_.comm = world_comm();
  }
  return *this;
}

SchurLevel::~SchurLevel() {
  release_comm();
}

void SchurLevel::resize_buffers() {
  if (cfg_.arity < 1 || cfg_.arity > max_arity) {
    throw std::runtime_error("SchurLevel arity out of range");
  }
  if (cfg_.prev_count < 1 || cfg_.owned.count < 0) {
    throw std::runtime_error("SchurLevel requires a non-empty previous line range and non-negative owned range");
  }
  const int max_owned = std::max(cfg_.prev_count, cfg_.owned.count);
  rows_.resize("schur_rows", static_cast<std::size_t>(row_width) * static_cast<std::size_t>(max_owned));
  values_.resize("schur_values", static_cast<std::size_t>(value_width) * static_cast<std::size_t>(max_owned));
  recovery_basis_.resize("schur_recovery_basis",
                         static_cast<std::size_t>(4 * cfg_.arity) * 5U * static_cast<std::size_t>(max_owned));

  child_ranges_.resize(static_cast<std::size_t>(cfg_.arity));
  row_send_counts_.resize(static_cast<std::size_t>(cfg_.arity));
  row_send_displs_.resize(static_cast<std::size_t>(cfg_.arity));
  row_recv_counts_.resize(static_cast<std::size_t>(cfg_.arity));
  row_recv_displs_.resize(static_cast<std::size_t>(cfg_.arity));
  value_send_counts_.resize(static_cast<std::size_t>(cfg_.arity));
  value_send_displs_.resize(static_cast<std::size_t>(cfg_.arity));
  value_recv_counts_.resize(static_cast<std::size_t>(cfg_.arity));
  value_recv_displs_.resize(static_cast<std::size_t>(cfg_.arity));

  for (int peer = 0; peer < cfg_.arity; ++peer) {
    child_ranges_[static_cast<std::size_t>(peer)] = split_range(peer, cfg_.prev_count, cfg_.arity);
    const auto peer_range = child_ranges_[static_cast<std::size_t>(peer)];
    row_send_counts_[static_cast<std::size_t>(peer)] = row_width * peer_range.count;
    row_send_displs_[static_cast<std::size_t>(peer)] = row_width * peer_range.first;
    row_recv_counts_[static_cast<std::size_t>(peer)] = row_width * cfg_.owned.count;
    row_recv_displs_[static_cast<std::size_t>(peer)] = peer * row_width * cfg_.owned.count;
    value_send_counts_[static_cast<std::size_t>(peer)] = value_width * cfg_.owned.count;
    value_send_displs_[static_cast<std::size_t>(peer)] = peer * value_width * cfg_.owned.count;
    value_recv_counts_[static_cast<std::size_t>(peer)] = value_width * peer_range.count;
    value_recv_displs_[static_cast<std::size_t>(peer)] = value_width * peer_range.first;
  }

  row_send_.resize("schur_row_send", static_cast<std::size_t>(row_width) * static_cast<std::size_t>(cfg_.prev_count));
  if (cfg_.exchange_mode == ExchangeMode::AllGather) {
    row_recv_.resize("schur_row_recv",
                     static_cast<std::size_t>(cfg_.arity) * static_cast<std::size_t>(row_width) *
                         static_cast<std::size_t>(cfg_.prev_count));
  } else {
    row_recv_.resize("schur_row_recv",
                     static_cast<std::size_t>(cfg_.arity) * static_cast<std::size_t>(row_width) *
                         static_cast<std::size_t>(cfg_.owned.count));
  }
  value_send_.resize("schur_value_send",
                     static_cast<std::size_t>(cfg_.arity) * static_cast<std::size_t>(value_width) *
                         static_cast<std::size_t>(cfg_.owned.count));
  value_recv_.resize("schur_value_recv",
                     static_cast<std::size_t>(value_width) * static_cast<std::size_t>(cfg_.prev_count));
}

int SchurLevel::solve_count(bool solve_redundant) const {
  return solve_redundant ? cfg_.prev_count : cfg_.owned.count;
}

void SchurLevel::release_comm() {
  if (cfg_.comm != world_comm()) {
    mpi_free(cfg_.comm);
    cfg_.comm = world_comm();
  }
}

void SchurLevel::pack_rows_for_exchange(const DeviceVector<Complex>& source_rows) {
  const auto expected = static_cast<std::size_t>(row_width) * static_cast<std::size_t>(cfg_.prev_count);
  if (source_rows.size() < expected) {
    throw std::runtime_error("SchurLevel::pack_rows_for_exchange source too small");
  }

  if (cfg_.exchange_mode == ExchangeMode::AllGather) {
    if (row_send_.size() != expected) row_send_.resize("schur_row_send", expected);
    auto source = source_rows.view();
    auto send = row_send_.view();
    Kokkos::parallel_for(
        "schur_pack_rows_allgather", Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(expected)),
        KOKKOS_LAMBDA(const int idx) { send(idx) = source(idx); });
    fence("schur_pack_rows_allgather");
    return;
  }

  auto source = source_rows.view();
  auto send = row_send_.view();
  for (int peer = 0; peer < cfg_.arity; ++peer) {
    const auto range = child_ranges_[static_cast<std::size_t>(peer)];
    const int displ = row_send_displs_[static_cast<std::size_t>(peer)];
    Kokkos::parallel_for(
        "schur_pack_rows", Kokkos::RangePolicy<ExecutionSpace>(0, row_width * range.count),
        KOKKOS_LAMBDA(const int idx) {
          const int line = idx / row_width;
          const int k = idx - line * row_width;
          send(displ + idx) = source(row_index(range.first + line, k));
        });
  }
  fence("schur_pack_rows");
}

void SchurLevel::exchange_rows() {
  if (cfg_.exchange_mode == ExchangeMode::AllGather) {
    allgather_complex_device(row_send_.data(), row_width * cfg_.prev_count, row_recv_.data(), cfg_.comm,
                             "schur_rows_allgather");
    return;
  }
  alltoallv_complex_device(row_send_.data(), row_send_counts_, row_send_displs_, row_recv_.data(), row_recv_counts_,
                           row_recv_displs_, cfg_.comm, "schur_rows_alltoallv");
}

void SchurLevel::compose_from_exchanged_rows(const DeviceVector<Complex>& exchanged_rows, bool solve_redundant) {
  const int nrows = reduced_row_count();
  const int nsolve = solve_count(solve_redundant);
  const int exchange_line_count =
      cfg_.exchange_mode == ExchangeMode::AllGather ? cfg_.prev_count : cfg_.owned.count;
  const auto expected = static_cast<std::size_t>(cfg_.arity) * static_cast<std::size_t>(row_width) *
                        static_cast<std::size_t>(exchange_line_count);
  if (exchanged_rows.size() < expected) {
    throw std::runtime_error("SchurLevel::compose_from_exchanged_rows input too small");
  }

  auto exchanged = exchanged_rows.view();
  auto rows = rows_.view();
  auto basis = recovery_basis_.view();
  const auto cfg = cfg_;
  Kokkos::parallel_for(
      "schur_compose_level", Kokkos::RangePolicy<ExecutionSpace>(0, nsolve), KOKKOS_LAMBDA(const int line) {
        Complex a[max_rows * (2 * bandwidth + 1)];
        Complex rhs[max_rows * 5];
        for (int i = 0; i < max_rows * (2 * bandwidth + 1); ++i) a[i] = Complex(0.0, 0.0);
        for (int i = 0; i < max_rows * 5; ++i) rhs[i] = Complex(0.0, 0.0);

        for (int child = 0; child < cfg.arity; ++child) {
          int rel_line = line;
          if (cfg.exchange_mode == ExchangeMode::AllGather) {
            rel_line = solve_redundant ? line : cfg.owned.first - cfg.prev_first + line;
          }
          const int offset = child * row_width * exchange_line_count + rel_line * row_width;
          const int row0 = 4 * child;
          for (int k = 0; k < 4; ++k) {
            const int row = row0 + k;
            a[band_index(row, bandwidth)] = Complex(1.0, 0.0);
            rhs[rhs_index(row, 0)] = exchanged(offset + k);
            if (child > 0) {
              int col = row0 - 2;
              a[band_index(row, bandwidth + col - row)] = -exchanged(offset + 4 + k);
              col = row0 - 1;
              a[band_index(row, bandwidth + col - row)] = -exchanged(offset + 8 + k);
            } else {
              rhs[rhs_index(row, 1)] = exchanged(offset + 4 + k);
              rhs[rhs_index(row, 2)] = exchanged(offset + 8 + k);
            }
            if (child < cfg.arity - 1) {
              int col = row0 + 4;
              a[band_index(row, bandwidth + col - row)] = -exchanged(offset + 12 + k);
              col = row0 + 5;
              a[band_index(row, bandwidth + col - row)] = -exchanged(offset + 16 + k);
            } else {
              rhs[rhs_index(row, 3)] = exchanged(offset + 12 + k);
              rhs[rhs_index(row, 4)] = exchanged(offset + 16 + k);
            }
          }
        }

        factor_banded_fixed(a, nrows);
        solve_factored_banded_multi_fixed(rhs, a, nrows, 5);

        for (int irhs = 0; irhs < 5; ++irhs) {
          for (int row = 0; row < nrows; ++row) {
            basis(recovery_index(line, irhs, row, nrows)) = rhs[rhs_index(row, irhs)];
          }
        }

        for (int k = 0; k < 4; ++k) {
          const int exposed_var = k < 2 ? k : 4 * (cfg.arity - 1) + k;
          rows(row_index(line, k)) = rhs[rhs_index(exposed_var, 0)];
          for (int ext_col = 0; ext_col < 4; ++ext_col) {
            rows(row_index(line, 4 * (ext_col + 1) + k)) = rhs[rhs_index(exposed_var, ext_col + 1)];
          }
        }
      });
  fence("schur_compose_level");
}

void SchurLevel::compose_from_exchanged_rows(bool solve_redundant) {
  compose_from_exchanged_rows(row_recv_, solve_redundant);
}

void SchurLevel::seed_root_values(bool solve_redundant) {
  const int nsolve = solve_count(solve_redundant);
  auto rows = rows_.view();
  auto values = values_.view();
  Kokkos::parallel_for(
      "schur_seed_root_values", Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>({0, 0}, {nsolve, value_width}),
      KOKKOS_LAMBDA(const int line, const int k) {
        values(value_index(line, k)) = k < 4 ? rows(row_index(line, k)) : Complex(0.0, 0.0);
      });
  fence("schur_seed_root_values");
}

void SchurLevel::pack_recovered_values(DeviceVector<Complex>& packed_values) const {
  const int nsolve = cfg_.owned.count;
  const int nrows = reduced_row_count();
  const auto expected = static_cast<std::size_t>(cfg_.arity) * static_cast<std::size_t>(value_width) *
                        static_cast<std::size_t>(nsolve);
  if (packed_values.size() < expected) {
    packed_values.resize("schur_packed_values", expected);
  }

  auto values = values_.view();
  auto basis = recovery_basis_.view();
  auto packed = packed_values.view();
  const auto cfg = cfg_;
  Kokkos::parallel_for(
      "schur_pack_recovered_values", Kokkos::RangePolicy<ExecutionSpace>(0, nsolve), KOKKOS_LAMBDA(const int line) {
        Complex recovered[max_rows];
        const Complex prev1 = values(value_index(line, 4));
        const Complex prev2 = values(value_index(line, 5));
        const Complex next1 = values(value_index(line, 6));
        const Complex next2 = values(value_index(line, 7));
        for (int k = 0; k < nrows; ++k) {
          recovered[k] = basis(recovery_index(line, 0, k, nrows)) +
                         basis(recovery_index(line, 1, k, nrows)) * prev1 +
                         basis(recovery_index(line, 2, k, nrows)) * prev2 +
                         basis(recovery_index(line, 3, k, nrows)) * next1 +
                         basis(recovery_index(line, 4, k, nrows)) * next2;
        }

        for (int child = 0; child < cfg.arity; ++child) {
          const int row0 = 4 * child;
          const int offset = child * value_width * nsolve + line * value_width;
          packed(offset + 0) = recovered[row0 + 0];
          packed(offset + 1) = recovered[row0 + 1];
          packed(offset + 2) = recovered[row0 + 2];
          packed(offset + 3) = recovered[row0 + 3];
          if (child > 0) {
            packed(offset + 4) = recovered[row0 - 2];
            packed(offset + 5) = recovered[row0 - 1];
          } else {
            packed(offset + 4) = prev1;
            packed(offset + 5) = prev2;
          }
          if (child < cfg.arity - 1) {
            packed(offset + 6) = recovered[row0 + 4];
            packed(offset + 7) = recovered[row0 + 5];
          } else {
            packed(offset + 6) = next1;
            packed(offset + 7) = next2;
          }
        }
      });
  fence("schur_pack_recovered_values");
}

void SchurLevel::pack_recovered_values() {
  pack_recovered_values(value_send_);
}

void SchurLevel::exchange_values() {
  if (cfg_.exchange_mode == ExchangeMode::AllGather) {
    throw std::runtime_error("SchurLevel::exchange_values is not used for allgather root recovery");
  }
  alltoallv_complex_device(value_send_.data(), value_send_counts_, value_send_displs_, value_recv_.data(),
                           value_recv_counts_, value_recv_displs_, cfg_.comm, "schur_values_alltoallv");
}

void SchurLevel::unpack_exchanged_values(DeviceVector<Complex>& child_values) const {
  const auto expected = static_cast<std::size_t>(value_width) * static_cast<std::size_t>(cfg_.prev_count);
  if (child_values.size() < expected) {
    child_values.resize("schur_child_values", expected);
  }
  auto recv = value_recv_.view();
  auto child = child_values.view();
  Kokkos::parallel_for(
      "schur_unpack_exchanged_values", Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(expected)),
      KOKKOS_LAMBDA(const int idx) { child(idx) = recv(idx); });
  fence("schur_unpack_exchanged_values");
}

void SchurLevel::recover_redundant_root_values(DeviceVector<Complex>& child_values) const {
  const int nsolve = cfg_.prev_count;
  const int nrows = reduced_row_count();
  const auto expected = static_cast<std::size_t>(value_width) * static_cast<std::size_t>(nsolve);
  if (child_values.size() < expected) {
    child_values.resize("schur_child_values", expected);
  }

  auto values = values_.view();
  auto basis = recovery_basis_.view();
  auto child_values_v = child_values.view();
  const auto cfg = cfg_;
  Kokkos::parallel_for(
      "schur_recover_redundant_root_values", Kokkos::RangePolicy<ExecutionSpace>(0, nsolve),
      KOKKOS_LAMBDA(const int line) {
        Complex recovered[max_rows];
        const Complex prev1 = values(value_index(line, 4));
        const Complex prev2 = values(value_index(line, 5));
        const Complex next1 = values(value_index(line, 6));
        const Complex next2 = values(value_index(line, 7));
        for (int k = 0; k < nrows; ++k) {
          recovered[k] = basis(recovery_index(line, 0, k, nrows)) +
                         basis(recovery_index(line, 1, k, nrows)) * prev1 +
                         basis(recovery_index(line, 2, k, nrows)) * prev2 +
                         basis(recovery_index(line, 3, k, nrows)) * next1 +
                         basis(recovery_index(line, 4, k, nrows)) * next2;
        }

        const int row0 = 4 * cfg.child_id;
        child_values_v(value_index(line, 0)) = recovered[row0 + 0];
        child_values_v(value_index(line, 1)) = recovered[row0 + 1];
        child_values_v(value_index(line, 2)) = recovered[row0 + 2];
        child_values_v(value_index(line, 3)) = recovered[row0 + 3];
        if (cfg.child_id > 0) {
          child_values_v(value_index(line, 4)) = recovered[row0 - 2];
          child_values_v(value_index(line, 5)) = recovered[row0 - 1];
        } else {
          child_values_v(value_index(line, 4)) = prev1;
          child_values_v(value_index(line, 5)) = prev2;
        }
        if (cfg.child_id < cfg.arity - 1) {
          child_values_v(value_index(line, 6)) = recovered[row0 + 4];
          child_values_v(value_index(line, 7)) = recovered[row0 + 5];
        } else {
          child_values_v(value_index(line, 6)) = next1;
          child_values_v(value_index(line, 7)) = next2;
        }
      });
  fence("schur_recover_redundant_root_values");
}

void SchurSolver::prepare(const SchurSolverConfig& cfg) {
  release();
  if (cfg.npy < 1) throw std::runtime_error("SchurSolver requires npy >= 1");
  if (cfg.ipy < 0 || cfg.ipy >= cfg.npy) throw std::runtime_error("SchurSolver received invalid y-rank topology");
  if (cfg.nlines < 1) throw std::runtime_error("SchurSolver requires nlines >= 1");
  if (!valid_schur_pass_sequence(cfg.npy, cfg.pass_counts)) {
    throw std::runtime_error("SchurSolver pass counts do not multiply to npy");
  }
  cfg_ = cfg;
  if (cfg.npy == 1) return;

  int span = 1;
  int prev_first = 0;
  int prev_count = cfg.nlines;
  const int pass_count = static_cast<int>(cfg.pass_counts.size());
  levels_.reserve(cfg.pass_counts.size());

  for (int ilevel = 0; ilevel < pass_count; ++ilevel) {
    const int arity = cfg.pass_counts[ilevel];
    const int child_span = span;
    span *= arity;
    const int rank_in_parent = cfg.ipy % span;
    const int child_id = rank_in_parent / child_span;
    const int child_pos = rank_in_parent % child_span;
    const int parent_group = cfg.ipy / span;
    const auto owned_rel = split_range(child_id, prev_count, arity);

    SchurLevelConfig level_cfg;
    level_cfg.level_index = ilevel;
    level_cfg.arity = arity;
    level_cfg.child_id = child_id;
    level_cfg.child_span = child_span;
    level_cfg.prev_first = prev_first;
    level_cfg.prev_count = prev_count;
    level_cfg.owned = {prev_first + owned_rel.first, owned_rel.count};
    level_cfg.exchange_mode = resolved_exchange_mode(cfg.exchange_mode, ilevel, pass_count, arity);
    level_cfg.comm = mpi_split(cfg.comm_y, parent_group * child_span + child_pos, child_id);

    levels_.emplace_back(std::move(level_cfg));
    prev_first = levels_.back().config().owned.first;
    prev_count = levels_.back().config().owned.count;
  }
}

void SchurSolver::release() {
  levels_.clear();
}

void SchurSolver::solve_from_leaf_rows(const DeviceVector<Complex>& leaf_rows, DeviceVector<Complex>& leaf_values) {
  const auto expected_rows = static_cast<std::size_t>(SchurLevel::row_width) * static_cast<std::size_t>(cfg_.nlines);
  if (leaf_rows.size() < expected_rows) {
    throw std::runtime_error("SchurSolver::solve_from_leaf_rows leaf rows too small");
  }
  const auto expected_values =
      static_cast<std::size_t>(SchurLevel::value_width) * static_cast<std::size_t>(cfg_.nlines);
  if (leaf_values.size() < expected_values) {
    leaf_values.resize("schur_leaf_values", expected_values);
  }
  if (levels_.empty()) {
    auto rows = leaf_rows.view();
    auto values = leaf_values.view();
    Kokkos::parallel_for(
        "schur_seed_single_rank_leaf_values",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>({0, 0}, {cfg_.nlines, SchurLevel::value_width}),
        KOKKOS_LAMBDA(const int line, const int k) {
          values(value_index(line, k)) = k < 4 ? rows(row_index(line, k)) : Complex(0.0, 0.0);
        });
    fence("schur_seed_single_rank_leaf_values");
    return;
  }

  const bool redundant_root = levels_.back().config().exchange_mode == ExchangeMode::AllGather;
  const DeviceVector<Complex>* current_rows = &leaf_rows;
  for (std::size_t ilevel = 0; ilevel < levels_.size(); ++ilevel) {
    auto& level = levels_[ilevel];
    level.pack_rows_for_exchange(*current_rows);
    level.exchange_rows();
    level.compose_from_exchanged_rows(redundant_root && ilevel + 1 == levels_.size());
    current_rows = &level.rows();
  }

  levels_.back().seed_root_values(redundant_root);
  for (std::size_t remaining = levels_.size(); remaining > 0; --remaining) {
    const std::size_t ilevel = remaining - 1;
    auto& level = levels_[ilevel];
    DeviceVector<Complex>& child_values = ilevel == 0 ? leaf_values : levels_[ilevel - 1].values();
    if (redundant_root && ilevel + 1 == levels_.size()) {
      level.recover_redundant_root_values(child_values);
      continue;
    }
    level.pack_recovered_values();
    level.exchange_values();
    level.unpack_exchanged_values(child_values);
  }
}

} // namespace channel
