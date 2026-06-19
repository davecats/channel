#include "channel/mean_correction.hpp"

#include "channel/runtime.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>

namespace channel {

namespace {

KOKKOS_INLINE_FUNCTION int coeff_index(int row, int offset) {
  return row * 5 + (offset + 2);
}

KOKKOS_INLINE_FUNCTION int line_index(int y) {
  return y + 1;
}

template <class Row>
KOKKOS_INLINE_FUNCTION double coeff(const Row& row, int offset) {
  return row[static_cast<std::size_t>(offset + 2)];
}

KOKKOS_INLINE_FUNCTION void factor_penta(double* mat, int n) {
  for (int i = 0; i < n; ++i) {
    const double piv = mat[coeff_index(i, 0)];
    mat[coeff_index(i, 0)] = 1.0 / piv;

    if (i + 1 < n) {
      const double factor = mat[coeff_index(i + 1, -1)] * mat[coeff_index(i, 0)];
      mat[coeff_index(i + 1, -1)] = factor;
      mat[coeff_index(i + 1, 0)] -= factor * mat[coeff_index(i, 1)];
      if (i + 2 < n) mat[coeff_index(i + 1, 1)] -= factor * mat[coeff_index(i, 2)];
    }

    if (i + 2 < n) {
      const double factor = mat[coeff_index(i + 2, -2)] * mat[coeff_index(i, 0)];
      mat[coeff_index(i + 2, -2)] = factor;
      mat[coeff_index(i + 2, -1)] -= factor * mat[coeff_index(i, 1)];
      mat[coeff_index(i + 2, 0)] -= factor * mat[coeff_index(i, 2)];
    }
  }
}

KOKKOS_INLINE_FUNCTION void solve_factored_penta(Complex* rhs, const double* mat, int n) {
  for (int i = 0; i < n; ++i) {
    if (i >= 2) rhs[i] -= mat[coeff_index(i, -2)] * rhs[i - 2];
    if (i >= 1) rhs[i] -= mat[coeff_index(i, -1)] * rhs[i - 1];
  }

  for (int i = n - 1; i >= 0; --i) {
    if (i + 1 < n) rhs[i] -= mat[coeff_index(i, 1)] * rhs[i + 1];
    if (i + 2 < n) rhs[i] -= mat[coeff_index(i, 2)] * rhs[i + 2];
    rhs[i] *= mat[coeff_index(i, 0)];
  }
}

std::vector<Complex> solve_full_line_dense_host(int ny,
                                                std::span<const double> matrix_rows,
                                                std::span<const Complex> rhs_line,
                                                const CompactBoundaryRows& boundaries) {
  const int n = ny + 3;
  std::vector<Complex> a(static_cast<std::size_t>(n * n), Complex(0.0, 0.0));
  std::vector<Complex> b(rhs_line.begin(), rhs_line.end());
  const auto dense = [&](int row, int col) -> Complex& {
    return a[static_cast<std::size_t>(row * n + col)];
  };
  const auto add_boundary = [&](int eq, int first_y, const std::array<double, 5>& row) {
    for (int k = 0; k < 5; ++k) {
      dense(eq, line_index(first_y + k)) = Complex(row[static_cast<std::size_t>(k)], 0.0);
    }
  };

  add_boundary(line_index(-1), -1, boundaries.lower_ghost);
  b[line_index(-1)] = boundaries.rhs_lower_ghost;
  add_boundary(line_index(0), -1, boundaries.lower);
  b[line_index(0)] = boundaries.rhs_lower;
  for (int iy = 1; iy <= ny - 1; ++iy) {
    for (int offset = -2; offset <= 2; ++offset) {
      dense(line_index(iy), line_index(iy + offset)) =
          Complex(matrix_rows[static_cast<std::size_t>(iy * 5 + offset + 2)], 0.0);
    }
  }
  add_boundary(line_index(ny), ny - 3, boundaries.upper);
  b[line_index(ny)] = boundaries.rhs_upper;
  add_boundary(line_index(ny + 1), ny - 3, boundaries.upper_ghost);
  b[line_index(ny + 1)] = boundaries.rhs_upper_ghost;

  for (int col = 0; col < n; ++col) {
    int pivot = col;
    const auto complex_abs = [](const Complex& value) {
      return std::hypot(static_cast<double>(value.real()), static_cast<double>(value.imag()));
    };
    double pivot_abs = complex_abs(dense(col, col));
    for (int row = col + 1; row < n; ++row) {
      const double candidate = complex_abs(dense(row, col));
      if (candidate > pivot_abs) {
        pivot = row;
        pivot_abs = candidate;
      }
    }
    if (pivot_abs == 0.0) {
      throw std::runtime_error("dense full-line correction solve found a singular pivot");
    }
    if (pivot != col) {
      for (int j = col; j < n; ++j) std::swap(dense(col, j), dense(pivot, j));
      std::swap(b[static_cast<std::size_t>(col)], b[static_cast<std::size_t>(pivot)]);
    }
    const Complex inv_pivot = Complex(1.0, 0.0) / dense(col, col);
    for (int j = col; j < n; ++j) dense(col, j) *= inv_pivot;
    b[static_cast<std::size_t>(col)] *= inv_pivot;
    for (int row = 0; row < n; ++row) {
      if (row == col) continue;
      const Complex factor = dense(row, col);
      if (factor == Complex(0.0, 0.0)) continue;
      for (int j = col; j < n; ++j) dense(row, j) -= factor * dense(col, j);
      b[static_cast<std::size_t>(row)] -= factor * b[static_cast<std::size_t>(col)];
    }
  }
  return b;
}

} // namespace

void FullLineCompactSolver::resize(int ny) {
  if (ny < 4) {
    throw std::runtime_error("FullLineCompactSolver requires ny >= 4");
  }
  ny_ = ny;
  matrix_rows_.resize("full_line_compact_matrix_rows", static_cast<std::size_t>(ny + 1) * 5);
  reduced_matrix_.resize("full_line_compact_reduced_matrix", static_cast<std::size_t>(ny - 1) * 5);
  line_.resize("full_line_compact_line", static_cast<std::size_t>(ny + 3));
  rhs_.resize("full_line_compact_rhs", static_cast<std::size_t>(ny - 1));
}

void FullLineCompactSolver::copy_system_from_host(std::span<const double> matrix_rows,
                                                  std::span<const Complex> rhs_line) {
  if (ny_ < 4) {
    throw std::runtime_error("FullLineCompactSolver is not initialized");
  }
  if (matrix_rows.size() != matrix_rows_.size()) {
    throw std::runtime_error("FullLineCompactSolver matrix row size mismatch");
  }
  if (rhs_line.size() != line_.size()) {
    throw std::runtime_error("FullLineCompactSolver line size mismatch");
  }
  matrix_rows_.copy_from_host(matrix_rows);
  line_.copy_from_host(rhs_line);
}

void FullLineCompactSolver::solve_on_root_and_broadcast(const CompactBoundaryRows& boundaries,
                                                        int root,
                                                        MpiComm comm,
                                                        const std::string& label) {
  if (ny_ < 4) {
    throw std::runtime_error("FullLineCompactSolver is not initialized");
  }
  const int rank = mpi_rank(comm);
  const int size = mpi_size(comm);
  if (root < 0 || root >= size) {
    throw std::runtime_error("FullLineCompactSolver received invalid root rank");
  }

  if (rank == root) {
    Kokkos::Array<double, 5> lower{};
    Kokkos::Array<double, 5> lower_ghost{};
    Kokkos::Array<double, 5> upper{};
    Kokkos::Array<double, 5> upper_ghost{};
    for (int i = 0; i < 5; ++i) {
      lower[static_cast<std::size_t>(i)] = boundaries.lower[static_cast<std::size_t>(i)];
      lower_ghost[static_cast<std::size_t>(i)] = boundaries.lower_ghost[static_cast<std::size_t>(i)];
      upper[static_cast<std::size_t>(i)] = boundaries.upper[static_cast<std::size_t>(i)];
      upper_ghost[static_cast<std::size_t>(i)] = boundaries.upper_ghost[static_cast<std::size_t>(i)];
    }

    const Complex rhs_lower = boundaries.rhs_lower;
    const Complex rhs_lower_ghost = boundaries.rhs_lower_ghost;
    const Complex rhs_upper = boundaries.rhs_upper;
    const Complex rhs_upper_ghost = boundaries.rhs_upper_ghost;
    const int ny = ny_;
    auto rows_v = matrix_rows_.view();
    auto mat_v = reduced_matrix_.view();
    auto line_v = line_.view();
    auto rhs_v = rhs_.view();

    Kokkos::parallel_for(
        label, Kokkos::RangePolicy<ExecutionSpace>(0, 1), KOKKOS_LAMBDA(const int) {
          Kokkos::Array<double, 5> lower_eq0{};
          Kokkos::Array<double, 5> upper_eqn{};
          const Complex lower_rhs0 = rhs_lower - rhs_lower_ghost * lower[0] / lower_ghost[0];
          const Complex upper_rhsn = rhs_upper - rhs_upper_ghost * upper[4] / upper_ghost[4];
          for (int i = 0; i < 5; ++i) {
            lower_eq0[static_cast<std::size_t>(i)] =
                lower[static_cast<std::size_t>(i)] -
                lower_ghost[static_cast<std::size_t>(i)] * lower[0] / lower_ghost[0];
            upper_eqn[static_cast<std::size_t>(i)] =
                upper[static_cast<std::size_t>(i)] -
                upper_ghost[static_cast<std::size_t>(i)] * upper[4] / upper_ghost[4];
          }
          lower_eq0[0] = 0.0;
          upper_eqn[4] = 0.0;

          const int active_n = ny - 1;
          for (int local_idx = 0; local_idx < active_n; ++local_idx) {
            const int iy = local_idx + 1;
            Complex value = line_v(line_index(iy));
            double row[5];
            for (int k = 0; k < 5; ++k) {
              row[k] = rows_v(static_cast<std::size_t>(iy * 5 + k));
            }

            if (local_idx == 0) {
              double fac = row[0] / lower_ghost[0];
              value -= rhs_lower_ghost * fac;
              for (int k = 0; k < 5; ++k) row[k] -= lower_ghost[static_cast<std::size_t>(k)] * fac;
              row[0] = 0.0;

              fac = row[1] / lower_eq0[1];
              value -= lower_rhs0 * fac;
              for (int k = 0; k < 5; ++k) row[k] -= lower_eq0[static_cast<std::size_t>(k)] * fac;
              row[1] = 0.0;
            } else if (local_idx == 1) {
              const double fac = row[0] / lower_eq0[1];
              value -= lower_rhs0 * fac;
              for (int k = 0; k < 4; ++k) row[k] -= lower_eq0[static_cast<std::size_t>(k + 1)] * fac;
              row[0] = 0.0;
            }

            if (local_idx == active_n - 2) {
              const double fac = row[4] / upper_eqn[3];
              value -= upper_rhsn * fac;
              for (int k = 1; k < 5; ++k) row[k] -= upper_eqn[static_cast<std::size_t>(k - 1)] * fac;
              row[4] = 0.0;
            } else if (local_idx == active_n - 1) {
              double fac = row[4] / upper_ghost[4];
              value -= rhs_upper_ghost * fac;
              for (int k = 0; k < 5; ++k) row[k] -= upper_ghost[static_cast<std::size_t>(k)] * fac;
              row[4] = 0.0;

              fac = row[3] / upper_eqn[3];
              value -= upper_rhsn * fac;
              for (int k = 0; k < 5; ++k) row[k] -= upper_eqn[static_cast<std::size_t>(k)] * fac;
              row[3] = 0.0;
            }

            rhs_v(local_idx) = value;
            for (int k = 0; k < 5; ++k) {
              mat_v(static_cast<std::size_t>(local_idx * 5 + k)) = row[k];
            }
          }

          factor_penta(mat_v.data(), active_n);
          solve_factored_penta(rhs_v.data(), mat_v.data(), active_n);

          for (int iy = 1; iy <= ny - 1; ++iy) {
            line_v(line_index(iy)) = rhs_v(iy - 1);
          }
          line_v(line_index(0)) =
              (lower_rhs0 - coeff(lower_eq0, 0) * line_v(line_index(1)) -
               coeff(lower_eq0, 1) * line_v(line_index(2)) -
               coeff(lower_eq0, 2) * line_v(line_index(3))) /
              coeff(lower_eq0, -1);
          line_v(line_index(-1)) =
              (rhs_lower_ghost - coeff(lower_ghost, -1) * line_v(line_index(0)) -
               coeff(lower_ghost, 0) * line_v(line_index(1)) -
               coeff(lower_ghost, 1) * line_v(line_index(2)) -
               coeff(lower_ghost, 2) * line_v(line_index(3))) /
              coeff(lower_ghost, -2);
          line_v(line_index(ny)) =
              (upper_rhsn - coeff(upper_eqn, -2) * line_v(line_index(ny - 3)) -
               coeff(upper_eqn, -1) * line_v(line_index(ny - 2)) -
               coeff(upper_eqn, 0) * line_v(line_index(ny - 1))) /
              coeff(upper_eqn, 1);
          line_v(line_index(ny + 1)) =
              (rhs_upper_ghost - coeff(upper_ghost, -2) * line_v(line_index(ny - 3)) -
               coeff(upper_ghost, -1) * line_v(line_index(ny - 2)) -
               coeff(upper_ghost, 0) * line_v(line_index(ny - 1)) -
               coeff(upper_ghost, 1) * line_v(line_index(ny))) /
              coeff(upper_ghost, 2);
        });
    fence(label.c_str());
  }

  bcast_complex_device(line_.view(), static_cast<int>(line_.size()), root, comm, label + "_line");
}

void DnsMeanCorrectionStep::prepare(const DnsState& state, DnsMeanCorrectionConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  [[maybe_unused]] const auto component_origin = state.index(cfg.component, 0, cfg.line);
  if (grid.ny > std::numeric_limits<int>::max() / 2) {
    throw std::runtime_error("DnsMeanCorrectionStep local y count exceeds int range");
  }
  comm_rank_ = mpi_rank(cfg.comm_y);
  comm_size_ = mpi_size(cfg.comm_y);
  if (cfg.root < 0 || cfg.root >= comm_size_) {
    throw std::runtime_error("DnsMeanCorrectionStep received invalid root rank");
  }

  local_counts_.assign(static_cast<std::size_t>(comm_size_), 0);
  const int local_n = grid.ny;
  const int err = MPI_Allgather(&local_n, 1, MPI_INT, local_counts_.data(), 1, MPI_INT, cfg.comm_y);
  if (err != MPI_SUCCESS) {
    throw std::runtime_error("MPI_Allgather failed in DnsMeanCorrectionStep::prepare");
  }

  int active_n = 0;
  first_global_y_ = 1;
  line_displs_.assign(static_cast<std::size_t>(comm_size_), 0);
  for (int rank = 0; rank < comm_size_; ++rank) {
    if (rank == comm_rank_) first_global_y_ = active_n + 1;
    line_displs_[static_cast<std::size_t>(rank)] = active_n + 2;
    active_n += local_counts_[static_cast<std::size_t>(rank)];
  }
  global_ny_ = active_n + 1;
  if (global_ny_ < 4) {
    throw std::runtime_error("DnsMeanCorrectionStep requires global_ny >= 4");
  }

  grid_ = grid;
  cfg_ = std::move(cfg);
  matrix_rows_.resize(cfg_.label + "_matrix_rows", static_cast<std::size_t>(global_ny_ + 1) * 5);
  local_line_.resize(cfg_.label + "_local_line", static_cast<std::size_t>(grid_.ny));
  full_line_solver_.resize(global_ny_);

  std::vector<double> zero_matrix(matrix_rows_.size(), 0.0);
  std::vector<Complex> zero_line(static_cast<std::size_t>(global_ny_ + 3), Complex(0.0, 0.0));
  matrix_rows_.copy_from_host(zero_matrix);
  full_line_solver_.copy_system_from_host(zero_matrix, zero_line);
}

void DnsMeanCorrectionStep::copy_matrix_from_host(std::span<const double> matrix_rows) {
  if (global_ny_ < 4) {
    throw std::runtime_error("DnsMeanCorrectionStep is not initialized");
  }
  if (matrix_rows.size() != matrix_rows_.size()) {
    throw std::runtime_error("DnsMeanCorrectionStep matrix size mismatch");
  }
  matrix_rows_.copy_from_host(matrix_rows);
  std::vector<Complex> zero_line(static_cast<std::size_t>(global_ny_ + 3), Complex(0.0, 0.0));
  full_line_solver_.copy_system_from_host(matrix_rows, zero_line);
}

void DnsMeanCorrectionStep::apply(DnsState& state) {
  check_grid(state, "DnsMeanCorrectionStep::apply");
  auto component = state.component_view_3d(cfg_.component);
  auto local = local_line_.view();
  const int line = cfg_.line;
  const int z = line / grid_.nx;
  const int x = line - z * grid_.nx;
  const int local_n = grid_.ny;
  {
    Kokkos::Profiling::ScopedRegion region("mean_correction pack_line");
    Kokkos::parallel_for(
        cfg_.label + "_pack", Kokkos::RangePolicy<ExecutionSpace>(0, local_n), KOKKOS_LAMBDA(const int y) {
          local(y) = component(y, z, x);
        });
    fence((cfg_.label + "_pack").c_str());
  }

  full_line_solver_.line().fill(Complex(0.0, 0.0));
  {
    Kokkos::Profiling::ScopedRegion region("mean_correction gather_line");
    gatherv_complex_device(local_line_.view(), local_n, full_line_solver_.line().view(),
                           local_counts_, line_displs_, cfg_.root, cfg_.comm_y, cfg_.label + "_gather");
  }
  {
    Kokkos::Profiling::ScopedRegion region("mean_correction solve_line");
    full_line_solver_.solve_on_root_and_broadcast(cfg_.boundaries, cfg_.root, cfg_.comm_y, cfg_.label + "_solve");
  }

  auto full_line = full_line_solver_.line().view();
  const int first_global_y = first_global_y_;
  {
    Kokkos::Profiling::ScopedRegion region("mean_correction scatter_line");
    Kokkos::parallel_for(
        cfg_.label + "_scatter", Kokkos::RangePolicy<ExecutionSpace>(0, local_n), KOKKOS_LAMBDA(const int y) {
          component(y, z, x) = full_line(line_index(first_global_y + y));
        });
    fence((cfg_.label + "_scatter").c_str());
  }
}

void DnsMeanCorrectionStep::check_grid(const DnsState& state, const char* caller) const {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DnsVelocityMeanCorrectionStep::prepare(const DnsState& state, DnsVelocityMeanCorrectionConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  [[maybe_unused]] const auto eta_origin = state.index(cfg.eta_component, 0, cfg.line);
  [[maybe_unused]] const auto w_origin = state.index(cfg.w_component, 0, cfg.line);
  comm_rank_ = mpi_rank(cfg.comm_y);
  comm_size_ = mpi_size(cfg.comm_y);
  if (cfg.root < 0 || cfg.root >= comm_size_) {
    throw std::runtime_error("DnsVelocityMeanCorrectionStep received invalid root rank");
  }
  if (grid.active_y_count < 1) {
    throw std::runtime_error("DnsVelocityMeanCorrectionStep requires active y metadata");
  }

  const int local_active_n = grid.active_y_count == 0 ? grid.ny : grid.active_y_count;
  int global_active_n = 0;
  MPI_Allreduce(&local_active_n, &global_active_n, 1, MPI_INT, MPI_SUM, cfg.comm_y);
  global_ny_ = global_active_n + 1;
  if (global_ny_ < 4) {
    throw std::runtime_error("DnsVelocityMeanCorrectionStep requires global_ny >= 4");
  }
  if (static_cast<int>(cfg.y.size()) != global_ny_ + 3) {
    throw std::runtime_error("DnsVelocityMeanCorrectionStep y-coordinate size mismatch");
  }

  grid_ = grid;
  cfg_ = std::move(cfg);
  matrix_rows_host_.assign(static_cast<std::size_t>(global_ny_ + 1) * 5, 0.0);
  correction_solver_.resize(global_ny_);
  local_line_.resize(cfg_.label + "_local_line", static_cast<std::size_t>(global_ny_ + 3));
  full_line_.resize(cfg_.label + "_full_line", static_cast<std::size_t>(global_ny_ + 3));
  corrected_eta_line_.resize(cfg_.label + "_corrected_eta_line", static_cast<std::size_t>(global_ny_ + 3));
  corrected_w_line_.resize(cfg_.label + "_corrected_w_line", static_cast<std::size_t>(global_ny_ + 3));
  std::vector<Complex> correction_rhs(static_cast<std::size_t>(global_ny_ + 3), Complex(0.0, 0.0));
  correction_solver_.copy_system_from_host(matrix_rows_host_, correction_rhs);
}

void DnsVelocityMeanCorrectionStep::copy_matrix_from_host(std::span<const double> matrix_rows) {
  if (global_ny_ < 4) {
    throw std::runtime_error("DnsVelocityMeanCorrectionStep is not initialized");
  }
  if (matrix_rows.size() != matrix_rows_host_.size()) {
    throw std::runtime_error("DnsVelocityMeanCorrectionStep matrix size mismatch");
  }
  std::copy(matrix_rows.begin(), matrix_rows.end(), matrix_rows_host_.begin());
}

double DnsVelocityMeanCorrectionStep::integrate_real_line(const std::vector<Complex>& line) const {
  double value = 0.0;
  const auto at = [](int iy) -> std::size_t { return static_cast<std::size_t>(iy + 1); };
  for (int iy = 1; iy <= global_ny_; iy += 2) {
    const double yp1 = cfg_.y[at(iy + 1)] - cfg_.y[at(iy)];
    const double ym1 = cfg_.y[at(iy - 1)] - cfg_.y[at(iy)];
    const double a1 = -ym1 / 3.0 + yp1 / 6.0 + yp1 * yp1 / (6.0 * ym1);
    const double a3 = yp1 / 3.0 - ym1 / 6.0 - ym1 * ym1 / (6.0 * yp1);
    const double a2 = yp1 - ym1 - a1 - a3;
    value += a1 * line[at(iy - 1)].real() + a2 * line[at(iy)].real() + a3 * line[at(iy + 1)].real();
  }
  return value;
}

void DnsVelocityMeanCorrectionStep::apply(DnsState& state) {
  check_grid(state, "DnsVelocityMeanCorrectionStep::apply");
  const int z = cfg_.line / grid_.nx;
  const int x = cfg_.line - z * grid_.nx;
  const int active_first = grid_.active_y_storage_first();
  const int active_count = grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count;
  auto eta = state.component_view_3d(cfg_.eta_component);
  auto w = state.component_view_3d(cfg_.w_component);
  auto local_line = local_line_.view();
  const int global_ny = global_ny_;
  const int grid_y_first = grid_.y_first;
  const int grid_ny = grid_.ny;

  {
    Kokkos::Profiling::ScopedRegion region("velocity_mean_correction pack_line");
    local_line_.fill(Complex(0.0, 0.0));
    Kokkos::parallel_for(
        cfg_.label + "_pack_line",
        Kokkos::RangePolicy<ExecutionSpace>(0, grid_ny),
        KOKKOS_LAMBDA(const int y_storage) {
          const int global_y = y_storage + grid_y_first;
          if (global_y < -1 || global_y > global_ny + 1) return;
          const bool physical_boundary = (global_y == -1 || global_y == 0 ||
                                          global_y == global_ny || global_y == global_ny + 1);
          const bool active_row = (y_storage >= active_first && y_storage < active_first + active_count);
          if (!physical_boundary && !active_row) return;
          local_line(line_index(global_y)) = eta(y_storage, z, x);
        });
    fence((cfg_.label + "_pack_line").c_str());
  }

  {
    Kokkos::Profiling::ScopedRegion region("velocity_mean_correction line_allreduce");
    allreduce_sum_complex_device(local_line_.view(),
                                 full_line_.view(),
                                 static_cast<int>(full_line_.size()),
                                 cfg_.comm_y,
                                 cfg_.label + "_line_allreduce");
  }

  std::vector<Complex> zero_mode_u(full_line_.size(), Complex(0.0, 0.0));
  std::vector<Complex> zero_mode_w(full_line_.size(), Complex(0.0, 0.0));
  {
    Kokkos::Profiling::ScopedRegion region("velocity_mean_correction host_solve");
    const auto full_line = full_line_.copy_to_host();
    const auto at = [](int iy) -> std::size_t { return static_cast<std::size_t>(iy + 1); };
    for (std::size_t i = 0; i < full_line.size(); ++i) {
      zero_mode_u[i] = Complex(full_line[i].real(), 0.0);
      zero_mode_w[i] = Complex(full_line[i].imag(), 0.0);
    }

    std::vector<Complex> correction_rhs(full_line.size(), Complex(0.0, 0.0));
    for (int iy = 1; iy <= global_ny_ - 1; ++iy) {
      correction_rhs[at(iy)] = Complex(1.0, 0.0);
    }
    const auto correction = solve_full_line_dense_host(global_ny_, matrix_rows_host_, correction_rhs, cfg_.boundaries);

    const double fr_u = integrate_real_line(zero_mode_u);
    const double fr_w = integrate_real_line(zero_mode_w);
    const double fr_c = integrate_real_line(correction);
    double corrpx = 0.0;
    double corrpz = 0.0;
    if (std::abs(cfg_.meanflowx) > 1.0e-7) corrpx = (cfg_.meanflowx - fr_u) / fr_c;
    if (std::abs(cfg_.meanflowz) > 1.0e-7) corrpz = (cfg_.meanflowz - fr_w) / fr_c;
    const char* debug_env = std::getenv("CHANNEL_DEBUG_TIMESTEP");
    if (comm_rank_ == cfg_.root && debug_env != nullptr && std::string(debug_env).find("ghost") != std::string::npos) {
      std::cout << "CPP_MEAN label=" << cfg_.label
                << " fr_u=" << fr_u
                << " fr_w=" << fr_w
                << " fr_c=" << fr_c
                << " corrpx=" << corrpx
                << " corrpz=" << corrpz
                << '\n';
    }

    for (std::size_t i = 0; i < full_line.size(); ++i) {
      zero_mode_u[i] = Complex(zero_mode_u[i].real() + corrpx * correction[i].real(), 0.0);
      zero_mode_w[i] = Complex(zero_mode_w[i].real() + corrpz * correction[i].real(), 0.0);
    }
  }
  {
    Kokkos::Profiling::ScopedRegion region("velocity_mean_correction copy_correction_to_device");
    corrected_eta_line_.copy_from_host(zero_mode_u);
    corrected_w_line_.copy_from_host(zero_mode_w);
  }

  auto corrected_eta = corrected_eta_line_.view();
  auto corrected_w = corrected_w_line_.view();
  {
    Kokkos::Profiling::ScopedRegion region("velocity_mean_correction scatter_line");
    Kokkos::parallel_for(
        cfg_.label + "_scatter_line",
        Kokkos::RangePolicy<ExecutionSpace>(0, grid_ny),
        KOKKOS_LAMBDA(const int y_storage) {
          const int global_y = y_storage + grid_y_first;
          if (global_y < -1 || global_y > global_ny + 1) return;
          eta(y_storage, z, x) = corrected_eta(line_index(global_y));
          w(y_storage, z, x) = corrected_w(line_index(global_y));
        });
    fence((cfg_.label + "_scatter_line").c_str());
  }
}

void DnsVelocityMeanCorrectionStep::check_grid(const DnsState& state, const char* caller) const {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components ||
      grid.y_first != grid_.y_first || grid.active_y_first != grid_.active_y_first ||
      grid.active_y_count != grid_.active_y_count) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

} // namespace channel
