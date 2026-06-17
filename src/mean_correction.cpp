#include "channel/mean_correction.hpp"

#include "channel/runtime.hpp"

#include <Kokkos_Core.hpp>

#include <algorithm>
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

  bcast_complex_device(line_.data(), static_cast<int>(line_.size()), root, comm, label + "_line");
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
  auto component = state.component_view(cfg_.component);
  auto local = local_line_.view();
  const int line = cfg_.line;
  const int line_count = static_cast<int>(grid_.line_count());
  const int local_n = grid_.ny;
  Kokkos::parallel_for(
      cfg_.label + "_pack", Kokkos::RangePolicy<ExecutionSpace>(0, local_n), KOKKOS_LAMBDA(const int y) {
        local(y) = component(static_cast<std::size_t>(y) * line_count + line);
      });
  fence((cfg_.label + "_pack").c_str());

  full_line_solver_.line().fill(Complex(0.0, 0.0));
  gatherv_complex_device(local_line_.data(), local_n, full_line_solver_.line().data(),
                         local_counts_, line_displs_, cfg_.root, cfg_.comm_y, cfg_.label + "_gather");
  full_line_solver_.solve_on_root_and_broadcast(cfg_.boundaries, cfg_.root, cfg_.comm_y, cfg_.label + "_solve");

  auto full_line = full_line_solver_.line().view();
  const int first_global_y = first_global_y_;
  Kokkos::parallel_for(
      cfg_.label + "_scatter", Kokkos::RangePolicy<ExecutionSpace>(0, local_n), KOKKOS_LAMBDA(const int y) {
        component(static_cast<std::size_t>(y) * line_count + line) = full_line(line_index(first_global_y + y));
      });
  fence((cfg_.label + "_scatter").c_str());
}

void DnsMeanCorrectionStep::check_grid(const DnsState& state, const char* caller) const {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

} // namespace channel
