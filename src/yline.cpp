#include "channel/yline.hpp"

#include "channel/runtime.hpp"

#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <mpi.h>

namespace channel {

namespace {

KOKKOS_INLINE_FUNCTION constexpr int compact_derivative_index_device(int y, int order, int offset) {
  return (y * 4 + order) * 5 + (offset + 2);
}

KOKKOS_INLINE_FUNCTION int penta_index(int row, int line, int batch_count) {
  return line + row * batch_count;
}

KOKKOS_INLINE_FUNCTION int endpoint_slot_for_row(int coupled_row,
                                                  int local_n,
                                                  bool has_left_interface,
                                                  bool has_right_interface,
                                                  int left_count) {
  (void)left_count;
  if (has_left_interface) {
    if (coupled_row == 0) return 0;
    if (coupled_row == 1) return 1;
  }
  if (has_right_interface) {
    if (coupled_row == local_n - 2) return 2;
    if (coupled_row == local_n - 1) return 3;
  }
  return -1;
}

KOKKOS_INLINE_FUNCTION int endpoint_iface_for_row(int coupled_row,
                                                   int local_n,
                                                   bool has_left_interface,
                                                   bool has_right_interface,
                                                   int left_count) {
  const int slot = endpoint_slot_for_row(coupled_row, local_n, has_left_interface, has_right_interface, left_count);
  if (slot < 0) return -1;
  if (slot < 2) return slot;
  return has_left_interface ? slot : slot - 2;
}

KOKKOS_INLINE_FUNCTION int local_row_for_endpoint_slot(int slot, int local_n) {
  return slot < 2 ? slot : local_n + slot - 4;
}

KOKKOS_INLINE_FUNCTION Complex coeff_for_offset(Complex ds, Complex dl, Complex d, Complex du, Complex dw, int offset) {
  switch (offset) {
    case -2: return ds;
    case -1: return dl;
    case 0: return d;
    case 1: return du;
    case 2: return dw;
    default: return Complex(0.0, 0.0);
  }
}

KOKKOS_INLINE_FUNCTION void set_coeff_for_offset(Complex coeff,
                                                  Complex& ds,
                                                  Complex& dl,
                                                  Complex& d,
                                                  Complex& du,
                                                  Complex& dw,
                                                  int offset) {
  switch (offset) {
    case -2: ds = coeff; break;
    case -1: dl = coeff; break;
    case 0: d = coeff; break;
    case 1: du = coeff; break;
    case 2: dw = coeff; break;
    default: break;
  }
}

KOKKOS_INLINE_FUNCTION void solve_small_dense(Complex* a, Complex* rhs, int n, int nrhs) {
  for (int k = 0; k < n; ++k) {
    const Complex pivot = a[k * 4 + k];
    for (int row = k + 1; row < n; ++row) {
      const Complex factor = a[row * 4 + k] / pivot;
      a[row * 4 + k] = factor;
      for (int col = k + 1; col < n; ++col) {
        a[row * 4 + col] -= factor * a[k * 4 + col];
      }
      for (int r = 0; r < nrhs; ++r) {
        rhs[row * 5 + r] -= factor * rhs[k * 5 + r];
      }
    }
  }

  for (int r = 0; r < nrhs; ++r) {
    for (int row = n - 1; row >= 0; --row) {
      for (int col = row + 1; col < n; ++col) {
        rhs[row * 5 + r] -= a[row * 4 + col] * rhs[col * 5 + r];
      }
      rhs[row * 5 + r] /= a[row * 4 + row];
    }
  }
}

bool debug_yline_enabled() {
  const char* value = std::getenv("CHANNEL_DEBUG_TIMESTEP");
  if (value == nullptr) return false;
  const std::string mode(value);
  return mode.find("yline") != std::string::npos || mode.find("worstline") != std::string::npos;
}

void print_yline_worst(const char* label,
                       int component,
                       int local_n,
                       int batch_count,
                       const DeviceVector<Complex>& values) {
  if (!debug_yline_enabled()) return;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const auto host = values.copy_to_host();
  int worst_row = -1;
  int worst_line = -1;
  int worst_bad = 0;
  double worst_mag = -1.0;
  for (int row = 0; row < local_n; ++row) {
    for (int line = 0; line < batch_count; ++line) {
      const auto& value = host[static_cast<std::size_t>(penta_index(row, line, batch_count))];
      const double re = static_cast<double>(value.real());
      const double im = static_cast<double>(value.imag());
      const bool bad = !std::isfinite(re) || !std::isfinite(im);
      const double mag = bad ? std::numeric_limits<double>::infinity() : std::hypot(re, im);
      if ((bad && worst_bad == 0) || (bad == (worst_bad != 0) && mag > worst_mag)) {
        worst_row = row;
        worst_line = line;
        worst_bad = bad ? 1 : 0;
        worst_mag = mag;
      }
    }
  }
  std::cout << std::scientific << std::setprecision(16)
            << "CPP_YLINE_WORST " << label
            << " rank=" << rank
            << " component=" << component
            << " row=" << worst_row
            << " line=" << worst_line
            << " bad=" << worst_bad
            << " mag=" << worst_mag
            << '\n';
}

void print_yline_diag(const char* label, int component, int local_n, int batch_count, const DeviceVector<Complex>& d) {
  if (!debug_yline_enabled()) return;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const auto host = d.copy_to_host();
  int min_row = -1;
  int min_line = -1;
  int nonfinite = 0;
  double min_abs = std::numeric_limits<double>::infinity();
  for (int row = 0; row < local_n; ++row) {
    for (int line = 0; line < batch_count; ++line) {
      const auto& value = host[static_cast<std::size_t>(penta_index(row, line, batch_count))];
      const double re = static_cast<double>(value.real());
      const double im = static_cast<double>(value.imag());
      if (!std::isfinite(re) || !std::isfinite(im)) {
        ++nonfinite;
        continue;
      }
      const double mag = std::hypot(re, im);
      if (mag < min_abs) {
        min_abs = mag;
        min_row = row;
        min_line = line;
      }
    }
  }
  std::cout << std::scientific << std::setprecision(16)
            << "CPP_YLINE_DIAG " << label
            << " rank=" << rank
            << " component=" << component
            << " row=" << min_row
            << " line=" << min_line
            << " nonfinite=" << nonfinite
            << " min_abs=" << min_abs
            << '\n';
}

void print_yline_coeffs(const char* label,
                        int component,
                        int local_n,
                        int batch_count,
                        const DeviceVector<Complex>& ds,
                        const DeviceVector<Complex>& dl,
                        const DeviceVector<Complex>& d,
                        const DeviceVector<Complex>& du,
                        const DeviceVector<Complex>& dw) {
  if (!debug_yline_enabled()) return;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const std::array<std::vector<Complex>, 5> bands{
      ds.copy_to_host(), dl.copy_to_host(), d.copy_to_host(), du.copy_to_host(), dw.copy_to_host()};
  int worst_band = -1;
  int worst_row = -1;
  int worst_line = -1;
  int worst_bad = 0;
  double worst_mag = -1.0;
  int nonfinite = 0;
  for (int band = 0; band < 5; ++band) {
    for (int row = 0; row < local_n; ++row) {
      for (int line = 0; line < batch_count; ++line) {
        const auto& value = bands[static_cast<std::size_t>(band)]
                                 [static_cast<std::size_t>(penta_index(row, line, batch_count))];
        const double re = static_cast<double>(value.real());
        const double im = static_cast<double>(value.imag());
        const bool bad = !std::isfinite(re) || !std::isfinite(im);
        if (bad) ++nonfinite;
        const double mag = bad ? std::numeric_limits<double>::infinity() : std::hypot(re, im);
        if ((bad && worst_bad == 0) || (bad == (worst_bad != 0) && mag > worst_mag)) {
          worst_band = band;
          worst_row = row;
          worst_line = line;
          worst_bad = bad ? 1 : 0;
          worst_mag = mag;
        }
      }
    }
  }
  std::cout << std::scientific << std::setprecision(16)
            << "CPP_YLINE_COEFFS " << label
            << " rank=" << rank
            << " component=" << component
            << " band=" << worst_band
            << " row=" << worst_row
            << " line=" << worst_line
            << " bad=" << worst_bad
            << " nonfinite=" << nonfinite
            << " mag=" << worst_mag
            << '\n';
}

void print_endpoint_vector(const char* label, int entry_width, int batch_count, const DeviceVector<Complex>& values) {
  if (!debug_yline_enabled()) return;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  const auto host = values.copy_to_host();
  int worst_line = -1;
  int worst_k = -1;
  int worst_bad = 0;
  int nonfinite = 0;
  double worst_mag = -1.0;
  for (int line = 0; line < batch_count; ++line) {
    for (int k = 0; k < entry_width; ++k) {
      const auto& value = host[static_cast<std::size_t>(line * entry_width + k)];
      const double re = static_cast<double>(value.real());
      const double im = static_cast<double>(value.imag());
      const bool bad = !std::isfinite(re) || !std::isfinite(im);
      if (bad) ++nonfinite;
      const double mag = bad ? std::numeric_limits<double>::infinity() : std::hypot(re, im);
      if ((bad && worst_bad == 0) || (bad == (worst_bad != 0) && mag > worst_mag)) {
        worst_line = line;
        worst_k = k;
        worst_bad = bad ? 1 : 0;
        worst_mag = mag;
      }
    }
  }
  std::cout << std::scientific << std::setprecision(16)
            << "CPP_ENDPOINT " << label
            << " rank=" << rank
            << " line=" << worst_line
            << " k=" << worst_k
            << " bad=" << worst_bad
            << " nonfinite=" << nonfinite
            << " mag=" << worst_mag;
  if (worst_line >= 0) {
    std::cout << " values=";
    for (int k = 0; k < entry_width; ++k) {
      const auto& value = host[static_cast<std::size_t>(worst_line * entry_width + k)];
      std::cout << " (" << value.real() << "," << value.imag() << ")";
    }
  }
  std::cout << '\n';
}

} // namespace

void YLineWorkspace::resize(int n, int batch_count) {
  if (n < 1 || batch_count < 1) {
    throw std::runtime_error("YLineWorkspace::resize requires n >= 1 and batch_count >= 1");
  }
  n_ = n;
  batch_count_ = batch_count;
  const auto count = static_cast<std::size_t>(n) * static_cast<std::size_t>(batch_count);
  ds_.resize("yline_ds", count);
  dl_.resize("yline_dl", count);
  d_.resize("yline_d", count);
  du_.resize("yline_du", count);
  dw_.resize("yline_dw", count);
  x_.resize("yline_x", count);
}

void YLineWorkspace::copy_from_host(const std::vector<Complex>& ds,
                                    const std::vector<Complex>& dl,
                                    const std::vector<Complex>& d,
                                    const std::vector<Complex>& du,
                                    const std::vector<Complex>& dw,
                                    const std::vector<Complex>& x) {
  const auto expected = static_cast<std::size_t>(n_) * static_cast<std::size_t>(batch_count_);
  if (ds.size() != expected || dl.size() != expected || d.size() != expected ||
      du.size() != expected || dw.size() != expected || x.size() != expected) {
    throw std::runtime_error("YLineWorkspace::copy_from_host size mismatch");
  }
  ds_.copy_from_host(ds);
  dl_.copy_from_host(dl);
  d_.copy_from_host(d);
  du_.copy_from_host(du);
  dw_.copy_from_host(dw);
  x_.copy_from_host(x);
}

void YLineWorkspace::solve(const std::string& label) {
  solve_batched_pentadiagonal(ds_, dl_, d_, du_, dw_, x_, n_, batch_count_, label);
}

void solve_batched_pentadiagonal(DeviceVector<Complex>& ds,
                                  DeviceVector<Complex>& dl,
                                  DeviceVector<Complex>& d,
                                  DeviceVector<Complex>& du,
                                  DeviceVector<Complex>& dw,
                                  DeviceVector<Complex>& x,
                                  int n,
                                  int batch_count,
                                  const std::string& label) {
  if (n < 1 || batch_count < 1) {
    throw std::runtime_error("solve_batched_pentadiagonal requires positive sizes");
  }
  const auto expected = static_cast<std::size_t>(n) * static_cast<std::size_t>(batch_count);
  if (ds.size() < expected || dl.size() < expected || d.size() < expected ||
      du.size() < expected || dw.size() < expected || x.size() < expected) {
    throw std::runtime_error("solve_batched_pentadiagonal input too small");
  }

  auto ds_v = ds.view();
  auto dl_v = dl.view();
  auto d_v = d.view();
  auto du_v = du.view();
  auto dw_v = dw.view();
  auto x_v = x.view();
  Kokkos::parallel_for(
      label, Kokkos::RangePolicy<ExecutionSpace>(0, batch_count), KOKKOS_LAMBDA(const int line) {
        if (n == 1) {
          const int p = line;
          const Complex inv = Complex(1.0, 0.0) / d_v(p);
          d_v(p) = inv;
          x_v(p) *= inv;
          return;
        }

        for (int i = 0; i <= n - 3; ++i) {
          const int p = line + i * batch_count;
          d_v(p) = Complex(1.0, 0.0) / d_v(p);

          const int p1 = p + batch_count;
          Complex factor = dl_v(p1) * d_v(p);
          dl_v(p1) = factor;
          d_v(p1) -= factor * du_v(p);
          du_v(p1) -= factor * dw_v(p);

          const int p2 = p + 2 * batch_count;
          factor = ds_v(p2) * d_v(p);
          ds_v(p2) = factor;
          dl_v(p2) -= factor * du_v(p);
          d_v(p2) -= factor * dw_v(p);
        }

        int p = line + (n - 2) * batch_count;
        d_v(p) = Complex(1.0, 0.0) / d_v(p);
        const int p1 = p + batch_count;
        const Complex factor = dl_v(p1) * d_v(p);
        dl_v(p1) = factor;
        d_v(p1) -= factor * du_v(p);
        d_v(p1) = Complex(1.0, 0.0) / d_v(p1);

        if (n >= 2) {
          p = line + batch_count;
          x_v(p) -= dl_v(p) * x_v(line);
        }
        for (int i = 2; i < n; ++i) {
          p = line + i * batch_count;
          x_v(p) -= ds_v(p) * x_v(line + (i - 2) * batch_count) +
                    dl_v(p) * x_v(line + (i - 1) * batch_count);
        }

        p = line + (n - 1) * batch_count;
        x_v(p) *= d_v(p);
        if (n >= 2) {
          p = line + (n - 2) * batch_count;
          x_v(p) = (x_v(p) - du_v(p) * x_v(p + batch_count)) * d_v(p);
        }
        for (int i = n - 3; i >= 0; --i) {
          p = line + i * batch_count;
          x_v(p) = (x_v(p) - du_v(p) * x_v(p + batch_count) -
                    dw_v(p) * x_v(p + 2 * batch_count)) *
                   d_v(p);
        }
      });
  fence(label.c_str());
}

int EndpointSchurYLineSolver::exposed_count() const {
  const bool has_left_interface = cfg_.ipy > 0;
  const bool has_right_interface = cfg_.ipy + 1 < cfg_.npy;
  return (has_left_interface ? 2 : 0) + (has_right_interface ? 2 : 0);
}

int EndpointSchurYLineSolver::interior_count() const {
  return cfg_.local_n - exposed_count();
}

void EndpointSchurYLineSolver::prepare(EndpointSchurYLineConfig cfg) {
  if (cfg.local_n < 1 || cfg.batch_count < 1) {
    throw std::runtime_error("EndpointSchurYLineSolver requires local_n >= 1 and batch_count >= 1");
  }
  if (cfg.npy < 1 || cfg.ipy < 0 || cfg.ipy >= cfg.npy) {
    throw std::runtime_error("EndpointSchurYLineSolver received invalid y-rank topology");
  }
  if (cfg.pass_counts.empty() && cfg.npy > 1) {
    cfg.pass_counts = default_schur_pass_counts(cfg.npy);
  }
  cfg_ = std::move(cfg);

  const int exposed = exposed_count();
  const int interior = interior_count();
  if (interior < 1) {
    throw std::runtime_error("EndpointSchurYLineSolver requires at least one interior row");
  }

  const int reduced_batch = cfg_.batch_count * (1 + exposed);
  const auto reduced_elems = static_cast<std::size_t>(interior) * static_cast<std::size_t>(reduced_batch);
  batch_ds_.resize("endpoint_schur_batch_ds", reduced_elems);
  batch_dl_.resize("endpoint_schur_batch_dl", reduced_elems);
  batch_d_.resize("endpoint_schur_batch_d", reduced_elems);
  batch_du_.resize("endpoint_schur_batch_du", reduced_elems);
  batch_dw_.resize("endpoint_schur_batch_dw", reduced_elems);
  batch_x_.resize("endpoint_schur_batch_x", reduced_elems);
  leaf_rows_.resize("endpoint_schur_leaf_rows",
                    static_cast<std::size_t>(SchurLevel::row_width) * static_cast<std::size_t>(cfg_.batch_count));
  leaf_values_.resize("endpoint_schur_leaf_values",
                      static_cast<std::size_t>(SchurLevel::value_width) * static_cast<std::size_t>(cfg_.batch_count));

  SchurSolverConfig schur_cfg;
  schur_cfg.npy = cfg_.npy;
  schur_cfg.ipy = cfg_.ipy;
  schur_cfg.nlines = cfg_.batch_count;
  schur_cfg.pass_counts = cfg_.pass_counts;
  schur_cfg.exchange_mode = cfg_.exchange_mode;
  schur_cfg.comm_y = cfg_.comm_y;
  schur_.prepare(schur_cfg);
}

void EndpointSchurYLineSolver::solve(DeviceVector<Complex>& ds,
                                      DeviceVector<Complex>& dl,
                                      DeviceVector<Complex>& d,
                                      DeviceVector<Complex>& du,
                                      DeviceVector<Complex>& dw,
                                      DeviceVector<Complex>& x) {
  const auto expected =
      static_cast<std::size_t>(cfg_.local_n) * static_cast<std::size_t>(cfg_.batch_count);
  if (ds.size() < expected || dl.size() < expected || d.size() < expected ||
      du.size() < expected || dw.size() < expected || x.size() < expected) {
    throw std::runtime_error("EndpointSchurYLineSolver::solve input too small");
  }

  const bool has_left_interface = cfg_.ipy > 0;
  const bool has_right_interface = cfg_.ipy + 1 < cfg_.npy;
  const int left_count = has_left_interface ? 2 : 0;
  const int exposed = exposed_count();
  const int interior_base = left_count;
  const int interior = interior_count();
  const int reduced_batch = cfg_.batch_count * (1 + exposed);

  if (exposed == 0) {
    Kokkos::Profiling::ScopedRegion region("endpoint_schur single_rank_pentadiagonal");
    solve_batched_pentadiagonal(ds, dl, d, du, dw, x, cfg_.local_n, cfg_.batch_count,
                                "endpoint_schur_single_rank_pentadiagonal");
    return;
  }

  auto ds_v = ds.view();
  auto dl_v = dl.view();
  auto d_v = d.view();
  auto du_v = du.view();
  auto dw_v = dw.view();
  auto x_v = x.view();
  auto bds = batch_ds_.view();
  auto bdl = batch_dl_.view();
  auto bd = batch_d_.view();
  auto bdu = batch_du_.view();
  auto bdw = batch_dw_.view();
  auto bx = batch_x_.view();
  const auto cfg = cfg_;

  {
    Kokkos::Profiling::ScopedRegion region("endpoint_schur pack_interior_systems");
    Kokkos::parallel_for(
        "endpoint_schur_pack_interior_systems",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>({0, 0}, {interior, reduced_batch}),
        KOKKOS_LAMBDA(const int local_i, const int sys) {
          int line = sys;
          int response_slot = -1;
          if (sys >= cfg.batch_count) {
            const int rel = sys - cfg.batch_count;
            response_slot = rel / cfg.batch_count;
            line = rel - response_slot * cfg.batch_count;
          }

          const int local_row = interior_base + local_i;
          const int actual_p = penta_index(local_row, line, cfg.batch_count);
          const Complex row_ds = ds_v(actual_p);
          const Complex row_dl = dl_v(actual_p);
          const Complex row_d = d_v(actual_p);
          const Complex row_du = du_v(actual_p);
          const Complex row_dw = dw_v(actual_p);

          Complex out_ds(0.0, 0.0);
          Complex out_dl(0.0, 0.0);
          Complex out_d(0.0, 0.0);
          Complex out_du(0.0, 0.0);
          Complex out_dw(0.0, 0.0);
          Complex rhs = response_slot < 0 ? x_v(actual_p) : Complex(0.0, 0.0);

          for (int offset = -2; offset <= 2; ++offset) {
            const Complex coeff = coeff_for_offset(row_ds, row_dl, row_d, row_du, row_dw, offset);
            if (coeff == Complex(0.0, 0.0)) continue;
            const int coupled_row = local_row + offset;
            if (coupled_row >= interior_base && coupled_row < interior_base + interior) {
              set_coeff_for_offset(coeff, out_ds, out_dl, out_d, out_du, out_dw, offset);
            } else if (response_slot >= 0) {
              const int exposed_iface = endpoint_iface_for_row(coupled_row, cfg.local_n, has_left_interface,
                                                               has_right_interface, left_count);
              if (exposed_iface == response_slot) rhs = coeff;
            }
          }

          const int p = penta_index(local_i, sys, reduced_batch);
          bds(p) = out_ds;
          bdl(p) = out_dl;
          bd(p) = out_d;
          bdu(p) = out_du;
          bdw(p) = out_dw;
          bx(p) = rhs;
        });
    fence("endpoint_schur_pack_interior_systems");
  }

  {
    Kokkos::Profiling::ScopedRegion region("endpoint_schur interior_responses");
    solve_batched_pentadiagonal(batch_ds_, batch_dl_, batch_d_, batch_du_, batch_dw_, batch_x_, interior, reduced_batch,
                                "endpoint_schur_interior_responses");
  }

  auto leaf_rows = leaf_rows_.view();
  {
    Kokkos::Profiling::ScopedRegion region("endpoint_schur build_leaf_rows");
    Kokkos::parallel_for(
        "endpoint_schur_build_leaf_rows", Kokkos::RangePolicy<ExecutionSpace>(0, cfg_.batch_count),
        KOKKOS_LAMBDA(const int line) {
        Complex s4[16];
        Complex rhs4[20];
        for (int i = 0; i < 16; ++i) s4[i] = Complex(0.0, 0.0);
        for (int i = 0; i < 20; ++i) rhs4[i] = Complex(0.0, 0.0);
        for (int k = 0; k < SchurLevel::row_width; ++k) {
          leaf_rows(line * SchurLevel::row_width + k) = Complex(0.0, 0.0);
        }

        for (int iface = 0; iface < exposed; ++iface) {
          const int slot = has_left_interface ? iface : iface + 2;
          const int local_row = local_row_for_endpoint_slot(slot, cfg.local_n);
          const int actual_p = penta_index(local_row, line, cfg.batch_count);
          const Complex row_ds = ds_v(actual_p);
          const Complex row_dl = dl_v(actual_p);
          const Complex row_d = d_v(actual_p);
          const Complex row_du = du_v(actual_p);
          const Complex row_dw = dw_v(actual_p);
          rhs4[iface * 5] = x_v(actual_p);

          for (int offset = -2; offset <= 2; ++offset) {
            const Complex coeff = coeff_for_offset(row_ds, row_dl, row_d, row_du, row_dw, offset);
            if (coeff == Complex(0.0, 0.0)) continue;
            const int coupled_row = local_row + offset;
            if (coupled_row >= interior_base && coupled_row < interior_base + interior) {
              const int j = coupled_row - interior_base;
              rhs4[iface * 5] -= coeff * bx(penta_index(j, line, reduced_batch));
              for (int exposed_slot = 0; exposed_slot < exposed; ++exposed_slot) {
                const int response_sys = cfg.batch_count * (1 + exposed_slot) + line;
                s4[iface * 4 + exposed_slot] -= coeff * bx(penta_index(j, response_sys, reduced_batch));
              }
            } else {
              const int exposed_iface = endpoint_iface_for_row(coupled_row, cfg.local_n, has_left_interface,
                                                               has_right_interface, left_count);
              if (exposed_iface >= 0) {
                s4[iface * 4 + exposed_iface] += coeff;
              } else {
                int ghost_col = -1;
                if (coupled_row == -2) ghost_col = 0;
                if (coupled_row == -1) ghost_col = 1;
                if (coupled_row == cfg.local_n) ghost_col = 2;
                if (coupled_row == cfg.local_n + 1) ghost_col = 3;
                if (ghost_col >= 0) rhs4[iface * 5 + ghost_col + 1] += coeff;
              }
            }
          }
        }

        solve_small_dense(s4, rhs4, exposed, 5);

        for (int iface = 0; iface < exposed; ++iface) {
          const int slot = has_left_interface ? iface : iface + 2;
          leaf_rows(line * SchurLevel::row_width + slot) = rhs4[iface * 5];
          for (int ghost_col = 0; ghost_col < 4; ++ghost_col) {
            leaf_rows(line * SchurLevel::row_width + 4 * (ghost_col + 1) + slot) =
                -rhs4[iface * 5 + ghost_col + 1];
          }
        }
        });
    fence("endpoint_schur_build_leaf_rows");
  }
  print_endpoint_vector("leaf_rows_before_schur", SchurLevel::row_width, cfg_.batch_count, leaf_rows_);

  {
    Kokkos::Profiling::ScopedRegion region("endpoint_schur solve_leaf_rows");
    schur_.solve_from_leaf_rows(leaf_rows_, leaf_values_);
  }
  print_endpoint_vector("leaf_values_after_schur", SchurLevel::value_width, cfg_.batch_count, leaf_values_);

  auto leaf_values = leaf_values_.view();
  {
    Kokkos::Profiling::ScopedRegion region("endpoint_schur reconstruct_solution");
    Kokkos::parallel_for(
        "endpoint_schur_reconstruct_solution",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>({0, 0}, {cfg_.local_n, cfg_.batch_count}),
        KOKKOS_LAMBDA(const int local_row, const int line) {
        Complex iface_values[4];
        for (int k = 0; k < 4; ++k) iface_values[k] = Complex(0.0, 0.0);
        iface_values[0] = leaf_values(line * SchurLevel::value_width + 0);
        iface_values[1] = leaf_values(line * SchurLevel::value_width + 1);
        iface_values[2] = leaf_values(line * SchurLevel::value_width + 2);
        iface_values[3] = leaf_values(line * SchurLevel::value_width + 3);

        const int exposed_slot = endpoint_slot_for_row(local_row, cfg.local_n, has_left_interface, has_right_interface,
                                                       left_count);
        const int p = penta_index(local_row, line, cfg.batch_count);
        if (exposed_slot >= 0) {
          x_v(p) = iface_values[exposed_slot];
          return;
        }

        const int local_i = local_row - interior_base;
        Complex value = bx(penta_index(local_i, line, reduced_batch));
        for (int iface = 0; iface < exposed; ++iface) {
          const int slot = has_left_interface ? iface : iface + 2;
          const int response_sys = cfg.batch_count * (1 + iface) + line;
          value -= bx(penta_index(local_i, response_sys, reduced_batch)) * iface_values[slot];
        }
        x_v(p) = value;
        });
    fence("endpoint_schur_reconstruct_solution");
  }
}

void DnsYLineComponentSolver::prepare(const DnsState& state, DnsYLineComponentConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  [[maybe_unused]] const auto component_origin = state.index(cfg.component, 0, 0);
  if (grid.line_count() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsYLineComponentSolver line count exceeds int range");
  }
  if (grid.active_values_per_component() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsYLineComponentSolver component size exceeds int range");
  }
  cfg_ = std::move(cfg);
  grid_ = grid;

  const auto count = grid_.active_values_per_component();
  ds_.resize("dns_yline_ds", count);
  dl_.resize("dns_yline_dl", count);
  d_.resize("dns_yline_d", count);
  du_.resize("dns_yline_du", count);
  dw_.resize("dns_yline_dw", count);
  work_ds_.resize("dns_yline_work_ds", count);
  work_dl_.resize("dns_yline_work_dl", count);
  work_d_.resize("dns_yline_work_d", count);
  work_du_.resize("dns_yline_work_du", count);
  work_dw_.resize("dns_yline_work_dw", count);
  rhs_.resize("dns_yline_rhs", count);
  const auto line_count = grid_.line_count();
  interface_send_.resize("dns_yline_interface_send", 2 * line_count);
  interface_recv_.resize("dns_yline_interface_recv", 2 * line_count);
  lower_ghost_rhs_.resize("dns_yline_lower_ghost_rhs", line_count);
  lower_boundary_rhs_.resize("dns_yline_lower_boundary_rhs", line_count);
  upper_boundary_rhs_.resize("dns_yline_upper_boundary_rhs", line_count);
  upper_ghost_rhs_.resize("dns_yline_upper_ghost_rhs", line_count);
  lower_ghost_rhs_.fill(Complex(0.0, 0.0));
  lower_boundary_rhs_.fill(Complex(0.0, 0.0));
  upper_boundary_rhs_.fill(Complex(0.0, 0.0));
  upper_ghost_rhs_.fill(Complex(0.0, 0.0));
  lower_boundary_enabled_ = false;
  upper_boundary_enabled_ = false;

  EndpointSchurYLineConfig endpoint_cfg;
  endpoint_cfg.local_n = grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count;
  endpoint_cfg.batch_count = static_cast<int>(grid_.line_count());
  endpoint_cfg.npy = cfg_.npy;
  endpoint_cfg.ipy = cfg_.ipy;
  endpoint_cfg.pass_counts = cfg_.pass_counts;
  endpoint_cfg.exchange_mode = cfg_.exchange_mode;
  endpoint_cfg.comm_y = cfg_.comm_y;
  endpoint_solver_.prepare(std::move(endpoint_cfg));
}

void DnsYLineComponentSolver::copy_coefficients_from_host(std::span<const Complex> ds,
                                                          std::span<const Complex> dl,
                                                          std::span<const Complex> d,
                                                          std::span<const Complex> du,
                                                          std::span<const Complex> dw) {
  const auto expected = grid_.active_values_per_component();
  if (ds.size() != expected || dl.size() != expected || d.size() != expected ||
      du.size() != expected || dw.size() != expected) {
    throw std::runtime_error("DnsYLineComponentSolver::copy_coefficients_from_host size mismatch");
  }
  ds_.copy_from_host(ds);
  dl_.copy_from_host(dl);
  d_.copy_from_host(d);
  du_.copy_from_host(du);
  dw_.copy_from_host(dw);
}

void DnsYLineComponentSolver::assemble_velocity_coefficients(const DeviceVector<double>& k2,
                                                             const DeviceVector<double>& derivatives,
                                                             int active_y_global_first,
                                                             int global_y_count,
                                                             double lambda,
                                                             double viscosity,
                                                             bool biharmonic) {
  const int lines = static_cast<int>(grid_.line_count());
  const int active_count = grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count;
  const auto expected = grid_.active_values_per_component();
  if (k2.size() != static_cast<std::size_t>(lines) ||
      derivatives.size() != static_cast<std::size_t>(active_count * 4 * 5) ||
      ds_.size() != expected) {
    throw std::runtime_error("DnsYLineComponentSolver::assemble_velocity_coefficients size mismatch");
  }

  auto ds = ds_.view();
  auto dl = dl_.view();
  auto d = d_.view();
  auto du = du_.view();
  auto dw = dw_.view();
  auto k2_v = k2.view();
  auto der = derivatives.view();
  Kokkos::parallel_for(
      "dns_yline_assemble_velocity_coefficients",
      Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(expected)),
      KOKKOS_LAMBDA(const int p) {
        const int y = p / lines;
        const int line = p - y * lines;
        const int global_y = active_y_global_first + y;
        if (global_y == 0 || global_y == global_y_count - 1) {
          ds(p) = Complex(0.0, 0.0);
          dl(p) = Complex(0.0, 0.0);
          d(p) = Complex(1.0, 0.0);
          du(p) = Complex(0.0, 0.0);
          dw(p) = Complex(0.0, 0.0);
          return;
        }

        const double k2_line = k2_v(line);
        Complex row[5];
        for (int offset = -2; offset <= 2; ++offset) {
          const double d0 = der(compact_derivative_index_device(y, 0, offset));
          const double d2 = der(compact_derivative_index_device(y, 2, offset));
          const double d4 = der(compact_derivative_index_device(y, 3, offset));
          double value = 0.0;
          if (biharmonic) {
            value = lambda * (d2 - k2_line * d0) -
                    viscosity * (d4 - 2.0 * k2_line * d2 + k2_line * k2_line * d0);
          } else {
            value = lambda * d0 - viscosity * (d2 - k2_line * d0);
          }
          row[static_cast<std::size_t>(offset + 2)] = Complex(value, 0.0);
        }
        ds(p) = row[0];
        dl(p) = row[1];
        d(p) = row[2];
        du(p) = row[3];
        dw(p) = row[4];
      });
  fence("dns_yline_assemble_velocity_coefficients");
}

void DnsYLineComponentSolver::copy_boundary_data_from_host(const DnsYLineBoundaryData& boundary) {
  const auto expected = grid_.line_count();
  const auto check_rhs = [&](std::span<const Complex> rhs, const char* name) {
    if (rhs.size() != expected) {
      throw std::runtime_error(std::string("DnsYLineComponentSolver::copy_boundary_data_from_host ") + name +
                               " size mismatch");
    }
  };

  lower_boundary_enabled_ = boundary.lower_enabled;
  upper_boundary_enabled_ = boundary.upper_enabled;
  if (lower_boundary_enabled_) {
    check_rhs(boundary.lower_ghost_rhs, "lower_ghost_rhs");
    check_rhs(boundary.lower_boundary_rhs, "lower_boundary_rhs");
    for (int i = 0; i < 5; ++i) {
      lower_ghost_eq_[static_cast<std::size_t>(i)] = boundary.lower_ghost_eq[static_cast<std::size_t>(i)];
      lower_boundary_eq_[static_cast<std::size_t>(i)] = boundary.lower_boundary_eq[static_cast<std::size_t>(i)];
    }
    lower_ghost_rhs_.copy_from_host(boundary.lower_ghost_rhs);
    lower_boundary_rhs_.copy_from_host(boundary.lower_boundary_rhs);
  }
  if (upper_boundary_enabled_) {
    check_rhs(boundary.upper_boundary_rhs, "upper_boundary_rhs");
    check_rhs(boundary.upper_ghost_rhs, "upper_ghost_rhs");
    for (int i = 0; i < 5; ++i) {
      upper_boundary_eq_[static_cast<std::size_t>(i)] = boundary.upper_boundary_eq[static_cast<std::size_t>(i)];
      upper_ghost_eq_[static_cast<std::size_t>(i)] = boundary.upper_ghost_eq[static_cast<std::size_t>(i)];
    }
    upper_boundary_rhs_.copy_from_host(boundary.upper_boundary_rhs);
    upper_ghost_rhs_.copy_from_host(boundary.upper_ghost_rhs);
  }
}

void DnsYLineComponentSolver::copy_boundary_data_from_device(const DnsYLineDeviceBoundaryData& boundary) {
  const auto expected = grid_.line_count();
  const auto check_rhs = [&](const DeviceVector<Complex>* rhs, const char* name) {
    if (rhs == nullptr || rhs->size() != expected) {
      throw std::runtime_error(std::string("DnsYLineComponentSolver::copy_boundary_data_from_device ") + name +
                               " size mismatch");
    }
  };

  lower_boundary_enabled_ = boundary.lower_enabled;
  upper_boundary_enabled_ = boundary.upper_enabled;
  if (lower_boundary_enabled_) {
    check_rhs(boundary.lower_ghost_rhs, "lower_ghost_rhs");
    check_rhs(boundary.lower_boundary_rhs, "lower_boundary_rhs");
    for (int i = 0; i < 5; ++i) {
      lower_ghost_eq_[static_cast<std::size_t>(i)] = boundary.lower_ghost_eq[static_cast<std::size_t>(i)];
      lower_boundary_eq_[static_cast<std::size_t>(i)] = boundary.lower_boundary_eq[static_cast<std::size_t>(i)];
    }
    lower_ghost_rhs_.copy_from(*boundary.lower_ghost_rhs);
    lower_boundary_rhs_.copy_from(*boundary.lower_boundary_rhs);
  }
  if (upper_boundary_enabled_) {
    check_rhs(boundary.upper_boundary_rhs, "upper_boundary_rhs");
    check_rhs(boundary.upper_ghost_rhs, "upper_ghost_rhs");
    for (int i = 0; i < 5; ++i) {
      upper_boundary_eq_[static_cast<std::size_t>(i)] = boundary.upper_boundary_eq[static_cast<std::size_t>(i)];
      upper_ghost_eq_[static_cast<std::size_t>(i)] = boundary.upper_ghost_eq[static_cast<std::size_t>(i)];
    }
    upper_boundary_rhs_.copy_from(*boundary.upper_boundary_rhs);
    upper_ghost_rhs_.copy_from(*boundary.upper_ghost_rhs);
  }
}

void DnsYLineComponentSolver::solve(DnsState& state) {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error("DnsYLineComponentSolver::solve grid mismatch");
  }
  const char* f90_label = "linsolve solve_scalar";
  if (cfg_.component == 0) {
    f90_label = "linsolve solve_eta";
  } else if (cfg_.component == 1) {
    f90_label = "linsolve solve_v";
  }
  Kokkos::Profiling::ScopedRegion f90_region(f90_label);
  auto component = state.component_view_3d(cfg_.component);
  auto rhs = rhs_.view();
  const int line_count = static_cast<int>(grid_.line_count());
  const int nx = grid_.nx;
  const int active_first = grid_.active_y_storage_first();
  {
    Kokkos::Profiling::ScopedRegion region("dns_yline pack_component_rhs");
    Kokkos::parallel_for(
        "dns_yline_pack_component_rhs",
        Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(grid_.active_values_per_component())),
        KOKKOS_LAMBDA(const int i) {
          const int y_active = i / line_count;
          const int line = i - y_active * line_count;
          const int z = line / nx;
          const int x = line - z * nx;
          const int y_storage = active_first + y_active;
          rhs(i) = component(y_storage, z, x);
        });
    fence("dns_yline_pack_component_rhs");
  }
  const int active_n = grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count;
  print_yline_worst("packed_rhs", cfg_.component, active_n, line_count, rhs_);

  {
    Kokkos::Profiling::ScopedRegion region("dns_yline copy_coefficients_to_work");
    work_ds_.copy_from(ds_);
    work_dl_.copy_from(dl_);
    work_d_.copy_from(d_);
    work_du_.copy_from(du_);
    work_dw_.copy_from(dw_);
  }
  {
    Kokkos::Profiling::ScopedRegion region("dns_yline boundary_elimination");
    eliminate_boundaries();
  }
  print_yline_diag("after_boundary_elimination", cfg_.component, active_n, line_count, work_d_);
  print_yline_coeffs("after_boundary_elimination",
                     cfg_.component,
                     active_n,
                     line_count,
                     work_ds_,
                     work_dl_,
                     work_d_,
                     work_du_,
                     work_dw_);
  print_yline_worst("after_boundary_elimination_rhs", cfg_.component, active_n, line_count, rhs_);
  {
    Kokkos::Profiling::ScopedRegion region("dns_yline endpoint_schur");
    endpoint_solver_.solve(work_ds_, work_dl_, work_d_, work_du_, work_dw_, rhs_);
  }
  print_yline_worst("after_endpoint_solve", cfg_.component, active_n, line_count, rhs_);

  {
    Kokkos::Profiling::ScopedRegion region("dns_yline unpack_component_solution");
    Kokkos::parallel_for(
        "dns_yline_unpack_component_solution",
        Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(grid_.active_values_per_component())),
        KOKKOS_LAMBDA(const int i) {
          const int y_active = i / line_count;
          const int line = i - y_active * line_count;
          const int z = line / nx;
          const int x = line - z * nx;
          const int y_storage = active_first + y_active;
          component(y_storage, z, x) = rhs(i);
        });
    fence("dns_yline_unpack_component_solution");
  }
  if (cfg_.fill_interface_ghosts && cfg_.npy > 1) {
    Kokkos::Profiling::ScopedRegion region("dns_yline exchange_interface_ghosts");
    exchange_interface_ghosts(state);
  }
  {
    Kokkos::Profiling::ScopedRegion region("dns_yline reconstruct_boundaries");
    reconstruct_boundaries(state);
  }
}

void DnsYLineComponentSolver::exchange_interface_ghosts(DnsState& state) {
  const int active_first = grid_.active_y_storage_first();
  const int active_n = grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count;
  const int line_count = static_cast<int>(grid_.line_count());
  if (active_n < 2) {
    throw std::runtime_error("DnsYLineComponentSolver interface ghost exchange requires at least two active rows");
  }
  if (active_first < 2 && cfg_.ipy > 0) {
    throw std::runtime_error("DnsYLineComponentSolver missing lower interface ghost storage");
  }
  if (active_first + active_n + 1 >= grid_.ny && cfg_.ipy + 1 < cfg_.npy) {
    throw std::runtime_error("DnsYLineComponentSolver missing upper interface ghost storage");
  }

  auto component = state.component_view_3d(cfg_.component);
  auto send = interface_send_.view();
  auto recv = interface_recv_.view();
  constexpr int lower_tag = 7101;
  constexpr int upper_tag = 7102;
  const int count = 2 * line_count;
  const int nx = grid_.nx;

  if (cfg_.ipy > 0) {
    const int peer = cfg_.ipy - 1;
    Kokkos::parallel_for(
        "dns_yline_pack_lower_interface",
        Kokkos::RangePolicy<ExecutionSpace>(0, count),
        KOKKOS_LAMBDA(const int i) {
          const int y_offset = i / line_count;
          const int line = i - y_offset * line_count;
          const int z = line / nx;
          const int x = line - z * nx;
          send(i) = component(active_first + y_offset, z, x);
        });
    fence("dns_yline_pack_lower_interface");
    sendrecv_complex_device(interface_send_.view(),
                            count,
                            peer,
                            lower_tag,
                            interface_recv_.view(),
                            count,
                            peer,
                            upper_tag,
                            cfg_.comm_y,
                            "dns_yline_lower_interface");
    Kokkos::parallel_for(
        "dns_yline_unpack_lower_interface",
        Kokkos::RangePolicy<ExecutionSpace>(0, count),
        KOKKOS_LAMBDA(const int i) {
          const int y_offset = i / line_count;
          const int line = i - y_offset * line_count;
          const int z = line / nx;
          const int x = line - z * nx;
          component(active_first - 2 + y_offset, z, x) = recv(i);
        });
    fence("dns_yline_unpack_lower_interface");
  }
  if (cfg_.ipy + 1 < cfg_.npy) {
    const int peer = cfg_.ipy + 1;
    Kokkos::parallel_for(
        "dns_yline_pack_upper_interface",
        Kokkos::RangePolicy<ExecutionSpace>(0, count),
        KOKKOS_LAMBDA(const int i) {
          const int y_offset = i / line_count;
          const int line = i - y_offset * line_count;
          const int z = line / nx;
          const int x = line - z * nx;
          send(i) = component(active_first + active_n - 2 + y_offset, z, x);
        });
    fence("dns_yline_pack_upper_interface");
    sendrecv_complex_device(interface_send_.view(),
                            count,
                            peer,
                            upper_tag,
                            interface_recv_.view(),
                            count,
                            peer,
                            lower_tag,
                            cfg_.comm_y,
                            "dns_yline_upper_interface");
    Kokkos::parallel_for(
        "dns_yline_unpack_upper_interface",
        Kokkos::RangePolicy<ExecutionSpace>(0, count),
        KOKKOS_LAMBDA(const int i) {
          const int y_offset = i / line_count;
          const int line = i - y_offset * line_count;
          const int z = line / nx;
          const int x = line - z * nx;
          component(active_first + active_n + y_offset, z, x) = recv(i);
        });
    fence("dns_yline_unpack_upper_interface");
  }
}

void DnsYLineComponentSolver::eliminate_boundaries() {
  const int n = grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count;
  const int batch = static_cast<int>(grid_.line_count());
  if (lower_boundary_enabled_ && n < 2) {
    throw std::runtime_error("DnsYLineComponentSolver lower boundary elimination requires at least two active rows");
  }
  if (upper_boundary_enabled_ && n < 2) {
    throw std::runtime_error("DnsYLineComponentSolver upper boundary elimination requires at least two active rows");
  }

  auto ds = work_ds_.view();
  auto dl = work_dl_.view();
  auto d = work_d_.view();
  auto du = work_du_.view();
  auto dw = work_dw_.view();
  auto rhs = rhs_.view();

  if (lower_boundary_enabled_) {
    const auto eqm1 = lower_ghost_eq_;
    const auto eq0_raw = lower_boundary_eq_;
    const auto lower_ghost = lower_ghost_rhs_.view();
    const auto lower_boundary = lower_boundary_rhs_.view();
    Kokkos::parallel_for(
        "dns_yline_eliminate_lower_boundary",
        Kokkos::RangePolicy<ExecutionSpace>(0, batch),
        KOKKOS_LAMBDA(const int line) {
          const Complex boundary_rhs0 = lower_boundary(line) -
                                        lower_ghost(line) * eq0_raw[0] / eqm1[0];
          Kokkos::Array<Complex, 4> eq0{};
          for (int i = 0; i < 4; ++i) {
            eq0[static_cast<std::size_t>(i)] =
                Complex(eq0_raw[static_cast<std::size_t>(i + 1)], 0.0) -
                Complex(eqm1[static_cast<std::size_t>(i + 1)], 0.0) * eq0_raw[0] / eqm1[0];
          }

          int p = line;
          Complex fac = ds(p) / eqm1[0];
          rhs(p) -= lower_ghost(line) * fac;
          ds(p) -= Complex(eqm1[0], 0.0) * fac;
          dl(p) -= Complex(eqm1[1], 0.0) * fac;
          d(p) -= Complex(eqm1[2], 0.0) * fac;
          du(p) -= Complex(eqm1[3], 0.0) * fac;
          dw(p) -= Complex(eqm1[4], 0.0) * fac;
          ds(p) = Complex(0.0, 0.0);

          fac = dl(p) / eq0[0];
          rhs(p) -= boundary_rhs0 * fac;
          dl(p) -= eq0[0] * fac;
          d(p) -= eq0[1] * fac;
          du(p) -= eq0[2] * fac;
          dw(p) -= eq0[3] * fac;
          dl(p) = Complex(0.0, 0.0);

          p = batch + line;
          fac = ds(p) / eq0[0];
          rhs(p) -= boundary_rhs0 * fac;
          ds(p) -= eq0[0] * fac;
          dl(p) -= eq0[1] * fac;
          d(p) -= eq0[2] * fac;
          du(p) -= eq0[3] * fac;
          ds(p) = Complex(0.0, 0.0);
        });
    fence("dns_yline_eliminate_lower_boundary");
  }

  if (upper_boundary_enabled_) {
    const auto eqn_raw = upper_boundary_eq_;
    const auto eqnp1 = upper_ghost_eq_;
    const auto upper_boundary = upper_boundary_rhs_.view();
    const auto upper_ghost = upper_ghost_rhs_.view();
    Kokkos::parallel_for(
        "dns_yline_eliminate_upper_boundary",
        Kokkos::RangePolicy<ExecutionSpace>(0, batch),
        KOKKOS_LAMBDA(const int line) {
          const Complex boundary_rhsn = upper_boundary(line) -
                                        upper_ghost(line) * eqn_raw[4] / eqnp1[4];
          Kokkos::Array<Complex, 4> eqn{};
          for (int i = 0; i < 4; ++i) {
            eqn[static_cast<std::size_t>(i)] =
                Complex(eqn_raw[static_cast<std::size_t>(i)], 0.0) -
                Complex(eqnp1[static_cast<std::size_t>(i)], 0.0) * eqn_raw[4] / eqnp1[4];
          }

          int p = (n - 2) * batch + line;
          Complex fac = dw(p) / eqn[3];
          rhs(p) -= boundary_rhsn * fac;
          dl(p) -= eqn[0] * fac;
          d(p) -= eqn[1] * fac;
          du(p) -= eqn[2] * fac;
          dw(p) -= eqn[3] * fac;
          dw(p) = Complex(0.0, 0.0);

          p = (n - 1) * batch + line;
          fac = dw(p) / eqnp1[4];
          rhs(p) -= upper_ghost(line) * fac;
          ds(p) -= Complex(eqnp1[0], 0.0) * fac;
          dl(p) -= Complex(eqnp1[1], 0.0) * fac;
          d(p) -= Complex(eqnp1[2], 0.0) * fac;
          du(p) -= Complex(eqnp1[3], 0.0) * fac;
          dw(p) -= Complex(eqnp1[4], 0.0) * fac;
          dw(p) = Complex(0.0, 0.0);

          fac = du(p) / eqn[3];
          rhs(p) -= boundary_rhsn * fac;
          ds(p) -= eqn[0] * fac;
          dl(p) -= eqn[1] * fac;
          d(p) -= eqn[2] * fac;
          du(p) -= eqn[3] * fac;
          du(p) = Complex(0.0, 0.0);
        });
    fence("dns_yline_eliminate_upper_boundary");
  }
}

void DnsYLineComponentSolver::reconstruct_boundaries(DnsState& state) {
  if (!lower_boundary_enabled_ && !upper_boundary_enabled_) return;
  const int n = grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count;
  const int batch = static_cast<int>(grid_.line_count());
  const int active_first = grid_.active_y_storage_first();
  const int nx = grid_.nx;
  auto component = state.component_view_3d(cfg_.component);

  if (lower_boundary_enabled_) {
    const auto eqm1 = lower_ghost_eq_;
    const auto eq0_raw = lower_boundary_eq_;
    const auto lower_ghost = lower_ghost_rhs_.view();
    const auto lower_boundary = lower_boundary_rhs_.view();
    Kokkos::parallel_for(
        "dns_yline_reconstruct_lower_boundary",
        Kokkos::RangePolicy<ExecutionSpace>(0, batch),
        KOKKOS_LAMBDA(const int line) {
          const int z = line / nx;
          const int x = line - z * nx;
          const Complex boundary_rhs0 = lower_boundary(line) -
                                        lower_ghost(line) * eq0_raw[0] / eqm1[0];
          Kokkos::Array<Complex, 4> eq0{};
          for (int i = 0; i < 4; ++i) {
            eq0[static_cast<std::size_t>(i)] =
                Complex(eq0_raw[static_cast<std::size_t>(i + 1)], 0.0) -
                Complex(eqm1[static_cast<std::size_t>(i + 1)], 0.0) * eq0_raw[0] / eqm1[0];
          }
          const auto y1 = component(active_first + 0, z, x);
          const auto y2 = component(active_first + 1, z, x);
          const auto y3 = component(active_first + 2, z, x);
          const Complex y0 = (boundary_rhs0 - eq0[1] * y1 - eq0[2] * y2 - eq0[3] * y3) / eq0[0];
          const Complex ym1 = (lower_ghost(line) - Complex(eqm1[1], 0.0) * y0 -
                               Complex(eqm1[2], 0.0) * y1 - Complex(eqm1[3], 0.0) * y2 -
                               Complex(eqm1[4], 0.0) * y3) /
                              eqm1[0];
          component(active_first - 1, z, x) = y0;
          component(active_first - 2, z, x) = ym1;
        });
    fence("dns_yline_reconstruct_lower_boundary");
  }

  if (upper_boundary_enabled_) {
    const auto eqn_raw = upper_boundary_eq_;
    const auto eqnp1 = upper_ghost_eq_;
    const auto upper_boundary = upper_boundary_rhs_.view();
    const auto upper_ghost = upper_ghost_rhs_.view();
    Kokkos::parallel_for(
        "dns_yline_reconstruct_upper_boundary",
        Kokkos::RangePolicy<ExecutionSpace>(0, batch),
        KOKKOS_LAMBDA(const int line) {
          const int z = line / nx;
          const int x = line - z * nx;
          const Complex boundary_rhsn = upper_boundary(line) -
                                        upper_ghost(line) * eqn_raw[4] / eqnp1[4];
          Kokkos::Array<Complex, 4> eqn{};
          for (int i = 0; i < 4; ++i) {
            eqn[static_cast<std::size_t>(i)] =
                Complex(eqn_raw[static_cast<std::size_t>(i)], 0.0) -
                Complex(eqnp1[static_cast<std::size_t>(i)], 0.0) * eqn_raw[4] / eqnp1[4];
          }
          const auto ynm3 = component(active_first + n - 3, z, x);
          const auto ynm2 = component(active_first + n - 2, z, x);
          const auto ynm1 = component(active_first + n - 1, z, x);
          const Complex yn = (boundary_rhsn - eqn[0] * ynm3 - eqn[1] * ynm2 - eqn[2] * ynm1) / eqn[3];
          const Complex ynp1 = (upper_ghost(line) - Complex(eqnp1[0], 0.0) * ynm3 -
                                Complex(eqnp1[1], 0.0) * ynm2 - Complex(eqnp1[2], 0.0) * ynm1 -
                                Complex(eqnp1[3], 0.0) * yn) /
                               eqnp1[4];
          component(active_first + n, z, x) = yn;
          component(active_first + n + 1, z, x) = ynp1;
        });
    fence("dns_yline_reconstruct_upper_boundary");
  }
}

void DnsImplicitYLineStep::prepare(const DnsState& state, DnsImplicitYLineStepConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  if (cfg.components.empty()) {
    throw std::runtime_error("DnsImplicitYLineStep requires at least one component");
  }
  for (std::size_t i = 0; i < cfg.components.size(); ++i) {
    [[maybe_unused]] const auto component_origin = state.index(cfg.components[i], 0, 0);
    for (std::size_t j = i + 1; j < cfg.components.size(); ++j) {
      if (cfg.components[i] == cfg.components[j]) {
        throw std::runtime_error("DnsImplicitYLineStep component list contains duplicates");
      }
    }
  }

  grid_ = grid;
  cfg_ = std::move(cfg);
  solvers_.clear();
  solvers_.reserve(cfg_.components.size());
  for (const int component : cfg_.components) {
    DnsYLineComponentConfig component_cfg;
    component_cfg.component = component;
    component_cfg.npy = cfg_.npy;
    component_cfg.ipy = cfg_.ipy;
    component_cfg.pass_counts = cfg_.pass_counts;
    component_cfg.exchange_mode = cfg_.exchange_mode;
    component_cfg.comm_y = cfg_.comm_y;
    component_cfg.fill_interface_ghosts = cfg_.fill_interface_ghosts;
    auto solver = std::make_unique<DnsYLineComponentSolver>();
    solver->prepare(state, std::move(component_cfg));
    solvers_.push_back(std::move(solver));
  }
}

void DnsImplicitYLineStep::copy_coefficients_from_host(int component,
                                                       std::span<const Complex> ds,
                                                       std::span<const Complex> dl,
                                                       std::span<const Complex> d,
                                                       std::span<const Complex> du,
                                                       std::span<const Complex> dw) {
  component_solver(component).copy_coefficients_from_host(ds, dl, d, du, dw);
}

void DnsImplicitYLineStep::assemble_velocity_coefficients(int component,
                                                          const DeviceVector<double>& k2,
                                                          const DeviceVector<double>& derivatives,
                                                          int active_y_global_first,
                                                          int global_y_count,
                                                          double lambda,
                                                          double viscosity,
                                                          bool biharmonic) {
  component_solver(component).assemble_velocity_coefficients(k2,
                                                             derivatives,
                                                             active_y_global_first,
                                                             global_y_count,
                                                             lambda,
                                                             viscosity,
                                                             biharmonic);
}

void DnsImplicitYLineStep::copy_boundary_data_from_host(int component, const DnsYLineBoundaryData& boundary) {
  component_solver(component).copy_boundary_data_from_host(boundary);
}

void DnsImplicitYLineStep::solve(DnsState& state) {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error("DnsImplicitYLineStep::solve grid mismatch");
  }
  for (auto& solver : solvers_) {
    solver->solve(state);
  }
}

DnsYLineComponentSolver& DnsImplicitYLineStep::component_solver(int component) {
  return *solvers_[static_cast<std::size_t>(solver_index(component))];
}

const DnsYLineComponentSolver& DnsImplicitYLineStep::component_solver(int component) const {
  return *solvers_[static_cast<std::size_t>(solver_index(component))];
}

int DnsImplicitYLineStep::solver_index(int component) const {
  for (std::size_t i = 0; i < cfg_.components.size(); ++i) {
    if (cfg_.components[i] == component) return static_cast<int>(i);
  }
  throw std::runtime_error("DnsImplicitYLineStep component was not prepared");
}

} // namespace channel
