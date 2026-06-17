#include "channel/yline.hpp"

#include "channel/runtime.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace channel {

namespace {

KOKKOS_INLINE_FUNCTION int penta_index(int row, int line, int batch_count) {
  return line + row * batch_count;
}

KOKKOS_INLINE_FUNCTION int endpoint_slot_for_row(int coupled_row,
                                                  int local_n,
                                                  bool has_left_interface,
                                                  bool has_right_interface,
                                                  int left_count) {
  if (has_left_interface) {
    if (coupled_row == 0) return 0;
    if (coupled_row == 1) return 1;
  }
  if (has_right_interface) {
    if (coupled_row == local_n - 2) return left_count;
    if (coupled_row == local_n - 1) return left_count + 1;
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

  solve_batched_pentadiagonal(batch_ds_, batch_dl_, batch_d_, batch_du_, batch_dw_, batch_x_, interior, reduced_batch,
                              "endpoint_schur_interior_responses");

  auto leaf_rows = leaf_rows_.view();
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

  schur_.solve_from_leaf_rows(leaf_rows_, leaf_values_);

  auto leaf_values = leaf_values_.view();
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

void DnsYLineComponentSolver::prepare(const DnsState& state, DnsYLineComponentConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  [[maybe_unused]] const auto component_origin = state.index(cfg.component, 0, 0);
  if (grid.line_count() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsYLineComponentSolver line count exceeds int range");
  }
  if (grid.values_per_component() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsYLineComponentSolver component size exceeds int range");
  }
  cfg_ = std::move(cfg);
  grid_ = grid;

  const auto count = grid_.values_per_component();
  ds_.resize("dns_yline_ds", count);
  dl_.resize("dns_yline_dl", count);
  d_.resize("dns_yline_d", count);
  du_.resize("dns_yline_du", count);
  dw_.resize("dns_yline_dw", count);
  rhs_.resize("dns_yline_rhs", count);

  EndpointSchurYLineConfig endpoint_cfg;
  endpoint_cfg.local_n = grid_.ny;
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
  const auto expected = grid_.values_per_component();
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

void DnsYLineComponentSolver::solve(DnsState& state) {
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error("DnsYLineComponentSolver::solve grid mismatch");
  }
  auto component = state.component_view(cfg_.component);
  auto rhs = rhs_.view();
  Kokkos::parallel_for(
      "dns_yline_pack_component_rhs",
      Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(grid_.values_per_component())),
      KOKKOS_LAMBDA(const int i) { rhs(i) = component(i); });
  fence("dns_yline_pack_component_rhs");

  endpoint_solver_.solve(ds_, dl_, d_, du_, dw_, rhs_);

  Kokkos::parallel_for(
      "dns_yline_unpack_component_solution",
      Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<int>(grid_.values_per_component())),
      KOKKOS_LAMBDA(const int i) { component(i) = rhs(i); });
  fence("dns_yline_unpack_component_solution");
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
