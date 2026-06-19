#include "channel/dns_solver.hpp"
#include "channel/decomposition.hpp"
#include "channel/ffts.hpp"
#include "channel/input.hpp"
#include "channel/io.hpp"
#include "channel/runge_kutta.hpp"
#include "channel/runtime.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <fftw3.h>
#include <mpi.h>

namespace {

int dns_index(int y, int line, int batch) {
  return line + y * batch;
}

double sum_allreduce(double local) {
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

double max_allreduce(double local) {
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global;
}

double magnitude(const channel::Complex& value) {
  return std::hypot(static_cast<double>(value.real()), static_cast<double>(value.imag()));
}

bool fft_fit(int value) {
  while ((value % 2) == 0) value /= 2;
  return value == 1 || value == 3;
}

int fft_fit_at_least(int value) {
  while (!fft_fit(value)) ++value;
  return value;
}

std::vector<channel::Complex> extract_slice(const channel::DnsState& state,
                                            const channel::LineRange& local_y,
                                            int component) {
  const auto line_count = state.line_count();
  std::vector<channel::Complex> slice(static_cast<std::size_t>(local_y.count) * line_count,
                                      channel::Complex(0.0, 0.0));
  const auto source = state.component_host(component);
  for (int y = 0; y < local_y.count; ++y) {
    const int source_y = local_y.first + y - state.grid().y_first;
    const auto source_start = static_cast<std::size_t>(source_y) * line_count;
    const auto target_start = static_cast<std::size_t>(y) * line_count;
    std::copy(source.begin() + source_start,
              source.begin() + source_start + line_count,
              slice.begin() + target_start);
  }
  return slice;
}

std::vector<channel::Complex> extract_slice_x_window(const channel::DnsState& state,
                                                     const channel::LineRange& local_y,
                                                     int x_first,
                                                     int x_count,
                                                     int component) {
  const auto source = state.component_host(component);
  const int source_nx = state.grid().nx;
  const int nz = state.grid().nz;
  const int source_lines = static_cast<int>(state.line_count());
  const int target_lines = x_count * nz;
  std::vector<channel::Complex> slice(static_cast<std::size_t>(local_y.count) *
                                          static_cast<std::size_t>(target_lines),
                                      channel::Complex(0.0, 0.0));
  for (int y = 0; y < local_y.count; ++y) {
    const int source_y = local_y.first + y - state.grid().y_first;
    for (int z = 0; z < nz; ++z) {
      for (int x = 0; x < x_count; ++x) {
        const auto source_index =
            static_cast<std::size_t>(source_y * source_lines + z * source_nx + x_first + x);
        const auto target_index = static_cast<std::size_t>(y * target_lines + z * x_count + x);
        slice[target_index] = source[source_index];
      }
    }
  }
  return slice;
}

std::vector<channel::Complex> extract_storage_window(const channel::DnsState& state,
                                                     int y_first,
                                                     int y_count,
                                                     int component) {
  const auto line_count = state.line_count();
  std::vector<channel::Complex> slice(static_cast<std::size_t>(y_count) * line_count,
                                      channel::Complex(0.0, 0.0));
  const auto source = state.component_host(component);
  for (int y = 0; y < y_count; ++y) {
    const int source_y = y_first + y - state.grid().y_first;
    const auto source_start = static_cast<std::size_t>(source_y) * line_count;
    const auto target_start = static_cast<std::size_t>(y) * line_count;
    std::copy(source.begin() + source_start,
              source.begin() + source_start + line_count,
              slice.begin() + target_start);
  }
  return slice;
}

std::vector<channel::Complex> extract_storage_x_window(const channel::DnsState& state,
                                                       int y_first,
                                                       int y_count,
                                                       int x_first,
                                                       int x_count,
                                                       int component) {
  const auto source = state.component_host(component);
  const int source_nx = state.grid().nx;
  const int nz = state.grid().nz;
  const int source_lines = static_cast<int>(state.line_count());
  const int target_lines = x_count * nz;
  std::vector<channel::Complex> slice(static_cast<std::size_t>(y_count) *
                                          static_cast<std::size_t>(target_lines),
                                      channel::Complex(0.0, 0.0));
  for (int y = 0; y < y_count; ++y) {
    const int source_y = y_first + y - state.grid().y_first;
    for (int z = 0; z < nz; ++z) {
      for (int x = 0; x < x_count; ++x) {
        const auto source_index =
            static_cast<std::size_t>(source_y * source_lines + z * source_nx + x_first + x);
        const auto target_index = static_cast<std::size_t>(y * target_lines + z * x_count + x);
        slice[target_index] = source[source_index];
      }
    }
  }
  return slice;
}

void ensure_regression_input_match(const channel::ChannelInput& cfg,
                                   const channel::DnsState& start,
                                   const channel::DnsState& expected) {
  if (cfg.mesh.nx + 1 != start.grid().nx || cfg.mesh.ny + 3 != start.grid().ny ||
      2 * cfg.mesh.nz + 1 != start.grid().nz) {
    throw std::runtime_error("Regression input/fixture grid mismatch for start field");
  }
  if (cfg.mesh.nx + 1 != expected.grid().nx || cfg.mesh.ny + 3 != expected.grid().ny ||
      2 * cfg.mesh.nz + 1 != expected.grid().nz) {
    throw std::runtime_error("Regression input/fixture grid mismatch for end field");
  }
  if (start.grid().components != 3) {
    throw std::runtime_error("Velocity regression start fixture must contain exactly three velocity components");
  }
  if (expected.grid().components != 3) {
    throw std::runtime_error("Velocity regression expected fixture must contain exactly three velocity components");
  }
  if (cfg.scalars.nphi != 0) {
    throw std::runtime_error("Velocity regression requires nphi=0");
  }
}

std::vector<double> y_coordinates(const channel::ChannelInput& input);
int derivative_index(int y, int order, int offset);

std::vector<double> benchmark_derivatives(int ny) {
  std::vector<double> derivatives(static_cast<std::size_t>(ny) *
                                      channel::DnsNonlinearVelocityRhsStage::derivative_orders *
                                      channel::DnsNonlinearVelocityRhsStage::derivative_stencil,
                                  0.0);
  auto index = [](int y, int order, int offset) {
    return (y * channel::DnsNonlinearVelocityRhsStage::derivative_orders + order) *
               channel::DnsNonlinearVelocityRhsStage::derivative_stencil +
           (offset + 2);
  };
  for (int y = 0; y < ny; ++y) {
    derivatives[static_cast<std::size_t>(index(y, 0, 0))] = 1.0;
    derivatives[static_cast<std::size_t>(index(y, 1, -1))] = -0.5;
    derivatives[static_cast<std::size_t>(index(y, 1, 1))] = 0.5;
    derivatives[static_cast<std::size_t>(index(y, 2, -1))] = 1.0;
    derivatives[static_cast<std::size_t>(index(y, 2, 0))] = -2.0;
    derivatives[static_cast<std::size_t>(index(y, 2, 1))] = 1.0;
    derivatives[static_cast<std::size_t>(index(y, 3, -2))] = 1.0;
    derivatives[static_cast<std::size_t>(index(y, 3, -1))] = -4.0;
    derivatives[static_cast<std::size_t>(index(y, 3, 0))] = 6.0;
    derivatives[static_cast<std::size_t>(index(y, 3, 1))] = -4.0;
    derivatives[static_cast<std::size_t>(index(y, 3, 2))] = 1.0;
  }
  return derivatives;
}

std::array<double, 5> solve_5x5(std::array<std::array<double, 5>, 5> a, std::array<double, 5> b) {
  for (int i = 4; i >= 1; --i) {
    if (a[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] == 0.0) {
      throw std::runtime_error("singular 5x5 derivative setup matrix");
    }
    const double pivot = 1.0 / a[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    a[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = pivot;
    for (int j = 0; j < i; ++j) {
      a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *= pivot;
    }
    for (int k = 0; k < i; ++k) {
      const double factor = a[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
      for (int j = 0; j < i; ++j) {
        a[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)] -=
            factor * a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      }
    }
  }
  if (a[0][0] == 0.0) {
    throw std::runtime_error("singular 5x5 derivative setup matrix");
  }
  a[0][0] = 1.0 / a[0][0];

  std::array<double, 5> x{};
  x[4] = b[4] * a[4][4];
  for (int i = 3; i >= 0; --i) {
    double sum = 0.0;
    for (int j = i + 1; j < 5; ++j) {
      sum += a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
             x[static_cast<std::size_t>(j)];
    }
    x[static_cast<std::size_t>(i)] =
        (b[static_cast<std::size_t>(i)] - sum) * a[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
  }
  for (int i = 1; i < 5; ++i) {
    double sum = 0.0;
    for (int j = 0; j < i; ++j) {
      sum += a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
             x[static_cast<std::size_t>(j)];
    }
    x[static_cast<std::size_t>(i)] -= sum;
  }
  return x;
}

std::array<std::array<double, 5>, 5> interpolation_matrix(const std::vector<double>& y, int iy) {
  const auto y_at = [&](int global_y) -> double {
    return y[static_cast<std::size_t>(global_y + 1)];
  };
  std::array<std::array<double, 5>, 5> m{};
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      m[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
          std::pow(y_at(iy - 2 + j) - y_at(iy), 4.0 - static_cast<double>(i));
    }
  }
  return m;
}

std::vector<double> compact_derivatives(const channel::ChannelInput& input, const channel::LineRange& local_y) {
  const auto y = y_coordinates(input);
  const auto y_at = [&](int global_y) -> double {
    return y[static_cast<std::size_t>(global_y + 1)];
  };
  std::vector<double> derivatives(static_cast<std::size_t>(local_y.count) *
                                      channel::DnsNonlinearVelocityRhsStage::derivative_orders *
                                      channel::DnsNonlinearVelocityRhsStage::derivative_stencil,
                                  0.0);
  auto set_derivative = [&](int local_row, int order, const std::array<double, 5>& coeffs) {
    for (int offset = -2; offset <= 2; ++offset) {
      derivatives[static_cast<std::size_t>(derivative_index(local_row, order, offset))] =
          coeffs[static_cast<std::size_t>(offset + 2)];
    }
  };

  for (int local_row = 0; local_row < local_y.count; ++local_row) {
    const int iy = local_y.first + local_row;
    if (iy < 1 || iy > input.mesh.ny - 1) {
      continue;
    }

    auto m = interpolation_matrix(y, iy);
    std::array<double, 5> t{};
    t[0] = 24.0;
    const auto d4 = solve_5x5(m, t);
    set_derivative(local_row, 3, d4);

    std::array<std::array<double, 5>, 5> d0_matrix{};
    for (int i = 0; i < 5; ++i) {
      const double factor = (5.0 - static_cast<double>(i)) * (6.0 - static_cast<double>(i)) *
                            (7.0 - static_cast<double>(i)) * (8.0 - static_cast<double>(i));
      for (int j = 0; j < 5; ++j) {
        d0_matrix[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
            factor * std::pow(y_at(iy - 2 + j) - y_at(iy), 4.0 - static_cast<double>(i));
      }
    }
    for (int i = 0; i < 5; ++i) {
      t[static_cast<std::size_t>(i)] = 0.0;
      for (int j = 0; j < 5; ++j) {
        t[static_cast<std::size_t>(i)] +=
            d4[static_cast<std::size_t>(j)] *
            std::pow(y_at(iy - 2 + j) - y_at(iy), 8.0 - static_cast<double>(i));
      }
    }
    const auto d0 = solve_5x5(d0_matrix, t);
    set_derivative(local_row, 0, d0);

    m = interpolation_matrix(y, iy);
    t = {};
    for (int i = 0; i <= 2; ++i) {
      for (int j = 0; j < 5; ++j) {
        t[static_cast<std::size_t>(i)] +=
            d0[static_cast<std::size_t>(j)] * (4.0 - static_cast<double>(i)) *
            (3.0 - static_cast<double>(i)) *
            std::pow(y_at(iy - 2 + j) - y_at(iy), 2.0 - static_cast<double>(i));
      }
    }
    const auto d2 = solve_5x5(m, t);
    set_derivative(local_row, 2, d2);

    m = interpolation_matrix(y, iy);
    t = {};
    for (int i = 0; i <= 3; ++i) {
      for (int j = 0; j < 5; ++j) {
        t[static_cast<std::size_t>(i)] +=
            d0[static_cast<std::size_t>(j)] * (4.0 - static_cast<double>(i)) *
            std::pow(y_at(iy - 2 + j) - y_at(iy), 3.0 - static_cast<double>(i));
      }
    }
    const auto d1 = solve_5x5(m, t);
    set_derivative(local_row, 1, d1);
  }
  static bool printed_debug_derivatives = false;
  const char* debug_env = std::getenv("CHANNEL_DEBUG_TIMESTEP");
  if (!printed_debug_derivatives && debug_env != nullptr && std::string(debug_env).find("ghost") != std::string::npos &&
      local_y.first == 1 && local_y.count == input.mesh.ny - 1) {
    printed_debug_derivatives = true;
    const int step = std::max(1, (input.mesh.ny - 2) / 2);
    for (int iy = 1; iy <= input.mesh.ny - 1; iy += step) {
      const int local_row = iy - local_y.first;
      for (int order : {0, 2, 3}) {
        std::cout << "CPP_DER iy=" << iy << " d" << order << "=";
        for (int offset = -2; offset <= 2; ++offset) {
          std::cout << " " << std::scientific << std::setprecision(16)
                    << derivatives[static_cast<std::size_t>(derivative_index(local_row, order, offset))];
        }
        std::cout << '\n';
      }
    }
  }
  return derivatives;
}

void make_line_wavenumbers(const channel::ChannelInput& input,
                           const channel::DnsGrid& grid,
                           int global_x_first,
                           std::vector<channel::Complex>& ialfa,
                           std::vector<channel::Complex>& ibeta,
                           std::vector<double>& k2) {
  const int nx = grid.nx;
  const int nz = grid.nz;
  const int lines = nx * nz;
  ialfa.assign(static_cast<std::size_t>(lines), channel::Complex(0.0, 0.0));
  ibeta.assign(static_cast<std::size_t>(lines), channel::Complex(0.0, 0.0));
  k2.assign(static_cast<std::size_t>(lines), 0.0);
  for (int iz = 0; iz < nz; ++iz) {
    const int signed_iz = iz <= input.mesh.nz ? iz : iz - (2 * input.mesh.nz + 1);
    for (int ix = 0; ix < nx; ++ix) {
      const int p = iz * nx + ix;
      const int global_ix = global_x_first + ix;
      ialfa[static_cast<std::size_t>(p)] = channel::Complex(0.0, input.mesh.alfa0 * global_ix);
      ibeta[static_cast<std::size_t>(p)] = channel::Complex(0.0, input.mesh.beta0 * signed_iz);
      k2[static_cast<std::size_t>(p)] =
          (input.mesh.alfa0 * global_ix) * (input.mesh.alfa0 * global_ix) +
          (input.mesh.beta0 * signed_iz) * (input.mesh.beta0 * signed_iz);
    }
  }
}

int derivative_index(int y, int order, int offset) {
  return (y * channel::DnsNonlinearVelocityRhsStage::derivative_orders + order) *
             channel::DnsNonlinearVelocityRhsStage::derivative_stencil +
         (offset + 2);
}

double derivative_value(const std::vector<double>& derivatives, int y, int order, int offset) {
  return derivatives[static_cast<std::size_t>(derivative_index(y, order, offset))];
}

struct BandedCoefficients {
  std::vector<channel::Complex> ds;
  std::vector<channel::Complex> dl;
  std::vector<channel::Complex> d;
  std::vector<channel::Complex> du;
  std::vector<channel::Complex> dw;
};

std::array<std::array<double, 5>, 5> boundary_matrix(const std::vector<double>& y,
                                                      int first_point,
                                                      int center) {
  const auto y_at = [&](int global_y) -> double {
    return y[static_cast<std::size_t>(global_y + 1)];
  };
  std::array<std::array<double, 5>, 5> m{};
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      m[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
          std::pow(y_at(first_point + j) - y_at(center), 4.0 - static_cast<double>(i));
    }
  }
  return m;
}

std::array<double, 5> boundary_derivative_coefficients(const std::vector<double>& y,
                                                       int first_point,
                                                       int center,
                                                       int derivative_order) {
  std::array<double, 5> t{};
  if (derivative_order == 1) {
    t[3] = 1.0;
  } else if (derivative_order == 2) {
    t[2] = 2.0;
  } else {
    throw std::runtime_error("unsupported boundary derivative order");
  }
  return solve_5x5(boundary_matrix(y, first_point, center), t);
}

std::array<double, 5> compact_row_d4(const channel::ChannelInput& input, const std::vector<double>& y, int iy) {
  std::array<double, 5> t{};
  t[0] = 24.0;
  return solve_5x5(interpolation_matrix(y, iy), t);
}

struct VelocityBoundaryData {
  channel::DnsYLineBoundaryData eta;
  channel::DnsYLineBoundaryData v;
  std::vector<channel::Complex> zero;
};

struct CompactDerivativeBoundaryData {
  channel::DnsYLineBoundaryData boundary;
  std::vector<channel::Complex> zero;
  std::array<double, 5> lower_boundary_rhs_coeff{};
  std::array<double, 5> lower_ghost_rhs_coeff{};
  std::array<double, 5> upper_boundary_rhs_coeff{};
  std::array<double, 5> upper_ghost_rhs_coeff{};
};

VelocityBoundaryData velocity_boundary_data(const channel::ChannelInput& input,
                                            const channel::LineRange& local_y,
                                            int line_count) {
  VelocityBoundaryData data;
  data.zero.assign(static_cast<std::size_t>(line_count), channel::Complex(0.0, 0.0));
  const bool has_lower = local_y.first == 1;
  const bool has_upper = (local_y.first + local_y.count - 1) == input.mesh.ny - 1;
  data.eta.lower_enabled = has_lower;
  data.eta.upper_enabled = has_upper;
  data.v.lower_enabled = has_lower;
  data.v.upper_enabled = has_upper;
  const auto y = y_coordinates(input);

  std::array<double, 5> d040{};
  d040[1] = 1.0;
  std::array<double, 5> d04n{};
  d04n[3] = 1.0;
  const auto d140 = boundary_derivative_coefficients(y, -1, 0, 1);
  const auto d14n = boundary_derivative_coefficients(y, input.mesh.ny - 3, input.mesh.ny, 1);
  const auto d4_lower = compact_row_d4(input, y, 1);
  const auto d4_upper = compact_row_d4(input, y, input.mesh.ny - 1);

  data.v.lower_boundary_eq = d040;
  data.v.lower_ghost_eq = d140;
  data.v.upper_boundary_eq = d04n;
  data.v.upper_ghost_eq = d14n;
  data.eta.lower_boundary_eq = d040;
  data.eta.lower_ghost_eq = d4_lower;
  data.eta.upper_boundary_eq = d04n;
  data.eta.upper_ghost_eq = d4_upper;

  const std::span<const channel::Complex> zero_span(data.zero.data(), data.zero.size());
  data.v.lower_boundary_rhs = zero_span;
  data.v.lower_ghost_rhs = zero_span;
  data.v.upper_boundary_rhs = zero_span;
  data.v.upper_ghost_rhs = zero_span;
  data.eta.lower_boundary_rhs = zero_span;
  data.eta.lower_ghost_rhs = zero_span;
  data.eta.upper_boundary_rhs = zero_span;
  data.eta.upper_ghost_rhs = zero_span;
  return data;
}

CompactDerivativeBoundaryData compact_derivative_boundary_data(const channel::ChannelInput& input,
                                                               const channel::LineRange& local_y,
                                                               int line_count) {
  CompactDerivativeBoundaryData data;
  data.zero.assign(static_cast<std::size_t>(line_count), channel::Complex(0.0, 0.0));
  const bool has_lower = local_y.first == 1;
  const bool has_upper = (local_y.first + local_y.count - 1) == input.mesh.ny - 1;
  data.boundary.lower_enabled = has_lower;
  data.boundary.upper_enabled = has_upper;

  data.boundary.lower_ghost_eq = {1.0, 0.0, 0.0, 0.0, 0.0};
  data.boundary.lower_boundary_eq = {0.0, 1.0, 0.0, 0.0, 0.0};
  data.boundary.upper_boundary_eq = {0.0, 0.0, 0.0, 1.0, 0.0};
  data.boundary.upper_ghost_eq = {0.0, 0.0, 0.0, 0.0, 1.0};

  const auto y = y_coordinates(input);
  data.lower_boundary_rhs_coeff = boundary_derivative_coefficients(y, -1, 0, 1);
  data.lower_ghost_rhs_coeff = boundary_derivative_coefficients(y, -1, -1, 1);
  data.upper_boundary_rhs_coeff = boundary_derivative_coefficients(y, input.mesh.ny - 3, input.mesh.ny, 1);
  data.upper_ghost_rhs_coeff = boundary_derivative_coefficients(y, input.mesh.ny - 3, input.mesh.ny + 1, 1);

  const std::span<const channel::Complex> zero_span(data.zero.data(), data.zero.size());
  data.boundary.lower_boundary_rhs = zero_span;
  data.boundary.lower_ghost_rhs = zero_span;
  data.boundary.upper_boundary_rhs = zero_span;
  data.boundary.upper_ghost_rhs = zero_span;
  return data;
}

BandedCoefficients assemble_velocity_operator(const channel::DnsState& state,
                                               const channel::LineRange& local_y,
                                               int global_ny,
                                               const std::vector<double>& k2,
                                               const std::vector<double>& derivatives,
                                               double lambda,
                                               double viscosity,
                                               bool biharmonic) {
  const auto& grid = state.grid();
  const int lines = static_cast<int>(grid.line_count());
  const int active_count = grid.active_y_count == 0 ? grid.ny : grid.active_y_count;
  BandedCoefficients coeffs;
  coeffs.ds.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.dl.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.d.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.du.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.dw.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));

  for (int y = 0; y < active_count; ++y) {
    const int global_y = local_y.first + y;
    for (int line = 0; line < lines; ++line) {
      const int p = dns_index(y, line, lines);
      if (global_y == 0 || global_y == global_ny - 1) {
        coeffs.d[static_cast<std::size_t>(p)] = channel::Complex(1.0, 0.0);
        continue;
      }

      std::array<channel::Complex, 5> row{};
      const double k2_line = k2[static_cast<std::size_t>(line)];
      for (int offset = -2; offset <= 2; ++offset) {
        const double d0 = derivative_value(derivatives, y, 0, offset);
        const double d2 = derivative_value(derivatives, y, 2, offset);
        const double d4 = derivative_value(derivatives, y, 3, offset);
        double value = 0.0;
        if (biharmonic) {
          value = lambda * (d2 - k2_line * d0) -
                  viscosity * (d4 - 2.0 * k2_line * d2 + k2_line * k2_line * d0);
        } else {
          value = lambda * d0 - viscosity * (d2 - k2_line * d0);
        }
        row[static_cast<std::size_t>(offset + 2)] = channel::Complex(value, 0.0);
      }
      coeffs.ds[static_cast<std::size_t>(p)] = row[0];
      coeffs.dl[static_cast<std::size_t>(p)] = row[1];
      coeffs.d[static_cast<std::size_t>(p)] = row[2];
      coeffs.du[static_cast<std::size_t>(p)] = row[3];
      coeffs.dw[static_cast<std::size_t>(p)] = row[4];
    }
  }
  return coeffs;
}

BandedCoefficients assemble_compact_derivative_operator(const channel::DnsState& state,
                                                        const std::vector<double>& derivatives) {
  const auto& grid = state.grid();
  const int lines = static_cast<int>(grid.line_count());
  const int active_count = grid.active_y_count == 0 ? grid.ny : grid.active_y_count;
  BandedCoefficients coeffs;
  coeffs.ds.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.dl.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.d.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.du.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));
  coeffs.dw.assign(grid.active_values_per_component(), channel::Complex(0.0, 0.0));

  for (int y = 0; y < active_count; ++y) {
    for (int line = 0; line < lines; ++line) {
      const int p = dns_index(y, line, lines);
      coeffs.ds[static_cast<std::size_t>(p)] = channel::Complex(derivative_value(derivatives, y, 0, -2), 0.0);
      coeffs.dl[static_cast<std::size_t>(p)] = channel::Complex(derivative_value(derivatives, y, 0, -1), 0.0);
      coeffs.d[static_cast<std::size_t>(p)] = channel::Complex(derivative_value(derivatives, y, 0, 0), 0.0);
      coeffs.du[static_cast<std::size_t>(p)] = channel::Complex(derivative_value(derivatives, y, 0, 1), 0.0);
      coeffs.dw[static_cast<std::size_t>(p)] = channel::Complex(derivative_value(derivatives, y, 0, 2), 0.0);
    }
  }
  return coeffs;
}

void add_velocity_yline(channel::DnsRungeKuttaTimestepper& stepper,
                        const channel::ChannelInput& input,
                        const channel::DnsState& state,
                        const channel::LineRange& local_y,
                        int global_ny,
                        const std::vector<double>& k2,
                        const std::vector<double>& derivatives,
                        const std::array<channel::RungeKuttaWeights, 3>& weights,
                        double dt,
                        double viscosity) {
  for (int stage = 0; stage < static_cast<int>(weights.size()); ++stage) {
    const double lambda = weights[static_cast<std::size_t>(stage)][0] / dt;
    auto eta = assemble_velocity_operator(state, local_y, global_ny, k2, derivatives, lambda, viscosity, false);
    auto v = assemble_velocity_operator(state, local_y, global_ny, k2, derivatives, lambda, viscosity, true);
    stepper.copy_yline_coefficients_from_host(stage, 0, eta.ds, eta.dl, eta.d, eta.du, eta.dw);
    stepper.copy_yline_coefficients_from_host(stage, 1, v.ds, v.dl, v.d, v.du, v.dw);
    const auto boundaries = velocity_boundary_data(input, local_y, static_cast<int>(state.line_count()));
    stepper.copy_yline_boundary_data_from_host(stage, 0, boundaries.eta);
    stepper.copy_yline_boundary_data_from_host(stage, 1, boundaries.v);
  }
}

channel::CompactBoundaryRows compact_boundaries_from_yline(const channel::DnsYLineBoundaryData& boundary) {
  channel::CompactBoundaryRows rows;
  rows.lower = boundary.lower_boundary_eq;
  rows.lower_ghost = boundary.lower_ghost_eq;
  rows.upper = boundary.upper_boundary_eq;
  rows.upper_ghost = boundary.upper_ghost_eq;
  rows.rhs_lower = channel::Complex(0.0, 0.0);
  rows.rhs_lower_ghost = channel::Complex(0.0, 0.0);
  rows.rhs_upper = channel::Complex(0.0, 0.0);
  rows.rhs_upper_ghost = channel::Complex(0.0, 0.0);
  return rows;
}

std::vector<double> assemble_velocity_mean_correction_matrix(const channel::ChannelInput& input,
                                                             double lambda,
                                                             double viscosity) {
  const channel::LineRange full_y{1, input.mesh.ny - 1};
  const auto derivatives = compact_derivatives(input, full_y);
  std::vector<double> matrix(static_cast<std::size_t>(input.mesh.ny + 1) * 5, 0.0);
  for (int iy = 1; iy <= input.mesh.ny - 1; ++iy) {
    const int local_y = iy - 1;
    for (int offset = -2; offset <= 2; ++offset) {
      const double d0 = derivative_value(derivatives, local_y, 0, offset);
      const double d2 = derivative_value(derivatives, local_y, 2, offset);
      matrix[static_cast<std::size_t>(iy * 5 + offset + 2)] = lambda * d0 - viscosity * d2;
    }
  }
  return matrix;
}

void add_velocity_mean_correction_matrix(channel::DnsRungeKuttaTimestepper& stepper,
                                         const channel::ChannelInput& input,
                                         const std::array<channel::RungeKuttaWeights, 3>& weights,
                                         double dt,
                                         double viscosity) {
  for (int stage = 0; stage < static_cast<int>(weights.size()); ++stage) {
    const double lambda = weights[static_cast<std::size_t>(stage)][0] / dt;
    const auto matrix = assemble_velocity_mean_correction_matrix(input, lambda, viscosity);
    stepper.copy_velocity_mean_correction_matrix_from_host(stage, matrix);
  }
}

std::vector<double> y_coordinates(const channel::ChannelInput& input) {
  std::vector<double> y(static_cast<std::size_t>(input.mesh.ny + 3), 0.0);
  const auto at = [&](int iy) -> std::size_t { return static_cast<std::size_t>(iy + 1); };
  for (int iy = -1; iy <= input.mesh.ny + 1; ++iy) {
    const double eta = 2.0 * static_cast<double>(iy) / static_cast<double>(input.mesh.ny) - 1.0;
    y[at(iy)] = input.mesh.ymin +
                0.5 * (input.mesh.ymax - input.mesh.ymin) *
                    (std::tanh(input.mesh.stretching * eta) / std::tanh(input.mesh.stretching) + 1.0);
  }
  return y;
}

double choose_regression_dt(const channel::ChannelInput& input,
                            const channel::DnsState& state,
                            const channel::LineRange& local_y,
                            const channel::Decomposition& decomp) {
  if (input.timestepping.dt > 0.0 || input.timestepping.cflmax <= 0.0) {
    return input.timestepping.dt;
  }

  // Legacy Fortran declares PI as double but initializes it from a default-real literal.
  constexpr double pi = 3.1415927410125732421875;
  const int nxd = fft_fit_at_least(3 * (input.mesh.nx + 1) / 2);
  const int nzd = fft_fit_at_least(3 * input.mesh.nz);
  const int physical_x = 2 * nxd;
  const int real_x_stride = 2 * (nxd + 1);
  const int padded_x_half = nxd + 1;
  const double dx = pi / (input.mesh.alfa0 * static_cast<double>(nxd));
  const double dz = 2.0 * pi / (input.mesh.beta0 * static_cast<double>(nzd));
  const auto y = y_coordinates(input);
  const auto y_at = [&](int iy) -> double { return y[static_cast<std::size_t>(iy + 1)]; };
  std::vector<double> dy_host(static_cast<std::size_t>(state.grid().ny), 0.0);
  for (int y_storage = 0; y_storage < state.grid().ny; ++y_storage) {
    const int y_global = state.grid().y_first + y_storage;
    if (y_global >= 1 && y_global <= input.mesh.ny - 1) {
      dy_host[static_cast<std::size_t>(y_storage)] = 0.5 * (y_at(y_global + 1) - y_at(y_global - 1));
    }
  }

  if (decomp.npxz > 1) {
    channel::DistributedDealiasedFft2DPlan cfl_fft;
    cfl_fft.configure(state.grid(),
                      {input.mesh.nx + 1,
                       decomp.nx0,
                       physical_x,
                       nzd,
                       decomp.npxz,
                       decomp.ipxz,
                       decomp.comm_x});
    channel::RealView3D u("regression_cfl_u", state.grid().ny, cfl_fft.z_count(), physical_x);
    channel::RealView3D v("regression_cfl_v", state.grid().ny, cfl_fft.z_count(), physical_x);
    channel::RealView3D w("regression_cfl_w", state.grid().ny, cfl_fft.z_count(), physical_x);
    cfl_fft.inverse_component_to_physical(state, 0, u, "regression_cfl_u");
    cfl_fft.inverse_component_to_physical(state, 1, v, "regression_cfl_v");
    cfl_fft.inverse_component_to_physical(state, 2, w, "regression_cfl_w");

    auto u_host = Kokkos::create_mirror_view(u);
    auto v_host = Kokkos::create_mirror_view(v);
    auto w_host = Kokkos::create_mirror_view(w);
    Kokkos::deep_copy(u_host, u);
    Kokkos::deep_copy(v_host, v);
    Kokkos::deep_copy(w_host, w);

    double local_cfl = 0.0;
    for (int y_storage = 0; y_storage < state.grid().ny; ++y_storage) {
      const int y_global = state.grid().y_first + y_storage;
      if (y_global < 1 || y_global > input.mesh.ny - 1) continue;
      for (int z = 0; z < cfl_fft.z_count(); ++z) {
        for (int x = 0; x < physical_x; ++x) {
          const double value =
              std::abs(u_host(y_storage, z, x)) / dx +
              std::abs(v_host(y_storage, z, x)) / dy_host[static_cast<std::size_t>(y_storage)] +
              std::abs(w_host(y_storage, z, x)) / dz;
          if (value > local_cfl) local_cfl = value;
        }
      }
    }

    (void)local_y;
    double global_cfl = 0.0;
    MPI_Allreduce(&local_cfl, &global_cfl, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (!(global_cfl > 0.0) || !std::isfinite(global_cfl)) {
      throw std::runtime_error("could not compute a finite adaptive distributed regression timestep");
    }
    return input.timestepping.cflmax / global_cfl;
  }

  const int storage_z = state.grid().nz;
  const int retained_z = (storage_z - 1) / 2;
  std::vector<std::complex<double>> vvdz(static_cast<std::size_t>(nzd) * static_cast<std::size_t>(state.grid().nx));
  std::vector<std::complex<double>> vvdx(static_cast<std::size_t>(padded_x_half) * static_cast<std::size_t>(nzd));
  std::vector<double> real_work(static_cast<std::size_t>(real_x_stride) * static_cast<std::size_t>(nzd));
  const int n_z[1] = {nzd};
  const int n_x[1] = {physical_x};
  const int inembed_z[1] = {nzd};
  const int inembed_x[1] = {padded_x_half};
  const int onembed_x[1] = {real_x_stride};
  fftw_plan z_ift = fftw_plan_many_dft(1,
                                       n_z,
                                       state.grid().nx,
                                       reinterpret_cast<fftw_complex*>(vvdz.data()),
                                       inembed_z,
                                       1,
                                       nzd,
                                       reinterpret_cast<fftw_complex*>(vvdz.data()),
                                       inembed_z,
                                       1,
                                       nzd,
                                       FFTW_BACKWARD,
                                       FFTW_PATIENT);
  fftw_plan x_rft = fftw_plan_many_dft_c2r(1,
                                           n_x,
                                           nzd,
                                           reinterpret_cast<fftw_complex*>(vvdx.data()),
                                           inembed_x,
                                           1,
                                           padded_x_half,
                                           real_work.data(),
                                           onembed_x,
                                           1,
                                           real_x_stride,
                                           FFTW_PATIENT);
  if (z_ift == nullptr || x_rft == nullptr) {
    throw std::runtime_error("could not create FFTW plans for regression CFL");
  }

  const auto transform_component = [&](int component_id) {
    const auto component = state.component_host(component_id);
    std::vector<double> physical(static_cast<std::size_t>(state.grid().ny) * static_cast<std::size_t>(nzd) *
                                     static_cast<std::size_t>(real_x_stride),
                                 0.0);
    for (int y_storage = 0; y_storage < state.grid().ny; ++y_storage) {
      std::fill(vvdz.begin(), vvdz.end(), std::complex<double>(0.0, 0.0));
      for (int x_storage = 0; x_storage < state.grid().nx; ++x_storage) {
        for (int z_storage = 0; z_storage < state.grid().nz; ++z_storage) {
          const int signed_z = z_storage <= retained_z ? z_storage : z_storage - storage_z;
          const int padded_z = signed_z >= 0 ? signed_z : signed_z + nzd;
          const auto source = static_cast<std::size_t>((y_storage * state.grid().nz + z_storage) * state.grid().nx +
                                                       x_storage);
          const auto target = static_cast<std::size_t>(padded_z + nzd * x_storage);
          vvdz[target] = std::complex<double>(component[source].real(), component[source].imag());
        }
      }
      fftw_execute_dft(z_ift, reinterpret_cast<fftw_complex*>(vvdz.data()), reinterpret_cast<fftw_complex*>(vvdz.data()));

      std::fill(vvdx.begin(), vvdx.end(), std::complex<double>(0.0, 0.0));
      for (int x_storage = 0; x_storage < state.grid().nx; ++x_storage) {
        for (int z_storage = 0; z_storage < nzd; ++z_storage) {
          vvdx[static_cast<std::size_t>(x_storage + padded_x_half * z_storage)] =
              vvdz[static_cast<std::size_t>(z_storage + nzd * x_storage)];
        }
      }
      std::fill(real_work.begin(), real_work.end(), 0.0);
      fftw_execute_dft_c2r(x_rft, reinterpret_cast<fftw_complex*>(vvdx.data()), real_work.data());

      const auto target = static_cast<std::size_t>(y_storage) * static_cast<std::size_t>(nzd) *
                          static_cast<std::size_t>(real_x_stride);
      std::copy(real_work.begin(), real_work.end(), physical.begin() + static_cast<std::ptrdiff_t>(target));
    }
    return physical;
  };

  const auto u_host = transform_component(0);
  const auto v_host = transform_component(1);
  const auto w_host = transform_component(2);
  fftw_destroy_plan(z_ift);
  fftw_destroy_plan(x_rft);

  const int y_first = state.grid().y_first;
  const int global_ny = input.mesh.ny;
  double local_cfl = 0.0;
  for (int y_storage = 0; y_storage < state.grid().ny; ++y_storage) {
    const int y_global = y_first + y_storage;
    if (y_global < 1 || y_global > global_ny - 1) continue;
    for (int z = 0; z < nzd; ++z) {
      for (int x = 0; x < physical_x; ++x) {
        const auto p = (static_cast<std::size_t>(y_storage) * static_cast<std::size_t>(nzd) +
                        static_cast<std::size_t>(z)) *
                           static_cast<std::size_t>(real_x_stride) +
                       static_cast<std::size_t>(x);
        const double value =
            std::abs(u_host[p]) / dx + std::abs(v_host[p]) / dy_host[static_cast<std::size_t>(y_storage)] +
            std::abs(w_host[p]) / dz;
        if (value > local_cfl) local_cfl = value;
      }
    }
  }

  const char* debug_env = std::getenv("CHANNEL_DEBUG_TIMESTEP");
  if (debug_env != nullptr && std::string(debug_env) == "ghosts") {
    double max_value = -1.0;
    int max_y = -999;
    int max_z = -999;
    int max_x = -999;
    double max_u = 0.0;
    double max_v = 0.0;
    double max_w = 0.0;
    double max_du = 0.0;
    double max_dv = 0.0;
    double max_dw = 0.0;
    const auto physical_index = [&](int y_storage, int z, int x) -> std::size_t {
      return (static_cast<std::size_t>(y_storage) * static_cast<std::size_t>(nzd) +
              static_cast<std::size_t>(z)) *
                 static_cast<std::size_t>(real_x_stride) +
             static_cast<std::size_t>(x);
    };
    for (int y_storage = 0; y_storage < state.grid().ny; ++y_storage) {
      const int y_global = y_first + y_storage;
      if (y_global < 1 || y_global > global_ny - 1) continue;
      for (int z = 0; z < nzd; ++z) {
        for (int x = 0; x < physical_x; ++x) {
          const auto p = physical_index(y_storage, z, x);
          const double value = std::abs(u_host[p]) / dx + std::abs(v_host[p]) / dy_host[y_storage] +
                               std::abs(w_host[p]) / dz;
          if (value > max_value) {
            max_value = value;
            max_y = y_global;
            max_z = z;
            max_x = x;
            max_u = u_host[p];
            max_v = v_host[p];
            max_w = w_host[p];
            max_du = std::abs(u_host[p]) / dx;
            max_dv = std::abs(v_host[p]) / dy_host[y_storage];
            max_dw = std::abs(w_host[p]) / dz;
          }
        }
      }
    }
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout << std::scientific << std::setprecision(16)
              << "CPP_CFL rank=" << rank
              << " y=" << max_y
              << " x=" << max_x
              << " z=" << max_z
              << " value=" << max_value
              << " u=" << max_u
              << " v=" << max_v
              << " w=" << max_w
              << " dx=" << dx
              << " dy=" << dy_host[static_cast<std::size_t>(max_y - y_first)]
              << " dz=" << dz
              << " du=" << max_du
              << " dv=" << max_dv
              << " dw=" << max_dw << '\n';
  }

  (void)local_y;

  double global_cfl = 0.0;
  MPI_Allreduce(&local_cfl, &global_cfl, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (!(global_cfl > 0.0) || !std::isfinite(global_cfl)) {
    throw std::runtime_error("could not compute a finite adaptive regression timestep");
  }
  return input.timestepping.cflmax / global_cfl;
}

int count_nonfinite_components(const channel::DnsState& state, int components) {
  int bad = 0;
  for (int component = 0; component < components; ++component) {
    const auto values = state.component_host(component);
    for (const auto& value : values) {
      if (!std::isfinite(static_cast<double>(value.real())) ||
          !std::isfinite(static_cast<double>(value.imag()))) {
        ++bad;
      }
    }
  }
  return bad;
}

int run_velocity_regression(const channel::Runtime& runtime,
                            const std::string& input_path,
                            const std::string& start_path,
                            const std::string& expected_path) {
  const auto input = channel::read_channel_input(input_path);
  const int requested_npy = input.parallel.npy_was_set ? input.parallel.npy : runtime.size();

  channel::DnsState start_field;
  const auto start_header = channel::read_legacy_fortran_full_spectrum_with_ghosts_for_tests(start_path, start_field);
  channel::DnsState expected_field;
  const auto expected_header = channel::read_legacy_fortran_full_spectrum_with_ghosts_for_tests(expected_path, expected_field);
  ensure_regression_input_match(input, start_field, expected_field);

  channel::Decomposition decomp;
  decomp.initialize({start_field.grid().nx,
                     input.mesh.nz,
                     input.mesh.ny,
                     fft_fit_at_least(3 * input.mesh.nz),
                     0,
                     false,
                     requested_npy},
                    channel::world_comm());

  constexpr int base_components = 3;
  constexpr int product_components_first = 5;
  constexpr int product_components_count = 6;
  constexpr int solver_components = 11;
  std::array<int, product_components_count> product_components{};
  for (int i = 0; i < product_components_count; ++i) {
    product_components[static_cast<std::size_t>(i)] = product_components_first + i;
  }

  auto local_y = channel::split_range(decomp.ipy, input.mesh.ny - 1, decomp.npy);
  local_y.first += 1;
  const int storage_y_first = local_y.first - 2;
  const int storage_y_count = local_y.count + 4;
  channel::DnsState state;
  state.resize({decomp.nxB,
                storage_y_count,
                start_field.grid().nz,
                solver_components,
                storage_y_first,
                local_y.first,
                local_y.count});
  state.fill(channel::Complex(0.0, 0.0));
  for (int component = 0; component < base_components; ++component) {
    state.copy_component_from_host(component,
                                   extract_storage_x_window(start_field,
                                                            storage_y_first,
                                                            storage_y_count,
                                                            decomp.nx0,
                                                            decomp.nxB,
                                                            component));
  }

  std::vector<channel::Complex> ialfa;
  std::vector<channel::Complex> ibeta;
  std::vector<double> k2;
  make_line_wavenumbers(input, state.grid(), decomp.nx0, ialfa, ibeta, k2);
  const auto derivatives = compact_derivatives(input, local_y);

  channel::DnsNonlinearVelocityRhsConfig nonlinear_cfg;
  nonlinear_cfg.u_component = 0;
  nonlinear_cfg.v_component = 1;
  nonlinear_cfg.w_component = 2;
  nonlinear_cfg.eta_rhs_component = 0;
  nonlinear_cfg.d2v_rhs_component = 1;
  nonlinear_cfg.product_components = product_components;
  nonlinear_cfg.viscosity = input.velocity.viscosity;
  nonlinear_cfg.mean_pressure = channel::Complex(input.velocity.meanpx, input.velocity.meanpz);

  const auto weights = channel::channel_rk3_weights();
  const int steps = std::max(1, input.timestepping.nstep);
  const bool adaptive_dt = input.timestepping.dt <= 0.0 && input.timestepping.cflmax > 0.0;
  double run_dt = choose_regression_dt(input, state, local_y, decomp);
  double next_run_dt = run_dt;
  channel::DnsRungeKuttaTimestepperConfig cfg;
  cfg.dt = run_dt;
  cfg.time = start_header.time;
  cfg.enable_nonlinear_product_transform = true;
  cfg.nonlinear_product_transform.u_component = 0;
  cfg.nonlinear_product_transform.v_component = 1;
  cfg.nonlinear_product_transform.w_component = 2;
  cfg.nonlinear_product_transform.product_components = nonlinear_cfg.product_components;
  const int nxd = fft_fit_at_least(3 * (input.mesh.nx + 1) / 2);
  const int nzd = fft_fit_at_least(3 * input.mesh.nz);
  cfg.nonlinear_product_transform.product_factor =
      1.0 / (2.0 * static_cast<double>(nxd) * static_cast<double>(nzd));
  cfg.nonlinear_product_transform.dealiased_physical_x = 2 * nxd;
  cfg.nonlinear_product_transform.dealiased_physical_z = nzd;
  if (decomp.npxz > 1) {
    cfg.nonlinear_product_transform.distributed_dealiased_fft = true;
    cfg.nonlinear_product_transform.distributed_spectral_x_total = start_field.grid().nx;
    cfg.nonlinear_product_transform.distributed_spectral_x_first = decomp.nx0;
    cfg.nonlinear_product_transform.distributed_npxz = decomp.npxz;
    cfg.nonlinear_product_transform.distributed_ipxz = decomp.ipxz;
    cfg.nonlinear_product_transform.distributed_comm_x = decomp.comm_x;
  }
  cfg.nonlinear_product_transform.enable_velocity_inverse_ffts = true;
  cfg.nonlinear_product_transform.enable_product_forward_ffts = true;
  cfg.nonlinear_product_transform.velocity_inverse_normalization = channel::FftNormalization::InverseLength;
  cfg.nonlinear_product_transform.product_forward_normalization = channel::FftNormalization::InverseLength;
  if (adaptive_dt) {
    cfg.after_product_transform_observer = [&](int stage, const channel::DnsState& observed_state) {
      if (stage == 2) {
        next_run_dt = choose_regression_dt(input, observed_state, local_y, decomp);
      }
    };
  }
  cfg.enable_nonlinear_velocity_rhs = true;
  cfg.nonlinear_velocity_rhs = nonlinear_cfg;
  cfg.enable_velocity_recovery = true;
  cfg.velocity_recovery.u_component = 0;
  cfg.velocity_recovery.v_component = 1;
  cfg.velocity_recovery.w_component = 2;
  cfg.velocity_recovery.eta_component = 0;
  cfg.velocity_recovery.dvdy_component = 2;
  cfg.velocity_recovery.enable_compact_dvdy = true;
  cfg.velocity_recovery.npy = decomp.npy;
  cfg.velocity_recovery.ipy = decomp.ipy;
  cfg.velocity_recovery.exchange_mode = channel::ExchangeMode::Auto;
  cfg.velocity_recovery.comm_y = decomp.comm_y;
  const auto derivative_boundary = compact_derivative_boundary_data(input, local_y, static_cast<int>(state.line_count()));
  cfg.velocity_recovery.lower_boundary_rhs_coeff = derivative_boundary.lower_boundary_rhs_coeff;
  cfg.velocity_recovery.lower_ghost_rhs_coeff = derivative_boundary.lower_ghost_rhs_coeff;
  cfg.velocity_recovery.upper_boundary_rhs_coeff = derivative_boundary.upper_boundary_rhs_coeff;
  cfg.velocity_recovery.upper_ghost_rhs_coeff = derivative_boundary.upper_ghost_rhs_coeff;
  cfg.stages.resize(weights.size());

  for (std::size_t stage = 0; stage < weights.size(); ++stage) {
    cfg.stages[stage].weights = weights[stage];
    cfg.stages[stage].linear.enable_implicit_yline = true;
    cfg.stages[stage].linear.implicit_yline.components = {0, 1};
    cfg.stages[stage].linear.implicit_yline.npy = decomp.npy;
    cfg.stages[stage].linear.implicit_yline.ipy = decomp.ipy;
    cfg.stages[stage].linear.implicit_yline.exchange_mode = channel::ExchangeMode::Auto;
    cfg.stages[stage].linear.implicit_yline.comm_y = decomp.comm_y;
    cfg.stages[stage].linear.implicit_yline.fill_interface_ghosts = true;
    cfg.stages[stage].enable_velocity_mean_correction = decomp.has_average;
    cfg.stages[stage].velocity_mean_correction.eta_component = 0;
    cfg.stages[stage].velocity_mean_correction.w_component = 2;
    cfg.stages[stage].velocity_mean_correction.line = 0;
    cfg.stages[stage].velocity_mean_correction.root = 0;
    cfg.stages[stage].velocity_mean_correction.comm_y = decomp.comm_y;
    cfg.stages[stage].velocity_mean_correction.meanflowx = input.velocity.meanflowx;
    cfg.stages[stage].velocity_mean_correction.meanflowz = input.velocity.meanflowz;
    cfg.stages[stage].velocity_mean_correction.y = y_coordinates(input);
    const auto boundaries = velocity_boundary_data(input, local_y, static_cast<int>(state.line_count()));
    cfg.stages[stage].velocity_mean_correction.boundaries = compact_boundaries_from_yline(boundaries.eta);
    cfg.stages[stage].velocity_mean_correction.label =
        "regression_velocity_mean_correction_stage_" + std::to_string(stage);
  }

  channel::DnsRungeKuttaTimestepper stepper;
  stepper.prepare(state, std::move(cfg));
  stepper.copy_nonlinear_line_wavenumbers_from_host(ialfa, ibeta, k2);
  stepper.copy_nonlinear_y_derivatives_from_host(derivatives);
  stepper.copy_recovery_line_wavenumbers_from_host(ialfa, ibeta, k2);
  stepper.copy_recovery_y_derivatives_from_host(derivatives);
  const auto compact_dvdy = assemble_compact_derivative_operator(state, derivatives);
  stepper.copy_recovery_compact_dvdy_coefficients_from_host(compact_dvdy.ds,
                                                            compact_dvdy.dl,
                                                            compact_dvdy.d,
                                                            compact_dvdy.du,
                                                            compact_dvdy.dw);
  stepper.copy_recovery_compact_dvdy_boundary_data_from_host(derivative_boundary.boundary);
  add_velocity_yline(stepper,
                     input,
                     state,
                     local_y,
                     input.mesh.ny + 1,
                     k2,
                     derivatives,
                     weights,
                     run_dt,
                     input.velocity.viscosity);
  if (decomp.has_average) {
    add_velocity_mean_correction_matrix(stepper, input, weights, run_dt, input.velocity.viscosity);
  }

  for (int step = 0; step < steps; ++step) {
    next_run_dt = run_dt;
    stepper.advance_one_step(state);
    const int local_bad_after_step = count_nonfinite_components(state, base_components);
    int bad_after_step = 0;
    MPI_Allreduce(&local_bad_after_step, &bad_after_step, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (bad_after_step != 0) {
      if (runtime.rank() == 0) {
        std::cerr << "Regression test FAILED after step " << (step + 1)
                  << ": nonfinite_values=" << bad_after_step
                  << " dt=" << run_dt << '\n';
      }
      return 1;
    }
    if (step + 1 < steps) {
      if (adaptive_dt) {
        run_dt = next_run_dt;
        stepper.set_dt(run_dt);
      }
      add_velocity_yline(stepper,
                         input,
                         state,
                         local_y,
                         input.mesh.ny + 1,
                         k2,
                         derivatives,
                         weights,
                         run_dt,
                         input.velocity.viscosity);
      if (decomp.has_average) {
        add_velocity_mean_correction_matrix(stepper, input, weights, run_dt, input.velocity.viscosity);
      }
    }
  }

  double local_error = 0.0;
  double local_max_diff = 0.0;
  int local_bad_values = 0;
  int local_max_component = -1;
  int local_max_y = -1;
  int local_max_line = -1;
  channel::Complex local_max_got(0.0, 0.0);
  channel::Complex local_max_want(0.0, 0.0);
  for (int component = 0; component < base_components; ++component) {
    const auto got = extract_slice(state, local_y, component);
    const auto want = extract_slice_x_window(expected_field, local_y, decomp.nx0, decomp.nxB, component);
    for (std::size_t i = 0; i < got.size(); ++i) {
      const auto diff = got[i] - want[i];
      const double diff_mag = magnitude(diff);
      const double got_mag = magnitude(got[i]);
      const double want_mag = magnitude(want[i]);
      if (!std::isfinite(diff_mag) || !std::isfinite(got_mag) || !std::isfinite(want_mag)) {
        ++local_bad_values;
        continue;
      }
      local_error += diff_mag * diff_mag;
      if (diff_mag > local_max_diff) {
        local_max_diff = diff_mag;
        local_max_component = component;
        local_max_y = local_y.first + static_cast<int>(i / state.line_count());
        local_max_line = static_cast<int>(i % state.line_count());
        local_max_got = got[i];
        local_max_want = want[i];
      }
    }
  }

  const double diffnorm = std::sqrt(sum_allreduce(local_error));
  const double max_diff = max_allreduce(local_max_diff);
  int bad_values = 0;
  MPI_Allreduce(&local_bad_values, &bad_values, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (bad_values != 0 || !std::isfinite(diffnorm) || diffnorm > 1.0e-12) {
    if (runtime.rank() == 0) {
      std::cerr << "Regression test FAILED, diffnorm=" << diffnorm
                << " max_diff=" << max_diff
                << " nonfinite_values=" << bad_values
                << " steps=" << steps
                << " dt=" << run_dt << '\n';
    }
    if (local_max_diff == max_diff) {
      const int z = local_max_line >= 0 ? local_max_line / state.grid().nx : -1;
      const int x = local_max_line >= 0 ? local_max_line - z * state.grid().nx : -1;
      std::cerr << std::scientific << std::setprecision(16)
                << "Regression local max rank=" << runtime.rank()
                << " component=" << local_max_component
                << " y=" << local_max_y
                << " z=" << z
                << " x=" << x
                << " got=(" << local_max_got.real() << "," << local_max_got.imag() << ")"
                << " want=(" << local_max_want.real() << "," << local_max_want.imag() << ")"
                << " diff=" << local_max_diff << '\n';
      if (local_max_component >= 0 && local_max_line >= 0) {
        const auto got = extract_slice(state, local_y, local_max_component);
        const auto want = extract_slice_x_window(expected_field, local_y, decomp.nx0, decomp.nxB, local_max_component);
        for (int y = 0; y < local_y.count; ++y) {
          const auto p = static_cast<std::size_t>(y) * state.line_count() +
                         static_cast<std::size_t>(local_max_line);
          const auto diff = got[p] - want[p];
          std::cerr << std::scientific << std::setprecision(16)
                    << "Regression worst line rank=" << runtime.rank()
                    << " component=" << local_max_component
                    << " y=" << (local_y.first + y)
                    << " z=" << z
                    << " x=" << x
                    << " got=(" << got[p].real() << "," << got[p].imag() << ")"
                    << " want=(" << want[p].real() << "," << want[p].imag() << ")"
                    << " diff=(" << diff.real() << "," << diff.imag() << ")"
                    << " diff_abs=" << magnitude(diff) << '\n';
        }
      }
    }
    return 1;
  }

  if (runtime.rank() == 0) {
    std::cout << "Regression test PASSED, diffnorm=" << diffnorm
              << " max_diff=" << max_diff
              << " steps=" << steps
              << " dt=" << run_dt << '\n';
  }
  return 0;
}

} // namespace

int main(int argc, char** argv) {
  channel::Runtime runtime(argc, argv);
  if (argc != 5 || std::string{argv[1]} != "--regression") {
    if (runtime.rank() == 0) {
      std::cerr << "usage: test_regression_cpp --regression <dns.in> <start_field> <end_field>\n";
    }
    return 1;
  }
  try {
    return run_velocity_regression(runtime, argv[2], argv[3], argv[4]);
  } catch (const std::exception& ex) {
    if (runtime.rank() == 0) {
      std::cerr << "Regression test FAILED: " << ex.what() << '\n';
    }
    return 1;
  }
}
