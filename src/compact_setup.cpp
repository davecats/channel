#include "channel/compact_setup.hpp"

#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>

namespace channel {

namespace {

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

std::array<double, 5> compact_row_d4(const ChannelInput& input, const std::vector<double>& y, int iy) {
  std::array<double, 5> t{};
  t[0] = 24.0;
  return solve_5x5(interpolation_matrix(y, iy), t);
}

int dns_index(int y, int line, int batch) {
  return line + y * batch;
}

} // namespace

namespace {

void attach_zero_rhs(DnsYLineBoundaryData& boundary, std::span<const Complex> zero) {
  boundary.lower_boundary_rhs = zero;
  boundary.lower_ghost_rhs = zero;
  boundary.upper_boundary_rhs = zero;
  boundary.upper_ghost_rhs = zero;
}

} // namespace

std::vector<double> y_coordinates(const ChannelInput& input) {
  std::vector<double> y(static_cast<std::size_t>(input.mesh.ny + 3), 0.0);
  const auto at = [&](int iy) -> std::size_t { return static_cast<std::size_t>(iy + 1); };
  const double tanh_stretching = std::tanh(input.mesh.stretching);
  for (int iy = -1; iy <= input.mesh.ny + 1; ++iy) {
    const double eta = 2.0 * static_cast<double>(iy) / static_cast<double>(input.mesh.ny) - 1.0;
    const double mapped = tanh_stretching == 0.0 ? eta : std::tanh(input.mesh.stretching * eta) / tanh_stretching;
    y[at(iy)] = input.mesh.ymin + 0.5 * (input.mesh.ymax - input.mesh.ymin) * (mapped + 1.0);
  }
  return y;
}

int compact_derivative_index(int y, int order, int offset) {
  return (y * DnsNonlinearVelocityRhsStage::derivative_orders + order) *
             DnsNonlinearVelocityRhsStage::derivative_stencil +
         (offset + 2);
}

double compact_derivative_value(const std::vector<double>& derivatives, int y, int order, int offset) {
  return derivatives[static_cast<std::size_t>(compact_derivative_index(y, order, offset))];
}

std::vector<double> compact_derivatives(const ChannelInput& input, const LineRange& local_y) {
  const auto y = y_coordinates(input);
  const auto y_at = [&](int global_y) -> double {
    return y[static_cast<std::size_t>(global_y + 1)];
  };
  std::vector<double> derivatives(static_cast<std::size_t>(local_y.count) *
                                      DnsNonlinearVelocityRhsStage::derivative_orders *
                                      DnsNonlinearVelocityRhsStage::derivative_stencil,
                                  0.0);
  auto set_derivative = [&](int local_row, int order, const std::array<double, 5>& coeffs) {
    for (int offset = -2; offset <= 2; ++offset) {
      derivatives[static_cast<std::size_t>(compact_derivative_index(local_row, order, offset))] =
          coeffs[static_cast<std::size_t>(offset + 2)];
    }
  };

  for (int local_row = 0; local_row < local_y.count; ++local_row) {
    const int iy = local_y.first + local_row;
    if (iy < 1 || iy > input.mesh.ny - 1) continue;

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
  return derivatives;
}

VelocityBoundaryData velocity_boundary_data(const ChannelInput& input, const LineRange& local_y, int line_count) {
  VelocityBoundaryData data;
  data.zero.assign(static_cast<std::size_t>(line_count), Complex(0.0, 0.0));
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

  refresh_zero_rhs(data);
  return data;
}

CompactDerivativeBoundaryData compact_derivative_boundary_data(const ChannelInput& input,
                                                               const LineRange& local_y,
                                                               int line_count) {
  CompactDerivativeBoundaryData data;
  data.zero.assign(static_cast<std::size_t>(line_count), Complex(0.0, 0.0));
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

  refresh_zero_rhs(data);
  return data;
}

void refresh_zero_rhs(VelocityBoundaryData& data) {
  const std::span<const Complex> zero_span(data.zero.data(), data.zero.size());
  attach_zero_rhs(data.eta, zero_span);
  attach_zero_rhs(data.v, zero_span);
}

void refresh_zero_rhs(CompactDerivativeBoundaryData& data) {
  const std::span<const Complex> zero_span(data.zero.data(), data.zero.size());
  attach_zero_rhs(data.boundary, zero_span);
}

BandedCoefficients assemble_velocity_operator(const DnsState& state,
                                              const LineRange& local_y,
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
  coeffs.ds.assign(grid.active_values_per_component(), Complex(0.0, 0.0));
  coeffs.dl.assign(grid.active_values_per_component(), Complex(0.0, 0.0));
  coeffs.d.assign(grid.active_values_per_component(), Complex(0.0, 0.0));
  coeffs.du.assign(grid.active_values_per_component(), Complex(0.0, 0.0));
  coeffs.dw.assign(grid.active_values_per_component(), Complex(0.0, 0.0));

  for (int y = 0; y < active_count; ++y) {
    const int global_y = local_y.first + y;
    for (int line = 0; line < lines; ++line) {
      const int p = dns_index(y, line, lines);
      if (global_y == 0 || global_y == global_ny - 1) {
        coeffs.d[static_cast<std::size_t>(p)] = Complex(1.0, 0.0);
        continue;
      }

      std::array<Complex, 5> row{};
      const double k2_line = k2[static_cast<std::size_t>(line)];
      for (int offset = -2; offset <= 2; ++offset) {
        const double d0 = compact_derivative_value(derivatives, y, 0, offset);
        const double d2 = compact_derivative_value(derivatives, y, 2, offset);
        const double d4 = compact_derivative_value(derivatives, y, 3, offset);
        double value = 0.0;
        if (biharmonic) {
          value = lambda * (d2 - k2_line * d0) -
                  viscosity * (d4 - 2.0 * k2_line * d2 + k2_line * k2_line * d0);
        } else {
          value = lambda * d0 - viscosity * (d2 - k2_line * d0);
        }
        row[static_cast<std::size_t>(offset + 2)] = Complex(value, 0.0);
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

BandedCoefficients assemble_compact_derivative_operator(const DnsState& state, const std::vector<double>& derivatives) {
  const auto& grid = state.grid();
  const int lines = static_cast<int>(grid.line_count());
  const int active_count = grid.active_y_count == 0 ? grid.ny : grid.active_y_count;
  BandedCoefficients coeffs;
  coeffs.ds.assign(grid.active_values_per_component(), Complex(0.0, 0.0));
  coeffs.dl.assign(grid.active_values_per_component(), Complex(0.0, 0.0));
  coeffs.d.assign(grid.active_values_per_component(), Complex(0.0, 0.0));
  coeffs.du.assign(grid.active_values_per_component(), Complex(0.0, 0.0));
  coeffs.dw.assign(grid.active_values_per_component(), Complex(0.0, 0.0));

  for (int y = 0; y < active_count; ++y) {
    for (int line = 0; line < lines; ++line) {
      const int p = dns_index(y, line, lines);
      coeffs.ds[static_cast<std::size_t>(p)] =
          Complex(compact_derivative_value(derivatives, y, 0, -2), 0.0);
      coeffs.dl[static_cast<std::size_t>(p)] =
          Complex(compact_derivative_value(derivatives, y, 0, -1), 0.0);
      coeffs.d[static_cast<std::size_t>(p)] =
          Complex(compact_derivative_value(derivatives, y, 0, 0), 0.0);
      coeffs.du[static_cast<std::size_t>(p)] =
          Complex(compact_derivative_value(derivatives, y, 0, 1), 0.0);
      coeffs.dw[static_cast<std::size_t>(p)] =
          Complex(compact_derivative_value(derivatives, y, 0, 2), 0.0);
    }
  }
  return coeffs;
}

CompactBoundaryRows compact_boundaries_from_yline(const DnsYLineBoundaryData& boundary) {
  CompactBoundaryRows rows;
  rows.lower = boundary.lower_boundary_eq;
  rows.lower_ghost = boundary.lower_ghost_eq;
  rows.upper = boundary.upper_boundary_eq;
  rows.upper_ghost = boundary.upper_ghost_eq;
  rows.rhs_lower = Complex(0.0, 0.0);
  rows.rhs_lower_ghost = Complex(0.0, 0.0);
  rows.rhs_upper = Complex(0.0, 0.0);
  rows.rhs_upper_ghost = Complex(0.0, 0.0);
  return rows;
}

std::vector<double> assemble_velocity_mean_correction_matrix(const ChannelInput& input,
                                                             double lambda,
                                                             double viscosity) {
  const LineRange full_y{1, input.mesh.ny - 1};
  const auto derivatives = compact_derivatives(input, full_y);
  std::vector<double> matrix(static_cast<std::size_t>(input.mesh.ny + 1) * 5, 0.0);
  for (int iy = 1; iy <= input.mesh.ny - 1; ++iy) {
    const int local_y = iy - 1;
    for (int offset = -2; offset <= 2; ++offset) {
      const double d0 = compact_derivative_value(derivatives, local_y, 0, offset);
      const double d2 = compact_derivative_value(derivatives, local_y, 2, offset);
      matrix[static_cast<std::size_t>(iy * 5 + offset + 2)] = lambda * d0 - viscosity * d2;
    }
  }
  return matrix;
}

} // namespace channel
