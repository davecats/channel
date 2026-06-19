#pragma once

#include "channel/dnsdata.hpp"
#include "channel/input.hpp"
#include "channel/mean_correction.hpp"
#include "channel/nonlinear.hpp"
#include "channel/types.hpp"
#include "channel/yline.hpp"

#include <array>
#include <vector>

namespace channel {

struct BandedCoefficients {
  std::vector<Complex> ds;
  std::vector<Complex> dl;
  std::vector<Complex> d;
  std::vector<Complex> du;
  std::vector<Complex> dw;
};

struct VelocityBoundaryData {
  DnsYLineBoundaryData eta;
  DnsYLineBoundaryData v;
  std::vector<Complex> zero;
};

struct CompactDerivativeBoundaryData {
  DnsYLineBoundaryData boundary;
  std::vector<Complex> zero;
  std::array<double, 5> lower_boundary_rhs_coeff{};
  std::array<double, 5> lower_ghost_rhs_coeff{};
  std::array<double, 5> upper_boundary_rhs_coeff{};
  std::array<double, 5> upper_ghost_rhs_coeff{};
};

[[nodiscard]] std::vector<double> y_coordinates(const ChannelInput& input);
[[nodiscard]] int compact_derivative_index(int y, int order, int offset);
[[nodiscard]] double compact_derivative_value(const std::vector<double>& derivatives, int y, int order, int offset);
[[nodiscard]] std::vector<double> compact_derivatives(const ChannelInput& input, const LineRange& local_y);
[[nodiscard]] VelocityBoundaryData velocity_boundary_data(const ChannelInput& input,
                                                          const LineRange& local_y,
                                                          int line_count);
[[nodiscard]] CompactDerivativeBoundaryData compact_derivative_boundary_data(const ChannelInput& input,
                                                                             const LineRange& local_y,
                                                                             int line_count);
void refresh_zero_rhs(VelocityBoundaryData& data);
void refresh_zero_rhs(CompactDerivativeBoundaryData& data);
[[nodiscard]] BandedCoefficients assemble_velocity_operator(const DnsState& state,
                                                            const LineRange& local_y,
                                                            int global_ny,
                                                            const std::vector<double>& k2,
                                                            const std::vector<double>& derivatives,
                                                            double lambda,
                                                            double viscosity,
                                                            bool biharmonic);
[[nodiscard]] BandedCoefficients assemble_compact_derivative_operator(const DnsState& state,
                                                                      const std::vector<double>& derivatives);
[[nodiscard]] CompactBoundaryRows compact_boundaries_from_yline(const DnsYLineBoundaryData& boundary);
[[nodiscard]] std::vector<double> assemble_velocity_mean_correction_matrix(const ChannelInput& input,
                                                                           double lambda,
                                                                           double viscosity);

} // namespace channel
