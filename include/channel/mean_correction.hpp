#pragma once

#include "channel/device_vector.hpp"
#include "channel/dnsdata.hpp"
#include "channel/mpi.hpp"
#include "channel/types.hpp"

#include <array>
#include <span>
#include <string>
#include <vector>

namespace channel {

struct CompactBoundaryRows {
  std::array<double, 5> lower{};
  std::array<double, 5> lower_ghost{};
  std::array<double, 5> upper{};
  std::array<double, 5> upper_ghost{};
  Complex rhs_lower = Complex(0.0, 0.0);
  Complex rhs_lower_ghost = Complex(0.0, 0.0);
  Complex rhs_upper = Complex(0.0, 0.0);
  Complex rhs_upper_ghost = Complex(0.0, 0.0);
};

class FullLineCompactSolver {
public:
  void resize(int ny);
  void copy_system_from_host(std::span<const double> matrix_rows, std::span<const Complex> rhs_line);
  void solve_on_root_and_broadcast(const CompactBoundaryRows& boundaries,
                                   int root,
                                   MpiComm comm = world_comm(),
                                   const std::string& label = "mean_correction_full_line");

  [[nodiscard]] int ny() const { return ny_; }
  [[nodiscard]] DeviceVector<Complex>& line() { return line_; }
  [[nodiscard]] const DeviceVector<Complex>& line() const { return line_; }
  [[nodiscard]] std::vector<Complex> line_host() const { return line_.copy_to_host(); }

private:
  int ny_ = 0;
  DeviceVector<double> matrix_rows_;
  DeviceVector<double> reduced_matrix_;
  DeviceVector<Complex> line_;
  DeviceVector<Complex> rhs_;
};

struct DnsMeanCorrectionConfig {
  int component = 0;
  int line = 0;
  int root = 0;
  MpiComm comm_y = world_comm();
  CompactBoundaryRows boundaries;
  std::string label = "dns_mean_correction";
};

struct DnsVelocityMeanCorrectionConfig {
  int eta_component = 0;
  int w_component = 2;
  int line = 0;
  int root = 0;
  MpiComm comm_y = world_comm();
  CompactBoundaryRows boundaries;
  double meanflowx = 0.0;
  double meanflowz = 0.0;
  std::vector<double> y;
  std::string label = "dns_velocity_mean_correction";
};

class DnsMeanCorrectionStep {
public:
  void prepare(const DnsState& state, DnsMeanCorrectionConfig cfg);
  void copy_matrix_from_host(std::span<const double> matrix_rows);
  void apply(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsMeanCorrectionConfig& config() const { return cfg_; }
  [[nodiscard]] int global_ny() const { return global_ny_; }
  [[nodiscard]] int first_global_y() const { return first_global_y_; }
  [[nodiscard]] FullLineCompactSolver& full_line_solver() { return full_line_solver_; }
  [[nodiscard]] const FullLineCompactSolver& full_line_solver() const { return full_line_solver_; }

private:
  void check_grid(const DnsState& state, const char* caller) const;

  DnsGrid grid_;
  DnsMeanCorrectionConfig cfg_;
  int comm_rank_ = 0;
  int comm_size_ = 1;
  int global_ny_ = 0;
  int first_global_y_ = 1;
  std::vector<int> local_counts_;
  std::vector<int> line_displs_;
  DeviceVector<double> matrix_rows_;
  DeviceVector<Complex> local_line_;
  FullLineCompactSolver full_line_solver_;
};

class DnsVelocityMeanCorrectionStep {
public:
  void prepare(const DnsState& state, DnsVelocityMeanCorrectionConfig cfg);
  void copy_matrix_from_host(std::span<const double> matrix_rows);
  void apply(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsVelocityMeanCorrectionConfig& config() const { return cfg_; }

private:
  [[nodiscard]] double integrate_real_line(const std::vector<Complex>& line) const;
  void check_grid(const DnsState& state, const char* caller) const;

  DnsGrid grid_;
  DnsVelocityMeanCorrectionConfig cfg_;
  int comm_rank_ = 0;
  int comm_size_ = 1;
  int global_ny_ = 0;
  std::vector<double> matrix_rows_host_;
  FullLineCompactSolver correction_solver_;
  DeviceVector<Complex> local_line_;
  DeviceVector<Complex> full_line_;
  DeviceVector<Complex> corrected_eta_line_;
  DeviceVector<Complex> corrected_w_line_;
};

} // namespace channel
