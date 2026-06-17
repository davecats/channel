#pragma once

#include "channel/device_vector.hpp"
#include "channel/dnsdata.hpp"
#include "channel/schur.hpp"
#include "channel/types.hpp"

#include <memory>
#include <span>
#include <vector>

namespace channel {

class YLineWorkspace {
public:
  void resize(int n, int batch_count);
  void copy_from_host(const std::vector<Complex>& ds,
                      const std::vector<Complex>& dl,
                      const std::vector<Complex>& d,
                      const std::vector<Complex>& du,
                      const std::vector<Complex>& dw,
                      const std::vector<Complex>& x);
  void solve(const std::string& label = "yline_pentadiagonal");

  [[nodiscard]] int n() const { return n_; }
  [[nodiscard]] int batch_count() const { return batch_count_; }
  [[nodiscard]] std::vector<Complex> x_host() const { return x_.copy_to_host(); }

private:
  int n_ = 0;
  int batch_count_ = 0;
  DeviceVector<Complex> ds_;
  DeviceVector<Complex> dl_;
  DeviceVector<Complex> d_;
  DeviceVector<Complex> du_;
  DeviceVector<Complex> dw_;
  DeviceVector<Complex> x_;
};

void solve_batched_pentadiagonal(DeviceVector<Complex>& ds,
                                  DeviceVector<Complex>& dl,
                                  DeviceVector<Complex>& d,
                                  DeviceVector<Complex>& du,
                                  DeviceVector<Complex>& dw,
                                  DeviceVector<Complex>& x,
                                  int n,
                                  int batch_count,
                                  const std::string& label);

struct EndpointSchurYLineConfig {
  int local_n = 0;
  int batch_count = 0;
  int npy = 1;
  int ipy = 0;
  std::vector<int> pass_counts;
  ExchangeMode exchange_mode = ExchangeMode::Auto;
  MpiComm comm_y = world_comm();
};

class EndpointSchurYLineSolver {
public:
  void prepare(EndpointSchurYLineConfig cfg);
  void solve(DeviceVector<Complex>& ds,
             DeviceVector<Complex>& dl,
             DeviceVector<Complex>& d,
             DeviceVector<Complex>& du,
             DeviceVector<Complex>& dw,
             DeviceVector<Complex>& x);

  [[nodiscard]] const EndpointSchurYLineConfig& config() const { return cfg_; }
  [[nodiscard]] int exposed_count() const;
  [[nodiscard]] int interior_count() const;
  [[nodiscard]] const SchurSolver& schur_solver() const { return schur_; }

private:
  EndpointSchurYLineConfig cfg_;
  SchurSolver schur_;
  DeviceVector<Complex> batch_ds_;
  DeviceVector<Complex> batch_dl_;
  DeviceVector<Complex> batch_d_;
  DeviceVector<Complex> batch_du_;
  DeviceVector<Complex> batch_dw_;
  DeviceVector<Complex> batch_x_;
  DeviceVector<Complex> leaf_rows_;
  DeviceVector<Complex> leaf_values_;
};

struct DnsYLineComponentConfig {
  int component = 0;
  int npy = 1;
  int ipy = 0;
  std::vector<int> pass_counts;
  ExchangeMode exchange_mode = ExchangeMode::Auto;
  MpiComm comm_y = world_comm();
};

class DnsYLineComponentSolver {
public:
  void prepare(const DnsState& state, DnsYLineComponentConfig cfg);
  void copy_coefficients_from_host(std::span<const Complex> ds,
                                   std::span<const Complex> dl,
                                   std::span<const Complex> d,
                                   std::span<const Complex> du,
                                   std::span<const Complex> dw);
  void solve(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsYLineComponentConfig& config() const { return cfg_; }
  [[nodiscard]] EndpointSchurYLineSolver& endpoint_solver() { return endpoint_solver_; }
  [[nodiscard]] const EndpointSchurYLineSolver& endpoint_solver() const { return endpoint_solver_; }

  [[nodiscard]] DeviceVector<Complex>& ds() { return ds_; }
  [[nodiscard]] DeviceVector<Complex>& dl() { return dl_; }
  [[nodiscard]] DeviceVector<Complex>& d() { return d_; }
  [[nodiscard]] DeviceVector<Complex>& du() { return du_; }
  [[nodiscard]] DeviceVector<Complex>& dw() { return dw_; }

private:
  DnsGrid grid_;
  DnsYLineComponentConfig cfg_;
  EndpointSchurYLineSolver endpoint_solver_;
  DeviceVector<Complex> ds_;
  DeviceVector<Complex> dl_;
  DeviceVector<Complex> d_;
  DeviceVector<Complex> du_;
  DeviceVector<Complex> dw_;
  DeviceVector<Complex> rhs_;
};

struct DnsImplicitYLineStepConfig {
  std::vector<int> components;
  int npy = 1;
  int ipy = 0;
  std::vector<int> pass_counts;
  ExchangeMode exchange_mode = ExchangeMode::Auto;
  MpiComm comm_y = world_comm();
};

class DnsImplicitYLineStep {
public:
  void prepare(const DnsState& state, DnsImplicitYLineStepConfig cfg);
  void copy_coefficients_from_host(int component,
                                   std::span<const Complex> ds,
                                   std::span<const Complex> dl,
                                   std::span<const Complex> d,
                                   std::span<const Complex> du,
                                   std::span<const Complex> dw);
  void solve(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsImplicitYLineStepConfig& config() const { return cfg_; }
  [[nodiscard]] int component_count() const { return static_cast<int>(solvers_.size()); }
  [[nodiscard]] DnsYLineComponentSolver& component_solver(int component);
  [[nodiscard]] const DnsYLineComponentSolver& component_solver(int component) const;

private:
  [[nodiscard]] int solver_index(int component) const;

  DnsGrid grid_;
  DnsImplicitYLineStepConfig cfg_;
  std::vector<std::unique_ptr<DnsYLineComponentSolver>> solvers_;
};

} // namespace channel
