#pragma once

#include "channel/device_vector.hpp"
#include "channel/dnsdata.hpp"
#include "channel/mpi.hpp"

#include <string>

namespace channel {

struct MpiTransposeConfig {
  int values_per_peer = 0;
  ExchangeMode exchange_mode = ExchangeMode::AllToAll;
  MpiComm comm = world_comm();
  std::string label = "mpi_transpose";
};

class MpiTransposePlan {
public:
  MpiTransposePlan() = default;
  explicit MpiTransposePlan(MpiTransposeConfig cfg);

  void configure(MpiTransposeConfig cfg);
  void copy_send_from_host(const std::vector<Complex>& values);
  void execute();

  [[nodiscard]] int peer_count() const { return peer_count_; }
  [[nodiscard]] int values_per_peer() const { return cfg_.values_per_peer; }
  [[nodiscard]] std::size_t buffer_size() const;
  [[nodiscard]] DeviceVector<Complex>& send_buffer() { return send_; }
  [[nodiscard]] DeviceVector<Complex>& recv_buffer() { return recv_; }
  [[nodiscard]] const DeviceVector<Complex>& recv_buffer() const { return recv_; }
  [[nodiscard]] std::vector<Complex> recv_host() const { return recv_.copy_to_host(); }

private:
  MpiTransposeConfig cfg_;
  int peer_count_ = 1;
  DeviceVector<Complex> send_;
  DeviceVector<Complex> recv_;
};

struct DnsComponentTransposeConfig {
  int component = 0;
  int values_per_peer = 0;
  ExchangeMode exchange_mode = ExchangeMode::AllToAll;
  MpiComm comm = world_comm();
  std::string label = "dns_component_transpose";
};

class DnsComponentTransposePlan {
public:
  void configure(const DnsState& state, DnsComponentTransposeConfig cfg);
  void execute(DnsState& state);

  [[nodiscard]] const DnsGrid& grid() const { return grid_; }
  [[nodiscard]] const DnsComponentTransposeConfig& config() const { return cfg_; }
  [[nodiscard]] MpiTransposePlan& transpose() { return transpose_; }
  [[nodiscard]] const MpiTransposePlan& transpose() const { return transpose_; }

private:
  DnsGrid grid_;
  DnsComponentTransposeConfig cfg_;
  MpiTransposePlan transpose_;
};

} // namespace channel
