#pragma once

#include "channel/device_vector.hpp"
#include "channel/mpi.hpp"
#include "channel/types.hpp"

#include <vector>

namespace channel {

struct SchurLevelConfig {
  int level_index = 0;
  int arity = 1;
  int child_id = 0;
  int child_span = 1;
  int prev_first = 0;
  int prev_count = 0;
  LineRange owned;
  ExchangeMode exchange_mode = ExchangeMode::Auto;
  MpiComm comm = world_comm();
};

class SchurLevel {
public:
  static constexpr int row_width = 20;
  static constexpr int value_width = 8;
  static constexpr int max_arity = 8;
  static constexpr int max_rows = 4 * max_arity;
  static constexpr int bandwidth = 5;

  SchurLevel() = default;
  explicit SchurLevel(SchurLevelConfig cfg);
  SchurLevel(const SchurLevel&) = delete;
  SchurLevel& operator=(const SchurLevel&) = delete;
  SchurLevel(SchurLevel&& other) noexcept;
  SchurLevel& operator=(SchurLevel&& other) noexcept;
  ~SchurLevel();

  [[nodiscard]] const SchurLevelConfig& config() const { return cfg_; }
  [[nodiscard]] int arity() const { return cfg_.arity; }
  [[nodiscard]] int owned_count() const { return cfg_.owned.count; }
  [[nodiscard]] int prev_count() const { return cfg_.prev_count; }
  [[nodiscard]] int level_index() const { return cfg_.level_index; }
  [[nodiscard]] int solve_count(bool solve_redundant = false) const;
  [[nodiscard]] int reduced_row_count() const { return 4 * cfg_.arity; }

  [[nodiscard]] DeviceVector<Complex>& rows() { return rows_; }
  [[nodiscard]] const DeviceVector<Complex>& rows() const { return rows_; }
  [[nodiscard]] DeviceVector<Complex>& values() { return values_; }
  [[nodiscard]] const DeviceVector<Complex>& values() const { return values_; }
  [[nodiscard]] DeviceVector<Complex>& recovery_basis() { return recovery_basis_; }
  [[nodiscard]] const DeviceVector<Complex>& recovery_basis() const { return recovery_basis_; }
  [[nodiscard]] DeviceVector<Complex>& row_send_buffer() { return row_send_; }
  [[nodiscard]] DeviceVector<Complex>& row_recv_buffer() { return row_recv_; }
  [[nodiscard]] const DeviceVector<Complex>& row_recv_buffer() const { return row_recv_; }
  [[nodiscard]] DeviceVector<Complex>& value_send_buffer() { return value_send_; }
  [[nodiscard]] DeviceVector<Complex>& value_recv_buffer() { return value_recv_; }
  [[nodiscard]] const DeviceVector<Complex>& value_recv_buffer() const { return value_recv_; }

  void resize_buffers();
  void release_comm();
  void pack_rows_for_exchange(const DeviceVector<Complex>& source_rows);
  void exchange_rows();
  void compose_from_exchanged_rows(const DeviceVector<Complex>& exchanged_rows, bool solve_redundant = false);
  void compose_from_exchanged_rows(bool solve_redundant = false);
  void seed_root_values(bool solve_redundant = false);
  void pack_recovered_values(DeviceVector<Complex>& packed_values) const;
  void pack_recovered_values();
  void exchange_values();
  void unpack_exchanged_values(DeviceVector<Complex>& child_values) const;
  void recover_redundant_root_values(DeviceVector<Complex>& child_values) const;

  [[nodiscard]] std::vector<Complex> rows_host() const { return rows_.copy_to_host(); }
  [[nodiscard]] std::vector<Complex> values_host() const { return values_.copy_to_host(); }
  [[nodiscard]] std::vector<Complex> recovery_basis_host() const { return recovery_basis_.copy_to_host(); }

private:
  SchurLevelConfig cfg_;
  DeviceVector<Complex> rows_;
  DeviceVector<Complex> values_;
  DeviceVector<Complex> recovery_basis_;
  DeviceVector<Complex> row_send_;
  DeviceVector<Complex> row_recv_;
  DeviceVector<Complex> value_send_;
  DeviceVector<Complex> value_recv_;
  std::vector<LineRange> child_ranges_;
  std::vector<int> row_send_counts_;
  std::vector<int> row_send_displs_;
  std::vector<int> row_recv_counts_;
  std::vector<int> row_recv_displs_;
  std::vector<int> value_send_counts_;
  std::vector<int> value_send_displs_;
  std::vector<int> value_recv_counts_;
  std::vector<int> value_recv_displs_;
};

struct SchurSolverConfig {
  int npy = 1;
  int ipy = 0;
  int nlines = 0;
  std::vector<int> pass_counts;
  ExchangeMode exchange_mode = ExchangeMode::Auto;
  MpiComm comm_y = world_comm();
};

class SchurSolver {
public:
  void prepare(const SchurSolverConfig& cfg);
  void release();
  void solve_from_leaf_rows(const DeviceVector<Complex>& leaf_rows, DeviceVector<Complex>& leaf_values);

  [[nodiscard]] const std::vector<SchurLevel>& levels() const { return levels_; }
  [[nodiscard]] std::vector<SchurLevel>& levels() { return levels_; }

private:
  SchurSolverConfig cfg_;
  std::vector<SchurLevel> levels_;
};

} // namespace channel
