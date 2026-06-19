#pragma once

#include "channel/config.hpp"

namespace channel {

class Runtime {
public:
  Runtime(int& argc, char**& argv);
  Runtime(const Runtime&) = delete;
  Runtime& operator=(const Runtime&) = delete;
  ~Runtime();

  [[nodiscard]] int rank() const { return rank_; }
  [[nodiscard]] int size() const { return size_; }

private:
  bool owns_mpi_ = false;
  bool owns_kokkos_ = false;
  int rank_ = 0;
  int size_ = 1;
};

void fence(const char* label = nullptr);

} // namespace channel
