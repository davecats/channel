#pragma once

#include "channel/config.hpp"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace channel {

using Complex = Kokkos::complex<double>;
static_assert(sizeof(Complex) == 2 * sizeof(double),
              "channel::Complex must be ABI-compatible with MPI_C_DOUBLE_COMPLEX payload size");
static_assert(alignof(Complex) >= alignof(double),
              "channel::Complex alignment must be compatible with double payloads");

enum class ExchangeMode : std::int32_t {
  Auto = 0,
  AllToAll = 1,
  AllGather = 2,
};

struct LineRange {
  int first = 0;
  int count = 0;
};

inline LineRange split_range(int rank, int nitems, int nranks) {
  const int base = nitems / nranks;
  const int remainder = nitems % nranks;
  const int count = base + (rank < remainder ? 1 : 0);
  const int first = rank * base + (rank < remainder ? rank : remainder);
  return {first, count};
}

inline std::vector<int> default_schur_pass_counts(int npy) {
  std::vector<int> passes;
  if (npy < 1) {
    throw std::runtime_error("default_schur_pass_counts requires npy >= 1");
  }
  int remaining = npy;
  while (remaining > 1) {
    int factor = remaining;
    if (remaining % 4 == 0) {
      factor = 4;
    } else if (remaining % 3 == 0) {
      factor = 3;
    } else if (remaining % 2 == 0) {
      factor = 2;
    }
    passes.push_back(factor);
    remaining /= factor;
  }
  return passes;
}

inline bool valid_schur_pass_sequence(int npy, const std::vector<int>& passes) {
  if (npy == 1 && passes.empty()) return true;
  int product = 1;
  for (const int pass : passes) {
    if (!(pass == 2 || pass == 3 || pass == 4 || pass == 6 || pass == 8)) return false;
    product *= pass;
  }
  return product == npy;
}

} // namespace channel
