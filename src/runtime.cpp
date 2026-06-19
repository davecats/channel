#include "channel/runtime.hpp"

#include "channel/memory.hpp"

#include <mpi.h>

namespace channel {

Runtime::Runtime(int& argc, char**& argv) {
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (!initialized) {
    MPI_Init(&argc, &argv);
    owns_mpi_ = true;
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);

  if (!Kokkos::is_initialized()) {
    Kokkos::initialize(argc, argv);
    owns_kokkos_ = true;
  }
  initialize_memory_pool();
}

Runtime::~Runtime() {
  if (owns_kokkos_ && Kokkos::is_initialized()) {
    Kokkos::finalize();
  }
  if (owns_mpi_) {
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) MPI_Finalize();
  }
}

void fence(const char* label) {
  if (label) {
    Kokkos::fence(label);
  } else {
    Kokkos::fence();
  }
}

} // namespace channel
