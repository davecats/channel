#include "channel/decomposition.hpp"

#include <algorithm>
#include <stdexcept>

namespace channel {

void Decomposition::initialize(const DecompositionConfig& cfg, MpiComm world) {
  finalize();

  nproc = mpi_size(world);
  iproc = mpi_rank(world);
  npy = cfg.npy;
  if (npy < 1) throw std::runtime_error("Decomposition requires npy >= 1");
  if (nproc % npy != 0) throw std::runtime_error("MPI rank count must be divisible by npy");

  npxz = nproc / npy;
  ipy = iproc / npxz;
  ipxz = iproc % npxz;

#if defined(CHANNEL_HAS_MPI)
  comm_x = mpi_split(world, ipy, ipxz);
  comm_y = mpi_split(world, ipxz, ipy);
#else
  comm_x = 0;
  comm_y = 0;
#endif

  ny0 = 1 + ipy * (cfg.ny - 1) / npy;
  nyN = (ipy + 1) * (cfg.ny - 1) / npy;

  nx0 = ipxz * cfg.nxpp / npxz;
  nxN = (ipxz + 1) * cfg.nxpp / npxz - 1;
  nxB = nxN - nx0 + 1;

  nz0 = ipxz * cfg.nzd / npxz;
  nzN = (ipxz + 1) * cfg.nzd / npxz - 1;
  nzB = nzN - nz0 + 1;

  if (cfg.nxpp % npxz != 0 || cfg.nzd % npxz != 0) {
    throw std::runtime_error("FFT transpose requires npxz to divide nxpp and nzd");
  }

  has_average = (nx0 == 0);
  fft_transpose_is_local = (nzB == cfg.nzd);
  sendcount = nxB * nzB * (nyN - ny0 + 5);

  const int fields = cfg.overlapping ? 2 : 1;
  const std::size_t elems = static_cast<std::size_t>(std::max(1, npxz * sendcount * fields));
  sendbuf.resize("transpose_sendbuf", elems);
  recvbuf.resize("transpose_recvbuf", elems);
}

void Decomposition::finalize() {
  if (sendbuf.empty() && recvbuf.empty()) return;
  sendbuf = DeviceVector<Complex>{};
  recvbuf = DeviceVector<Complex>{};
  mpi_free(comm_x);
  mpi_free(comm_y);
}

} // namespace channel
