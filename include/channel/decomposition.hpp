#pragma once

#include "channel/device_vector.hpp"
#include "channel/mpi.hpp"

namespace channel {

struct DecompositionConfig {
  int nxpp = 1;
  int nz = 0;
  int ny = 1;
  int nzd = 1;
  int nphi = 0;
  bool overlapping = false;
  int npy = 1;
};

struct Decomposition {
  int nproc = 1;
  int iproc = 0;
  int npy = 1;
  int npxz = 1;
  int ipy = 0;
  int ipxz = 0;
  int nx0 = 0;
  int nxN = 0;
  int nxB = 1;
  int nz0 = 0;
  int nzN = 0;
  int nzB = 1;
  int ny0 = 1;
  int nyN = 0;
  int sendcount = 0;
  bool has_average = true;
  bool fft_transpose_is_local = true;
  MpiComm comm_x = world_comm();
  MpiComm comm_y = world_comm();
  DeviceVector<Complex> sendbuf;
  DeviceVector<Complex> recvbuf;

  void initialize(const DecompositionConfig& cfg, MpiComm world = world_comm());
  void finalize();
};

} // namespace channel
