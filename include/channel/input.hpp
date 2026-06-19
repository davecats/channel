#pragma once

#include <string>
#include <vector>

namespace channel {

struct MeshInput {
  int nx = 0;
  int ny = 0;
  int nz = 0;
  double alfa0 = 0.0;
  double beta0 = 0.0;
  double stretching = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
};

struct VelocityInput {
  double reynolds = 0.0;
  double viscosity = 0.0;
  double meanpx = 0.0;
  double meanpz = 0.0;
  double meanflowx = 0.0;
  double meanflowz = 0.0;
  double u0 = 0.0;
  double uN = 0.0;
};

struct ScalarInput {
  int nphi = 0;
  std::vector<double> inverse_prandtl;
};

struct TimesteppingInput {
  double dt = 0.0;
  double cflmax = 0.0;
  double time = 0.0;
  double t_max = 0.0;
  int nstep = 0;
};

struct ParallelInput {
  int npy = 1;
  bool npy_was_set = false;
};

struct ChannelInput {
  MeshInput mesh;
  VelocityInput velocity;
  ScalarInput scalars;
  TimesteppingInput timestepping;
  ParallelInput parallel;
};

ChannelInput read_channel_input(const std::string& path);

} // namespace channel
