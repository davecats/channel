# channel C++/Kokkos Port

This branch contains an experimental C++20/Kokkos port of the core velocity-solve path from `channel`.

For the original Fortran project documentation, use the main branch:

https://github.com/davecats/channel/tree/main

## Status

This port was produced by Codex with oversight from Jona Neuhauser. It has not been reviewed in detail. Treat it as a working MVP and benchmarking artifact, not as reviewed production code.

What currently works:

- C++/Kokkos build with MPI, Kokkos, Kokkos-FFT, Umpire, and umpire-kokkos.
- CUDA build tested locally with NVHPC 25.5 and CUDA 12.9.
- Native `channel --run dns.in --steps N` path for the velocity solve.
- Runge-Kutta nonlinear velocity pipeline.
- Local dealiased nonlinear FFT path.
- Assembled implicit velocity y-line solve and velocity recovery.
- Mean correction in the velocity path.
- `npy=1` and `npy=2` runs have been exercised during the port.
- Regression harnesses exist for the accepted MVP subset, with pressure/convvelo intentionally excluded.

Known caveats:

- This is an MVP for the velocity solve only. Pressure/projection, convvelo, and production restart/output compatibility are not complete.
- The Kokkos port is slower than the Fortran GPU implementation on the current benchmark. On the local `dns.in`, `2` ranks, `npy=2`, `npxz=1`, `5` timestep CUDA profile, the Kokkos separable-FFT path is about `2.3 s/step`; the Fortran baseline is about `1.6-1.7 s/step`.
- The velocity y-line solve is roughly at parity with Fortran. The main slowdown is the local nonlinear FFT path.
- GPU-aware MPI is not fully resolved. Recent traces show a CUDA peer-to-peer memcpy row, but OpenMPI still prints `cuMemHostRegister ... Registration cache: smcuda`, so the MPI path should not be considered cleanly verified.
- Memory behavior is still not optimized. Kokkos view initialization, memset, D2D copies, and KokkosFFT scratch/transpose behavior are still visible in profiles.
- `CHANNEL_USE_DEALIASED_2D_FFT=1` compiles and runs, and removes KokkosFFT transpose kernels in a one-step trace, but it is slower in wall time on the current benchmark (`~5.0 s/step` for `5` steps). It is experimental and not the default.
- The code has accumulated debug hooks and profiling markers from the porting process. Cleanup is still required before review.

## Submodules

Third-party Kokkos-related dependencies are tracked as submodules:

- `third_party/kokkos`
- `third_party/kokkos-fft`
- `third_party/kokkos-tools`
- `third_party/umpire`
- `third_party/umpire-kokkos`

Initialize them after cloning or checking out this branch:

```bash
git submodule update --init --recursive
```

## Build: CUDA/NVHPC

The instructions assume `module load toolkits/nvhpc/25.5` provides the correct `mpicxx`, `mpirun`, `nvc++`, and `nsys` on `PATH`. CUDA 12.9 is selected explicitly.

```bash
. /etc/profile.d/lmod.sh
module load toolkits/nvhpc/25.5

export CUDA_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9

cmake -S . -B build-cuda-129-ampere86-release-gcc \
  -DCMAKE_CXX_COMPILER="$PWD/third_party/kokkos/bin/nvcc_wrapper" \
  -DCMAKE_CUDA_COMPILER="$CUDA_ROOT/bin/nvcc" \
  -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++ \
  -DCMAKE_CUDA_ARCHITECTURES=86 \
  -DMPI_CXX_COMPILER=mpicxx \
  -DCUDAToolkit_ROOT="$CUDA_ROOT" \
  -DCHANNEL_ENABLE_CUDA=ON \
  -DCHANNEL_ENABLE_HIP=OFF \
  -DKokkos_ARCH_AMPERE86=ON

cmake --build build-cuda-129-ampere86-release-gcc --target channel -j 4
```

## Run

Default separable local dealiased FFT path:

```bash
. /etc/profile.d/lmod.sh
module load toolkits/nvhpc/25.5

export CUDA_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9
export CHANNEL_UMPIRE_RESOURCE=DEVICE

mpirun -np 2 build-cuda-129-ampere86-release-gcc/channel-bin/channel \
  --run dns.in --steps 5
```

Experimental local dealiased 2D KokkosFFT path:

```bash
export CHANNEL_USE_DEALIASED_2D_FFT=1

mpirun -np 2 build-cuda-129-ampere86-release-gcc/channel-bin/channel \
  --run dns.in --steps 5
```

Do not use `CHANNEL_USE_DEALIASED_2D_FFT=1` as the default path yet. It currently runs slower than the separable path.

## Build Profiling Module

The Kokkos profiling module used for Nsight ranges is `kp_nvtx_connector` from the `third_party/kokkos-tools` submodule.

```bash
. /etc/profile.d/lmod.sh
module load toolkits/nvhpc/25.5

export CUDA_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9

cmake -S third_party/kokkos-tools -B build-kokkos-tools-nvtx \
  -DCMAKE_CXX_COMPILER=nvc++ \
  -DCUDAToolkit_ROOT="$CUDA_ROOT" \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkosTools_ENABLE_MPI=OFF \
  -DKokkosTools_ENABLE_TESTS=OFF \
  -DKokkosTools_ENABLE_EXAMPLES=OFF

cmake --build build-kokkos-tools-nvtx --target kp_nvtx_connector -j 4
```

The resulting shared library is:

```bash
build-kokkos-tools-nvtx/profiling/nvtx-connector/libkp_nvtx_connector.so
```

## Profile

Run with `KOKKOS_TOOLS_LIBS` pointing at the profiling module:

```bash
. /etc/profile.d/lmod.sh
module load toolkits/nvhpc/25.5

export CHANNEL_UMPIRE_RESOURCE=DEVICE
export KOKKOS_TOOLS_LIBS="$PWD/build-kokkos-tools-nvtx/profiling/nvtx-connector/libkp_nvtx_connector.so"

mpirun -np 2 nsys profile \
  --trace=cuda,nvtx,mpi \
  --sample=none \
  --cpuctxsw=none \
  --force-overwrite=true \
  --output profiles/cpp_dns_npy2_rank%q{OMPI_COMM_WORLD_RANK} \
  build-cuda-129-ampere86-release-gcc/channel-bin/channel \
  --run dns.in --steps 5
```

Export summaries:

```bash
nsys stats --force-export=true --report nvtx_sum --format csv \
  --output profiles/cpp_dns_npy2_rank1_stats_nvtx_sum \
  profiles/cpp_dns_npy2_rank1.nsys-rep

nsys stats --force-export=true --report cuda_gpu_kern_sum --format csv \
  --output profiles/cpp_dns_npy2_rank1_stats_cuda_gpu_kern_sum \
  profiles/cpp_dns_npy2_rank1.nsys-rep

nsys stats --force-export=true --report cuda_gpu_mem_time_sum --format csv \
  --output profiles/cpp_dns_npy2_rank1_stats_cuda_gpu_mem_time_sum \
  profiles/cpp_dns_npy2_rank1.nsys-rep
```
