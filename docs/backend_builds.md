# Backend Build Notes

The port requires MPI, Kokkos, Kokkos-FFT, Umpire, and UmpireSpace.  The local
CPU build uses the bundled third-party checkouts and OpenMP/Serial Kokkos
backends.

## CPU/OpenMP

```bash
cmake -S . -B build-kokkos-umpire-fft-mpi \
  -DCHANNEL_ENABLE_CUDA=OFF \
  -DCHANNEL_ENABLE_HIP=OFF
cmake --build build-kokkos-umpire-fft-mpi -j 4
ctest --test-dir build-kokkos-umpire-fft-mpi --output-on-failure
```

## NVHPC CPU/OpenMP

Use the local NVHPC module when checking the NVIDIA compiler stack on CPU.
With this module loaded, `mpicxx` resolves to `nvc++` and `nvfortran` is also
available for Fortran-side compatibility checks.

```bash
set +u
source ~/.bashrc
module load toolkits/nvhpc/25.5
set -u

cmake -S . -B build-nvhpc-cpu \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCHANNEL_ENABLE_CUDA=OFF \
  -DCHANNEL_ENABLE_HIP=OFF
cmake --build build-nvhpc-cpu -j 4
ctest --test-dir build-nvhpc-cpu --output-on-failure
```

## NVIDIA CUDA

Use Kokkos' `nvcc_wrapper` as the C++ compiler and point CMake's CUDA language
support at a toolkit new enough for the host compiler. On the local RTX A6000
nodes, `module load toolkits/nvhpc/25.5` puts an older CUDA 11.8 `nvcc` first
on `PATH`; the working build uses the module's CUDA 12.9 toolkit explicitly.

```bash
set +u
source ~/.bashrc
module load toolkits/nvhpc/25.5
export CUDA_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9
export PATH=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9/bin:$PATH
set -u

cmake -S . -B build-cuda-129-ampere86 \
  -DCMAKE_CXX_COMPILER="$PWD/third_party/kokkos/bin/nvcc_wrapper" \
  -DCMAKE_CUDA_COMPILER=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9/bin/nvcc \
  -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++ \
  -DCMAKE_CUDA_ARCHITECTURES=86 \
  -DMPI_CXX_COMPILER=mpicxx \
  -DCUDAToolkit_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9 \
  -DCHANNEL_ENABLE_CUDA=ON \
  -DCHANNEL_ENABLE_HIP=OFF \
  -DKokkos_ARCH_AMPERE86=ON
cmake --build build-cuda-129-ampere86 -j 4
CHANNEL_UMPIRE_RESOURCE=DEVICE ctest --test-dir build-cuda-129-ampere86 --output-on-failure
```

`CHANNEL_ENABLE_CUDA=ON` enables `Kokkos_ENABLE_CUDA`, requires a CUDA-enabled
Kokkos build, enables Umpire CUDA support, and lets Kokkos-FFT select cuFFT.
`CUDA_ROOT` is required here because Kokkos' `nvcc_wrapper` selects `nvcc`
from `CUDA_ROOT` before searching `PATH`.

## AMD HIP

Use `hipcc`/`amdclang++` from a ROCm environment.  Pick the `Kokkos_ARCH_*` flag
for the target GPU.  This host does not currently have ROCm tools or an AMD GPU
visible: `hipcc`, `amdclang++`, `hipconfig`, and `rocminfo` are not on `PATH`,
and `module avail rocm`, `module avail amd`, and `module avail cray` did not
expose matching modules. The recipe below is the required target command, not a
locally verified build.

```bash
cmake -S . -B build-hip \
  -DCMAKE_CXX_COMPILER=hipcc \
  -DCHANNEL_ENABLE_HIP=ON \
  -DKokkos_ARCH_AMD_GFX90A=ON
cmake --build build-hip -j 4
CHANNEL_UMPIRE_RESOURCE=DEVICE ctest --test-dir build-hip --output-on-failure
```

`CHANNEL_ENABLE_HIP=ON` enables `Kokkos_ENABLE_HIP`, requires a HIP-enabled
Kokkos build, enables Umpire HIP support, and lets Kokkos-FFT select hipFFT.

CUDA and HIP are intentionally mutually exclusive in this project configure.
