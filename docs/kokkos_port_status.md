# C++/Kokkos Port Status

This worktree is the C++ port branch. The current implementation is a bootstrap,
not a full DNS replacement yet.

Implemented:

- C++20 `channel_core` target.
- MPI is required; the no-MPI fallback path has been removed from the build.
- Vendored Kokkos checkout in `third_party/kokkos` (`601d7ac`) with OpenMP/Serial CPU backends in the local build.
- Vendored Umpire checkout in `third_party/umpire` (`f544027`) plus `third_party/umpire-kokkos`.
- Explicit `CHANNEL_ENABLE_CUDA` / `CHANNEL_ENABLE_HIP` CMake switches drive
  the Kokkos, Umpire, and Kokkos-FFT backend configuration. See
  `docs/backend_builds.md`.
- `DefaultMemorySpace` is `UmpireSpace<Kokkos::DefaultExecutionSpace::memory_space>`.
- All `DeviceVector` storage uses the Umpire-backed default memory space, so DNS, Schur, y-line, MPI buffers, FFT scratch, and workspace allocations go through Umpire.
- Vendored Kokkos-FFT checkout in `third_party/kokkos-fft` (`106cb9d`) linked through `KokkosFFT::fft`; the local CPU build uses FFTW OpenMP/Serial backends.
- `DnsState` with component-separated, y-line/interleaved storage.
- `DnsState` component subviews for device-side per-component kernels without
  whole-state host staging.
- Disabled `convvelo` and `pressure` hooks for this milestone.
- `WorkspaceArena` with scoped leases over Umpire-backed `DeviceVector` storage.
- GPU-aware MPI wrapper functions that keep device pointers as the primary path.
- GPU-aware `MPI_Gatherv` wrapper for root-side device gathers of distributed
  full-line data.
- GPU-aware `MPI_Bcast` wrapper for device-resident full-line broadcasts.
- GPU-aware `MPI_Alltoallv` wrapper for uneven Schur line partitions.
- `Decomposition` with the current `npxz/npy` topology.
- `MpiTransposePlan` with device-resident send/recv buffers and collective-ready layout.
- MPI transpose coverage now exercises both `AllToAll` and `AllGather` with
  rank-routed payloads, including a dedicated two-rank CTest path, and asserts
  that send/recv buffers use the Umpire-backed `DefaultMemorySpace`.
- `LocalFftPlan` with explicit normalization and Kokkos-FFT C2C dispatch.
- `DnsComponentFftPlan` bridge from `DnsState` component storage to Kokkos-FFT,
  using Umpire-backed scratch and device-side pack/unpack of the selected
  component.
- `DnsComponentTransposePlan` bridge from `DnsState` component storage to the
  GPU-aware MPI transpose buffers, covering both `AllToAll` and `AllGather`
  exchange semantics.
- Portable Kokkos batched pentadiagonal y-line solver.
- `SchurSolver` with one `SchurLevel` object per configured level.
- `SchurSolver` validates y-rank topology, non-empty line counts, and pass
  sequence compatibility before constructing level communicators.
- Packed endpoint Schur driver:
  upward row packing/exchange, level composition, root seeding,
  downward value recovery/exchange, and leaf value unpack.
- Multilevel Schur row packing handles active row prefixes from overallocated
  level buffers, including the root `AllGather` path.
- Multilevel Schur value unpack also copies active prefixes from overallocated
  level buffers during downward recovery.
- `EndpointSchurYLineSolver` bridge from local pentadiagonal coefficients/RHS
  to reduced Schur leaf rows and reconstructed local y-line solutions.
- Single-rank `EndpointSchurYLineSolver` path dispatches to the local Kokkos
  batched pentadiagonal solve without Schur exchange.
- `DnsYLineComponentSolver` bridge from `DnsState` component storage to the
  endpoint Schur y-line solve, with Umpire-backed coefficient/RHS buffers and
  device-side pack/unpack of the selected component.
- `DnsImplicitYLineStep` orchestration over multiple DNS components, with one
  component solver per selected component and a single `solve(DnsState&)`
  entry point for the implicit y-line stage.
- `DnsLinearStep` orchestration over DNS component FFT, GPU-aware MPI transpose,
  implicit y-line, and inverse FFT stages from one `advance(DnsState&)` entry
  point.
- Endpoint Schur manufactured-solution coverage includes both uniform local
  y slabs and rank-dependent uneven local y slabs.
- Endpoint Schur manufactured-solution coverage includes small local y slabs
  where interior ranks have exactly one interior row after endpoint exposure.
- Focused Schur tests for isolated composition/recovery and a two-rank uneven
  `Alltoallv` solve path.
- Four-rank two-level Schur test covering first-level `Alltoallv` plus
  automatic redundant root `AllGather`.
- Six-rank two-level Schur test covering a non-binary `{3,2}` pass tree,
  including ternary first-level composition and binary redundant root recovery.
- Six-rank two-level Schur test covering a `{2,3}` pass tree whose final/root
  level uses `Alltoallv`, exercising non-redundant multilevel root recovery.
- Schur levels support zero-owned parent line ranges, which occur in uneven
  multilevel partitions when the active line count is smaller than the y-rank
  count.
- One-rank manufactured-solution endpoint y-line test for the direct local
  solver path.
- Two-rank manufactured-solution endpoint y-line test that compares the
  reconstructed distributed solution against a known global pentadiagonal
  solution in both forced `Alltoallv` and automatic root `AllGather` modes.
- Three-rank manufactured-solution endpoint y-line test covering the default
  ternary Schur root pass.
- Four-rank manufactured-solution endpoint y-line test over a `{2,2}` Schur
  pass tree, covering the full endpoint bridge through multilevel recovery.
- Six-rank manufactured-solution endpoint y-line tests over default `{3,2}` and
  explicit `{2,3}` pass trees, including ranks with zero owned root-level lines.
- DNS component y-line manufactured-solution tests over one, two, four, and
  six ranks verify the `DnsState` component bridge, Schur exchange variants,
  untouched neighboring components, and device-resident CUDA/Umpire execution.
- DNS implicit y-line step manufactured-solution tests solve multiple DNS
  components with distinct operators, preserve inactive components, and cover
  one-, two-, four-, and six-rank Schur exchange paths.
- DNS FFT/transpose tests verify component-level FFT roundtrips and impulse
  spectra, inactive component preservation, one-rank transpose behavior, and
  two-rank GPU-aware MPI payload routing through `DnsState`.
- DNS linear-step tests chain FFT roundtrip, MPI transpose, and implicit y-line
  stages through one DNS-level step on one and two MPI ranks.
- The `channel` executable now has a native deterministic
  `--linear-fixture [alltoall|allgather]` path that drives `DnsLinearStep`
  through FFT, GPU-aware MPI transpose, implicit Schur y-line, and
  mean-correction stages and fails on numerical mismatch.
- `FullLineCompactSolver` ports the mean-correction full-line compact solve:
  y-root performs the boundary elimination, pentadiagonal factor/backsolve, and
  ghost reconstruction in Kokkos over Umpire-backed storage, then broadcasts the
  complete ghosted line through the GPU-aware MPI path.
- `DnsMeanCorrectionStep` packs a selected `DnsState` component/line on device,
  gathers distributed y slabs to y-root, runs the full-line compact solve,
  broadcasts the ghosted result, and scatters the corrected local rows back
  while preserving neighboring lines and inactive components.
- `DnsLinearStep` can now own mean-correction stages between the implicit
  y-line solve and inverse FFT stages.
- Mean-correction tests compare the Kokkos compact solve against an independent
  dense full-system reference and cover one-rank solve, two-rank gather/broadcast,
  `DnsState` scatterback, and `DnsLinearStep` orchestration.
- `DnsNonlinearVelocityRhsStage` ports the velocity RHS algebra from
  `buildrhs_prepare` / `buildrhs`: it owns line wavenumbers, y-derivative
  stencils, eta/D2v old-RHS history, optional velocity-product construction,
  and accumulates the six velocity product contributions in Kokkos.
- `DnsNonlinearProductTransformStage` wires the nonlinear product leg into the
  solver graph: optional inverse velocity component FFTs, optional GPU-aware
  MPI component transposes, device-side construction of `uu`, `vv`, `ww`,
  `uv`, `vw`, `uw`, optional product FFTs, and optional product transposes.
- `DnsVelocityRecoveryStage` ports the nonzero-mode velocity recovery after
  the implicit solve, computing `dv/dy` and reconstructing `u,w` from eta and
  `dv/dy` while preserving the zero mode for the mean-correction path.
- `DnsVelocityStep` can now chain nonlinear velocity RHS assembly before the
  existing linear FFT/transpose/y-line/mean-correction/inverse stages, and can
  run product transforms before RHS assembly plus velocity recovery after the
  linear solve.
- `DnsRungeKuttaTimestepper` owns the channel RK3 substep weights, advances
  substep time with the Fortran `2/RK_rai(1,i)*deltat` rule, keeps nonlinear
  old-RHS history shared across substages, runs the nonlinear product transform
  and RHS stages per substep, runs one linear step configuration per RK
  substage, and optionally recovers velocity for the next substep.
- The `channel --run [dns.in] [--steps N]` executable path parses the native
  input file, initializes a deterministic velocity state, and runs an
  end-to-end RK velocity timestep through product transforms, nonlinear RHS
  assembly, assembled implicit eta/v y-line solves, and velocity recovery.
  `--identity-yline` remains as a solver-bypass diagnostic.
- The parsed runner now uploads per-stage velocity y-line coefficients by
  assembling benchmark finite-difference Helmholtz/biharmonic operators from
  the RK implicit weight, timestep, viscosity, and local line wavenumbers.
- Focused nonlinear/RK tests compare velocity RHS accumulation against an
  independent host reference, check velocity-product construction and product
  FFT wiring, verify `DnsVelocityStep` orchestration, verify velocity recovery,
  and catch RK old-RHS history sharing across substages.
- Executable timestep coverage includes one-rank and two-rank `channel --run`
  smoke tests on `tests/cpp/tiny_dns.in`, so the MPI Schur y-line path is now
  exercised through the parsed timestep driver.

Not implemented yet:

- Production-equivalent nonlinear transform layout: the current C++ product
  transform stage uses the generic component FFT/transpose plans, but still
  needs the exact dealiased x/z physical layout, half/full FFT choices, and
  pack/unpack ordering used by the Fortran `transform_to_physical` and
  `transform_back_and_build_rhs` path.
- Production compact derivative/boundary setup in the parsed runner. The
  current executable timestep uses generated benchmark finite-difference
  derivative stencils and assembled y-line operators, not the exact Fortran
  compact matrices and boundary elimination.
- Mean-correction integration in the parsed RK runner. `DnsMeanCorrectionStep`
  is ported and tested independently, but the `--run` state layout still needs
  the production active-y/ghost convention before zero-mode correction should
  be wired into the executable timestep.
- Pressure/projection physics.
- End-to-end nonlinear timestep regression around convective transforms,
  pressure/projection, boundary conditions, scalar RHS, and output hooks.
- Full native C++ regression fixtures replacing old restart files and covering
  nonlinear timestep/projection behavior.
- AMD HIP acceptance builds/runs. HIP remains unverified here because
  `hipcc`, `amdclang++`, `hipconfig`, and `rocminfo` are not on `PATH`, and
  `module avail rocm`, `module avail amd`, and `module avail cray` did not
  expose matching modules.

Checks run in this worktree:

- `cmake -S . -B build-kokkos-umpire-fft-mpi`
- `cmake --build build-kokkos-umpire-fft-mpi -j 4`
- `ctest --test-dir build-kokkos-umpire-fft-mpi -R 'test_input_cpp|test_nonlinear_product_transform_cpp|test_runge_kutta_cpp' --output-on-failure`
- `ctest --test-dir build-kokkos-umpire-fft-mpi -R 'test_(nonlinear_product_transform|runge_kutta|velocity_recovery)_cpp' --output-on-failure`
- `ctest --test-dir build-kokkos-umpire-fft-mpi -R test_channel_tiny_timestep_cpp --output-on-failure`
- `ctest --test-dir build-kokkos-umpire-fft-mpi -R 'test_input_cpp|test_nonlinear_product_transform_cpp|test_runge_kutta_cpp|test_channel_tiny_timestep_cpp' --output-on-failure`
- `ctest --test-dir build-kokkos-umpire-fft-mpi --output-on-failure`
- `ctest --test-dir build-kokkos-umpire-fft-mpi -R 'test_nonlinear_product_transform_cpp|test_channel_tiny_timestep.*cpp' --output-on-failure`
- `module load toolkits/nvhpc/25.5`
- `cmake -S . -B build-nvhpc-cpu -DCMAKE_CXX_COMPILER=mpicxx -DCHANNEL_ENABLE_CUDA=OFF -DCHANNEL_ENABLE_HIP=OFF`
- `cmake --build build-nvhpc-cpu -j 4`
- `ctest --test-dir build-nvhpc-cpu --output-on-failure`
- `cmake -S . -B build-cuda-129-ampere86-v2 -DCMAKE_CXX_COMPILER="$PWD/third_party/kokkos/bin/nvcc_wrapper" -DCMAKE_CUDA_COMPILER=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9/bin/nvcc -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++ -DCMAKE_CUDA_ARCHITECTURES=86 -DMPI_CXX_COMPILER=mpicxx -DCUDAToolkit_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9 -DCHANNEL_ENABLE_CUDA=ON -DCHANNEL_ENABLE_HIP=OFF -DKokkos_ARCH_AMPERE86=ON`
- `cmake --build build-cuda-129-ampere86-v2 -j 4`
- `CHANNEL_UMPIRE_RESOURCE=DEVICE ctest --test-dir build-cuda-129-ampere86-v2 --output-on-failure`
- `CUDA_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9 cmake -S . -B build-cuda-129-ampere86-timestep -DCMAKE_CXX_COMPILER="$PWD/third_party/kokkos/bin/nvcc_wrapper" -DCMAKE_CUDA_COMPILER=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9/bin/nvcc -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++ -DCMAKE_CUDA_ARCHITECTURES=86 -DMPI_CXX_COMPILER=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/comm_libs/12.9/openmpi4/latest/bin/mpicxx -DCUDAToolkit_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9 -DCHANNEL_ENABLE_CUDA=ON -DCHANNEL_ENABLE_HIP=OFF -DKokkos_ARCH_AMPERE86=ON`
- `CUDA_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9 cmake --build build-cuda-129-ampere86-timestep -j 4`
- `CUDA_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9 CHANNEL_UMPIRE_RESOURCE=DEVICE ctest --test-dir build-cuda-129-ampere86-timestep -R "test_channel_tiny_timestep.*cpp" --output-on-failure`
- `CUDA_ROOT=/opt/Nvidia/nvhpc/Linux_x86_64/25.5/cuda/12.9 CHANNEL_UMPIRE_RESOURCE=DEVICE mpirun -np 1 build-cuda-129-ampere86-timestep/channel-bin/channel --run tests/cpp/tiny_dns.in --steps 3`

The sandboxed CTest path still fails before test code starts because Open MPI
cannot create PMIx/listener sockets there (`socket() failed with errno=1`).
The same CTest command passes outside the sandbox: 38/38 focused tests.

After `module load toolkits/nvhpc/25.5`, the NVHPC module check found
`nvfortran` and `nvc++` 25.5, with the module's `mpicxx` wrapper resolving to
`nvc++`. The `build-nvhpc-cpu` configure, build, and CTest run passed with
Kokkos OpenMP/Serial, Umpire, Kokkos-FFT, and MPI: 38/38 focused tests.

The NVIDIA CUDA build was verified locally on RTX A6000 GPUs with NVHPC 25.5,
CUDA 12.9, Kokkos CUDA `AMPERE86`, Umpire CUDA `DEVICE` allocations,
Kokkos-FFT/cuFFT, and GPU-aware MPI device buffers: 38/38 focused tests.
The fresh CUDA timestep build also passes the one-rank and two-rank parsed
velocity timestep smoke tests and a direct three-step GPU run.
