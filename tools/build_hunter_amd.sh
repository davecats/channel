ld channel on HLRS Hunter with the TheRock/amdflang toolchain.
#
# hipfort is built and installed only when it is missing or when the
# toolchain / hipfort source revision changed (or HIPFORT_REBUILD=1).
# channel is (re)configured only when there is no CMake cache yet, when
# extra CMake arguments are passed, or when the cached compiler/launcher no
# longer matches the environment; otherwise only the incremental build runs.
set -eo pipefail

repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
# Shared with interactive shells / job scripts: source it there before ctest.
source "$repo_dir/tools/hunter_env.sh"
set -u

: "${ROCM_PATH:?The therock module must define ROCM_PATH}"
hipfort_source="${HIPFORT_SOURCE:-$HOME/hipfort}"
hipfort_build="${HIPFORT_BUILD_DIR:-$hipfort_source/build-amdflang}"
hipfort_prefix="${HIPFORT_PREFIX:-$HOME/hipfort-afar}"
hipfort_dir="$hipfort_prefix/lib/cmake/hipfort"
hipfort_stamp="$hipfort_prefix/.channel-build-stamp"
build_dir="${CHANNEL_BUILD_DIR:-$repo_dir/build-amd}"
jobs="${BUILD_JOBS:-8}"

if [[ ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "BUILD_JOBS must be a positive integer." >&2
    exit 1
fi

mpi_fc="$(command -v mpif90)"
# Use the same Fortran frontend as MPI, without introducing MPI into hipfort.
read -r hipfort_fc _ <<< "$("$mpi_fc" -show)"
if [[ ! -x "$hipfort_fc" || "$hipfort_fc" != *amdflang ]]; then
    echo "Expected mpif90 -show to start with an absolute amdflang path; got: $hipfort_fc" >&2
    exit 1
fi

mpiexec_bin="$(command -v mpiexec || true)"
if [[ -z "$mpiexec_bin" || "$mpiexec_bin" == /opt/pbs/* ]]; then
    echo "mpiexec resolves to '${mpiexec_bin:-nothing}', not the Cray PALS launcher; load cray-pals." >&2
    exit 1
fi

launcher="$repo_dir/tools/hunter_mpiexec.sh"  # binds rank i to APU i (cores + GPU)

printf 'MPI compiler: %s\nMPI launcher: %s\nhipfort compiler: %s\nhipfort install: %s\n' \
    "$mpi_fc" "$mpiexec_bin" "$hipfort_fc" "$hipfort_prefix"

# ---------------------------------------------------------------- hipfort
hipfort_rev="$(git -C "$hipfort_source" rev-parse HEAD 2>/dev/null || echo unknown)"
expected_stamp="$(printf 'fc=%s\nrocm=%s\nrev=%s\n' "$hipfort_fc" "$ROCM_PATH" "$hipfort_rev")"

if [[ "${HIPFORT_REBUILD:-0}" != 1 && -f "$hipfort_dir/hipfort-config.cmake" \
      && -f "$hipfort_stamp" && "$(cat "$hipfort_stamp")" == "$expected_stamp" ]]; then
    echo "hipfort: up to date, skipping (HIPFORT_REBUILD=1 to force)."
else
    if [[ ! -f "$hipfort_source/CMakeLists.txt" ]]; then
        echo "hipfort source not found: $hipfort_source (set HIPFORT_SOURCE)." >&2
        exit 1
    fi
    echo "hipfort: building and installing."
    # The module exports HIP-specific CXXFLAGS unsuitable for /usr/bin/c++.
    cmake -S "$hipfort_source" -B "$hipfort_build" \
        -DCMAKE_Fortran_COMPILER="$hipfort_fc" \
        -DCMAKE_C_COMPILER=/usr/bin/cc \
        -DCMAKE_CXX_COMPILER=/usr/bin/c++ \
        -DCMAKE_C_FLAGS= -DCMAKE_CXX_FLAGS= \
        -DCMAKE_INSTALL_PREFIX="$hipfort_prefix" \
        -DCMAKE_PREFIX_PATH="$ROCM_PATH" \
        -DBUILD_TESTING=OFF
    cmake --build "$hipfort_build" --parallel "$jobs"
    cmake --install "$hipfort_build"
    printf '%s\n' "$expected_stamp" > "$hipfort_stamp"
fi

# ---------------------------------------------------------------- channel
cache="$build_dir/CMakeCache.txt"
cached() { sed -n "s/^$1:[A-Z]*=//p" "$cache"; }

reconfigure=0
if [[ ! -f "$cache" ]]; then
    reconfigure=1
elif (( $# > 0 )); then
    reconfigure=1
elif [[ "$(cached CMAKE_Fortran_COMPILER)" != "$mpi_fc" \
        || "$(cached MPIEXEC_EXECUTABLE)" != "$launcher" \
        || "$(cached hipfort_DIR)" != "$hipfort_dir" ]]; then
    echo "channel: cached compiler/launcher/hipfort differ from environment, reconfiguring."
    reconfigure=1
fi

if (( reconfigure )); then
    cmake -S "$repo_dir" -B "$build_dir" \
        -DCMAKE_Fortran_COMPILER="$mpi_fc" \
        -DMPIEXEC_EXECUTABLE="$launcher" \
        -Dhipfort_DIR="$hipfort_dir" "$@"
else
    echo "channel: existing configuration is current, building incrementally."
fi
# cmake --build re-runs configure by itself if any CMakeLists.txt changed.
cmake --build "$build_dir" --parallel "$jobs"

