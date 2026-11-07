#!/bin/bash
# Build channel on HLRS Hunter with the HPE Cray compiler (CCE/crayftn).
#
# hipfort is built and installed only when it is missing or when the
# toolchain / hipfort source revision changed (or HIPFORT_REBUILD=1).  Cray
# .mod files are tied to the CCE version, so this uses its own hipfort prefix
# and build directory.
# channel is (re)configured only when there is no CMake cache yet, when
# extra CMake arguments are passed, or when the cached compiler/launcher no
# longer matches the environment; otherwise only the incremental build runs.
set -eo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
# Shared with interactive shells / job scripts: source it there before ctest.
source "$repo_dir/tools/hunter_cray_env.sh"
set -u
: "${ROCM_PATH:?The rocm module must define ROCM_PATH}"
build_dir="${CHANNEL_BUILD_DIR:-$repo_dir/build-cray}"
jobs="${BUILD_JOBS:-16}"
if [[ ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "BUILD_JOBS must be a positive integer." >&2
    exit 1
fi
# The Cray wrappers add MPI, ROCm and the GTL library themselves.
fc="$(command -v ftn)"
cc_bin="$(command -v cc)"
cxx_bin="$(command -v CC)"
if [[ "$fc" != /opt/cray/pe/craype/* ]]; then
    echo "ftn resolves to '$fc', not the Cray compiler wrapper." >&2
    exit 1
fi
fc_version="$("$fc" --version 2>&1 | head -n 1)"
if [[ "$fc_version" != "Cray Fortran"* ]]; then
    echo "ftn is not Cray Fortran: $fc_version" >&2
    exit 1
fi
mpiexec_bin="$(command -v mpiexec || true)"
if [[ -z "$mpiexec_bin" || "$mpiexec_bin" == /opt/pbs/* ]]; then
    echo "mpiexec resolves to '${mpiexec_bin:-nothing}', not the Cray PALS launcher; load cray-pals." >&2
    exit 1
fi
launcher="$repo_dir/tools/hunter_mpiexec.sh"  # binds rank i to APU i (cores + GPU)
# Cray .mod files are tied to the CCE major version: one hipfort per version.
if [[ "$fc_version" =~ Version\ ([0-9]+)\. ]]; then
    cce_major="${BASH_REMATCH[1]}"
else
    echo "Cannot parse the CCE version from: $fc_version" >&2
    exit 1
fi
hipfort_source="${HIPFORT_SOURCE:-$HOME/hipfort}"
hipfort_build="${HIPFORT_BUILD_DIR:-$hipfort_source/build-cce$cce_major}"
hipfort_prefix="${HIPFORT_PREFIX:-$HOME/hipfort-cce$cce_major}"
hipfort_stamp="$hipfort_prefix/.channel-build-stamp"
printf 'Fortran compiler: %s (%s)\nMPI launcher: %s (via %s)\nhipfort install: %s\n' \
    "$fc" "$fc_version" "$launcher" "$mpiexec_bin" "$hipfort_prefix"
# ---------------------------------------------------------------- hipfort
find_hipfort_dir() {
    local cfg
    cfg="$(find "$hipfort_prefix" -name hipfort-config.cmake -print -quit 2>/dev/null || true)"
    [[ -n "$cfg" ]] && dirname "$cfg"
}
hipfort_rev="$(git -C "$hipfort_source" rev-parse HEAD 2>/dev/null || echo unknown)"
expected_stamp="$(printf 'fc=%s\nversion=%s\nrocm=%s\nrev=%s\n' "$fc" "$fc_version" "$ROCM_PATH" "$hipfort_rev")"
hipfort_dir="$(find_hipfort_dir || true)"
if [[ "${HIPFORT_REBUILD:-0}" != 1 && -n "$hipfort_dir" \
      && -f "$hipfort_stamp" && "$(cat "$hipfort_stamp")" == "$expected_stamp" ]]; then
    echo "hipfort: up to date, skipping (HIPFORT_REBUILD=1 to force)."
else
    if [[ ! -f "$hipfort_source/CMakeLists.txt" ]]; then
        echo "hipfort source not found: $hipfort_source (set HIPFORT_SOURCE)." >&2
        exit 1
    fi
    echo "hipfort: building and installing."
    cmake -S "$hipfort_source" -B "$hipfort_build" \
        -DCMAKE_Fortran_COMPILER="$fc" \
        -DCMAKE_C_COMPILER="$cc_bin" \
        -DCMAKE_CXX_COMPILER="$cxx_bin" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$hipfort_prefix" \
        -DCMAKE_PREFIX_PATH="$ROCM_PATH" \
        -DBUILD_TESTING=OFF
    cmake --build "$hipfort_build" --parallel "$jobs"
    cmake --install "$hipfort_build"
    hipfort_dir="$(find_hipfort_dir || true)"
    if [[ -z "$hipfort_dir" ]]; then
        echo "hipfort-config.cmake not found under $hipfort_prefix after install." >&2
        exit 1
    fi
    printf '%s\n' "$expected_stamp" > "$hipfort_stamp"
fi
# ---------------------------------------------------------------- channel
cache="$build_dir/CMakeCache.txt"
cached() { sed -n "s/^$1:[A-Z]*=//p" "$cache"; }
reconfigure=0
if [[ -f "$cache" && "$(cached CMAKE_HOME_DIRECTORY)" != "$repo_dir" ]]; then
    echo "channel: build dir was configured from another source dir, starting a fresh cache."
    rm -rf "$cache" "$build_dir/CMakeFiles" "$build_dir/_source"
fi
if [[ ! -f "$cache" ]]; then
    reconfigure=1
elif (( $# > 0 )); then
    reconfigure=1
elif [[ "$(cached CMAKE_Fortran_COMPILER)" != "$fc" \
        || "$(cached MPIEXEC_EXECUTABLE)" != "$launcher" \
        || "$(cached hipfort_DIR)" != "$hipfort_dir" ]]; then
    echo "channel: cached compiler/launcher/hipfort differ from environment, reconfiguring."
    reconfigure=1
fi
if (( reconfigure )); then
    cmake -S "$repo_dir" -B "$build_dir" \
        -DCMAKE_Fortran_COMPILER="$fc" \
        -DCMAKE_C_COMPILER="$cc_bin" \
        -DMPIEXEC_EXECUTABLE="$launcher" \
        -Dhipfort_DIR="$hipfort_dir" "$@"
else
    echo "channel: existing configuration is current, building incrementally."
fi
# cmake --build re-runs configure by itself if any CMakeLists.txt changed.
cmake --build "$build_dir" --parallel "$jobs"
