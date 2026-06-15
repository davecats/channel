#!/usr/bin/env bash
#PBS -N schur_alltoallv_repro
#PBS -l select=1:node_type=mi300a:mpiprocs=4
#PBS -l walltime=00:20:00

set -euo pipefail

SRC="${SRC:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
BUILD_DIR="${BUILD_DIR:-${SRC}/build-schur-alltoallv-repro}"
EXE="${EXE:-${BUILD_DIR}/schur_alltoallv_repro}"

FC="${FC:-}"
FFLAGS="${FFLAGS:-}"
LDFLAGS="${LDFLAGS:-}"

DO_BUILD="${DO_BUILD:-1}"
DO_RUN="${DO_RUN:-1}"
PROFILE="${PROFILE:-0}"

NP="${NP:-4}"
PPN="${PPN:-4}"
NX="${NX:-511}"
NY="${NY:-300}"
NZ="${NZ:-513}"
NPY="${NPY:-4}"
ITERS="${ITERS:-20}"
WARMUP="${WARMUP:-2}"
WRAP_USE_DEVICE_ADDR="${WRAP_USE_DEVICE_ADDR:-0}"
EXCHANGE="${EXCHANGE:-alltoallv}"
ALLOC="${ALLOC:-omp}"

CPU_BIND="${CPU_BIND:-list:0-23:24-47:48-71:72-95}"
GPU_BIND="${GPU_BIND:-list:0:1:2:3}"

if [[ -z "${FC}" ]]; then
  if command -v ftn >/dev/null 2>&1; then
    FC=ftn
  elif command -v mpifort >/dev/null 2>&1; then
    FC=mpifort
  else
    FC=mpif90
  fi
fi

if [[ -z "${FFLAGS}" ]]; then
  case "$(basename "${FC}")" in
    ftn)
      # Match the Cray branch in CMakeLists.txt:
      #   -e Z -O2 -hfp0 -homp -lrocprofiler-sdk-roctx
      # plus the target/MPI definitions normally supplied by CMake.
      FFLAGS="-e Z -O2 -hfp0 -homp -lrocprofiler-sdk-roctx -DHAVE_MPI -DHAVE_HIP"
      ;;
    *)
      FFLAGS="-cpp -O2 -fopenmp -DHAVE_MPI"
      ;;
  esac
fi

if [[ -z "${LDFLAGS}" && "${FFLAGS}" == *"-DHAVE_HIP"* ]]; then
  # Standalone addition for the repro's direct hipMalloc/hipFree calls.
  # The main CMake build gets HIP runtime linkage through hipfort targets.
  LDFLAGS="-lamdhip64"
fi

log_cmd() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
  "$@"
}

if [[ "${DO_BUILD}" == "1" ]]; then
  mkdir -p "${BUILD_DIR}"
  log_cmd "${FC}" ${FFLAGS} -I"${SRC}" "${SRC}/schur_alltoallv_repro.F90" ${LDFLAGS} -o "${EXE}"
fi

if [[ "${DO_RUN}" != "1" ]]; then
  echo "Built ${EXE}"
  exit 0
fi

if [[ ! -x "${EXE}" ]]; then
  echo "Missing executable: ${EXE}" >&2
  exit 1
fi

export NX NY NZ NPY ITERS WARMUP WRAP_USE_DEVICE_ADDR EXCHANGE ALLOC

launcher=(
  mpiexec
  -ppn "${PPN}"
  -np "${NP}"
  --line-buffer
  --cpu-bind "${CPU_BIND}"
  --gpu-bind "${GPU_BIND}"
)

app=("${EXE}")
if [[ "${PROFILE}" == "1" ]]; then
  app=(
    rocprofv3
    --sys-trace
    --kernel-trace
    --marker-trace
    --stats
    --output-format=pftrace
    --
    "${EXE}"
  )
fi

echo "schur_alltoallv_repro run:"
echo "  exe=${EXE}"
echo "  NP=${NP} PPN=${PPN} NPY=${NPY}"
echo "  NX=${NX} NY=${NY} NZ=${NZ}"
echo "  ITERS=${ITERS} WARMUP=${WARMUP} WRAP_USE_DEVICE_ADDR=${WRAP_USE_DEVICE_ADDR}"
echo "  EXCHANGE=${EXCHANGE} ALLOC=${ALLOC}"
echo "  PROFILE=${PROFILE}"

log_cmd "${launcher[@]}" "${app[@]}"
