#!/usr/bin/env bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --time=08:00:00
#SBATCH --gres=gpu:4

set -euo pipefail

SRC="${SRC:-$(pwd)}"
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
INPUT_FILE="${INPUT_FILE:-${SUBMIT_DIR}/dns.in}"
RUN_ROOT="${RUN_ROOT:-${SUBMIT_DIR}/channel_comm_b2b_$(date +%Y%m%d_%H%M%S)}"
BUILD_DIR="${BUILD_DIR:-${SRC}/build}"
CHANNEL_EXE="${CHANNEL_EXE:-${BUILD_DIR}/channel}"
NVHPC_MODULE="${NVHPC_MODULE:-toolkits/nvhpc/25.5}"

NP="${NP:-${SLURM_NTASKS:-8}}"
AUTOTUNE_REPEATS="${AUTOTUNE_REPEATS:-3}"
RUN_FORCED="${RUN_FORCED:-1}"
PROFILE="${PROFILE:-1}"
NSYS_TRACE="${NSYS_TRACE:-nvtx,cuda}"
LU_BATCHES="${LU_BATCHES:-4}"

load_toolchain() {
  if [[ -f /etc/profile.d/lmod.sh ]]; then
    # shellcheck disable=SC1091
    . /etc/profile.d/lmod.sh
  fi
  module load "${NVHPC_MODULE}"
}

run_case() {
  local label="$1"
  local comm_backend="$2"
  local np="$3"
  shift 3
  local env_vars=("$@")
  local out_dir="${RUN_ROOT}/${label}_${comm_backend}"

  mkdir -p "${out_dir}"
  cp "${INPUT_FILE}" "${out_dir}/dns.in"

  echo "==> ${label}_${comm_backend}"
  (
    cd "${out_dir}"
    rm -f Dati.cart.out
    env \
      CHANNEL_COMM="${comm_backend}" \
      CHANNEL_DISABLE_RESTART_WRITE=1 \
      "${env_vars[@]}" \
      mpirun -np "${np}" "${CHANNEL_EXE}" \
      > stdout.log 2> stderr.log
  )
}

run_profile_case() {
  local label="$1"
  local comm_backend="$2"
  local np="$3"
  shift 3
  local env_vars=("$@")
  local out_dir="${RUN_ROOT}/${label}_${comm_backend}_nsys"

  mkdir -p "${out_dir}"
  cp "${INPUT_FILE}" "${out_dir}/dns.in"

  echo "==> ${label}_${comm_backend}_nsys"
  (
    cd "${out_dir}"
    rm -f Dati.cart.out
    env \
      CHANNEL_COMM="${comm_backend}" \
      CHANNEL_DISABLE_RESTART_WRITE=1 \
      "${env_vars[@]}" \
      mpirun -np "${np}" bash -lc '
        set -euo pipefail
        if [[ -f /etc/profile.d/lmod.sh ]]; then
          . /etc/profile.d/lmod.sh
        fi
        module load "'"${NVHPC_MODULE}"'"
        host=$(hostname -s)
        rank=${OMPI_COMM_WORLD_RANK:-${PMIX_RANK:-${PMI_RANK:-${SLURM_PROCID:-0}}}}
        exec nsys profile \
          --trace="'"${NSYS_TRACE}"'" \
          --sample=none \
          --cpuctxsw=none \
          --stats=true \
          --force-overwrite=true \
          --export=sqlite \
          -o "channel_${host}_rank_${rank}" \
          "'"${CHANNEL_EXE}"'"
      ' > stdout.log 2> nsys.log
  )
}

run_case_with_optional_profile() {
  local label="$1"
  local comm_backend="$2"
  local np="$3"
  shift 3

  run_case "${label}" "${comm_backend}" "${np}" "$@"
  if [[ "${PROFILE}" == "1" ]]; then
    run_profile_case "${label}" "${comm_backend}" "${np}" "$@"
  fi
}

mkdir -p "${RUN_ROOT}"
load_toolchain

{
  echo "run_root=${RUN_ROOT}"
  echo "src=${SRC}"
  echo "input=${INPUT_FILE}"
  echo "build_dir=${BUILD_DIR}"
  echo "channel_exe=${CHANNEL_EXE}"
  echo "nvhpc_module=${NVHPC_MODULE}"
  echo "np=${NP}"
  echo "autotune_repeats=${AUTOTUNE_REPEATS}"
  echo "run_forced=${RUN_FORCED}"
  echo "profile=${PROFILE}"
  echo "lu_batches=${LU_BATCHES}"
  date
} > "${RUN_ROOT}/run_config.txt"

cd "${SRC}"
{
  git status --short
  git rev-parse HEAD
  git diff
} > "${RUN_ROOT}/source_state.log"

mkdir -p "${BUILD_DIR}"
cmake -S "${SRC}" -B "${BUILD_DIR}" > "${RUN_ROOT}/cmake.log" 2>&1
cmake --build "${BUILD_DIR}" -j > "${RUN_ROOT}/build.log" 2>&1

if [[ ! -x "${CHANNEL_EXE}" ]]; then
  echo "Missing executable: ${CHANNEL_EXE}" >&2
  exit 1
fi

common_autotune_env=(
  CHANNEL_MPI_AUTOTUNE=report
  CHANNEL_MPI_AUTOTUNE_REPEATS="${AUTOTUNE_REPEATS}"
  CHANNEL_EXIT_AFTER_MPI_AUTOTUNE=1
  CHANNEL_Y_PIPELINE_TIMING=1
)

run_case_with_optional_profile autotune mpi "${NP}" "${common_autotune_env[@]}"
run_case_with_optional_profile autotune nccl "${NP}" "${common_autotune_env[@]}"

if [[ "${RUN_FORCED}" == "1" ]]; then
  run_case_with_optional_profile forced_lu_b${LU_BATCHES} mpi "${NP}" \
    CHANNEL_NPXZ=1 CHANNEL_NPY="${NP}" CHANNEL_Y_SOLVER=pipelined_lu CHANNEL_Y_PIPELINE_BATCHES="${LU_BATCHES}"
  run_case_with_optional_profile forced_lu_b${LU_BATCHES} nccl "${NP}" \
    CHANNEL_NPXZ=1 CHANNEL_NPY="${NP}" CHANNEL_Y_SOLVER=pipelined_lu CHANNEL_Y_PIPELINE_BATCHES="${LU_BATCHES}"

  run_case_with_optional_profile forced_xz mpi "${NP}" \
    CHANNEL_NPXZ="${NP}" CHANNEL_NPY=1
  run_case_with_optional_profile forced_xz nccl "${NP}" \
    CHANNEL_NPXZ="${NP}" CHANNEL_NPY=1

  run_case_with_optional_profile forced_schur_alltoall mpi "${NP}" \
    CHANNEL_NPXZ=1 CHANNEL_NPY="${NP}" CHANNEL_Y_SOLVER=schur CHANNEL_Y_SCHUR_EXCHANGE=alltoall
  run_case_with_optional_profile forced_schur_alltoall nccl "${NP}" \
    CHANNEL_NPXZ=1 CHANNEL_NPY="${NP}" CHANNEL_Y_SOLVER=schur CHANNEL_Y_SCHUR_EXCHANGE=alltoall
fi

echo "Wrote runs to ${RUN_ROOT}"
