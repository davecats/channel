#!/usr/bin/env bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:4

set -euo pipefail

SRC="${SRC:-$(ws_find convvelo)/channel}"
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
INPUT_FILE="${INPUT_FILE:-${SUBMIT_DIR}/dns.in}"
RUN_ROOT="${RUN_ROOT:-${SUBMIT_DIR}/test_2xa100_lu_probe_$(date +%Y%m%d_%H%M%S)}"
BUILD_DIR="${BUILD_DIR:-${SRC}/build}"
CHANNEL_EXE="${CHANNEL_EXE:-${BUILD_DIR}/channel}"
NVHPC_MODULE="${NVHPC_MODULE:-toolkit/nvidia-hpc-sdk/25.3}"

NP="${NP:-8}"
AUTOTUNE_REPEATS="${AUTOTUNE_REPEATS:-1}"
RUN_DNS_MAIN="${RUN_DNS_MAIN:-0}"
RUN_STANDALONE_PROBES="${RUN_STANDALONE_PROBES:-1}"

read_mesh_value() {
  local key="$1"
  local default="$2"
  local value
  value="$(
    awk -F= -v key="${key}" '
      BEGIN { in_mesh = 0 }
      /^\[/ { in_mesh = ($0 == "[mesh]") }
      in_mesh {
        lhs = $1
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", lhs)
        if (lhs == key) {
          rhs = $2
          gsub(/^[[:space:]]+|[[:space:]]+$/, "", rhs)
          print rhs
          exit
        }
      }
    ' "${INPUT_FILE}" 2>/dev/null || true
  )"
  echo "${value:-${default}}"
}

PROBE_NX="${PROBE_NX:-$(read_mesh_value nx 255)}"
PROBE_NY="${PROBE_NY:-$(read_mesh_value ny 250)}"
PROBE_NZ="${PROBE_NZ:-$(read_mesh_value nz 257)}"

load_nvhpc() {
  if [[ -f /etc/profile.d/lmod.sh ]]; then
    # shellcheck disable=SC1091
    . /etc/profile.d/lmod.sh
  fi
  module load "${NVHPC_MODULE}"
}

mkdir -p "${RUN_ROOT}"

{
  echo "run_root=${RUN_ROOT}"
  echo "src=${SRC}"
  echo "input=${INPUT_FILE}"
  echo "build_dir=${BUILD_DIR}"
  echo "channel_exe=${CHANNEL_EXE}"
  echo "nvhpc_module=${NVHPC_MODULE}"
  echo "np=${NP}"
  echo "probe_nx=${PROBE_NX}"
  echo "probe_ny=${PROBE_NY}"
  echo "probe_nz=${PROBE_NZ}"
  echo "autotune_repeats=${AUTOTUNE_REPEATS}"
  echo "run_dns_main=${RUN_DNS_MAIN}"
  echo "run_standalone_probes=${RUN_STANDALONE_PROBES}"
  echo "submit_host=$(hostname)"
  date
} > "${RUN_ROOT}/run_config.txt"

load_nvhpc

cd "${SRC}"
{
  git status
  git rev-parse HEAD
  git diff
} > "${RUN_ROOT}/source_state.log"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"
cmake "${SRC}" > "${RUN_ROOT}/cmake.log" 2>&1
make -j > "${RUN_ROOT}/build.log" 2>&1
ctest --output-on-failure > "${RUN_ROOT}/ctest.log" 2>&1 || true

if [[ ! -x "${CHANNEL_EXE}" ]]; then
  echo "Missing executable: ${CHANNEL_EXE}" >&2
  exit 1
fi

if [[ "${RUN_DNS_MAIN}" == "1" ]]; then
  cd "${SUBMIT_DIR}"
  if [[ -x scripts/run_a100_new.sh ]]; then
    bash scripts/run_a100_new.sh "${INPUT_FILE}" "${RUN_ROOT}/dns_main"
  else
    echo "Skipping DNS benchmark: scripts/run_a100_new.sh not found in ${SUBMIT_DIR}" >&2
  fi
fi

run_channel_autotune_plain() {
  local label="$1"
  local repeats="$2"
  shift 2

  local out_dir="${RUN_ROOT}/${label}"
  mkdir -p "${out_dir}"
  cp "${INPUT_FILE}" "${out_dir}/dns.in"

  echo "==> ${label}"
  (
    cd "${out_dir}"
    rm -f Dati.cart.out
    env \
      CHANNEL_MPI_AUTOTUNE=report \
      CHANNEL_MPI_AUTOTUNE_REPEATS="${repeats}" \
      CHANNEL_EXIT_AFTER_MPI_AUTOTUNE=1 \
      CHANNEL_Y_PIPELINE_TIMING=1 \
      CHANNEL_DISABLE_RESTART_WRITE=1 \
      "$@" \
      mpirun -np "${NP}" "${CHANNEL_EXE}" \
      > stdout.log 2> stderr.log
  )
}

run_channel_autotune_nsys_per_rank() {
  local label="$1"
  local repeats="$2"
  local trace="$3"
  shift 3

  local out_dir="${RUN_ROOT}/${label}"
  mkdir -p "${out_dir}"
  cp "${INPUT_FILE}" "${out_dir}/dns.in"

  echo "==> ${label}"
  (
    cd "${out_dir}"
    rm -f Dati.cart.out
    env \
      CHANNEL_MPI_AUTOTUNE=report \
      CHANNEL_MPI_AUTOTUNE_REPEATS="${repeats}" \
      CHANNEL_EXIT_AFTER_MPI_AUTOTUNE=1 \
      CHANNEL_Y_PIPELINE_TIMING=1 \
      CHANNEL_DISABLE_RESTART_WRITE=1 \
      "$@" \
      mpirun -np "${NP}" bash -lc '
        set -euo pipefail
        if [[ -f /etc/profile.d/lmod.sh ]]; then
          . /etc/profile.d/lmod.sh
        fi
        module load "'"${NVHPC_MODULE}"'"
        host=$(hostname -s)
        rank=${OMPI_COMM_WORLD_RANK:-${PMIX_RANK:-${PMI_RANK:-${SLURM_PROCID:-0}}}}
        exec nsys profile \
          --trace="'"${trace}"'" \
          --sample=none \
          --cpuctxsw=none \
          --stats=true \
          --force-overwrite=true \
          --export=sqlite \
          -o "autotune_${host}_rank_${rank}" \
          "'"${CHANNEL_EXE}"'"
      ' > stdout.log 2> nsys.log
  )
}

run_channel_autotune_single_report() {
  local label="$1"
  local repeats="$2"
  local trace="$3"
  shift 3

  local out_dir="${RUN_ROOT}/${label}"
  mkdir -p "${out_dir}"
  cp "${INPUT_FILE}" "${out_dir}/dns.in"

  echo "==> ${label}"
  (
    cd "${out_dir}"
    rm -f Dati.cart.out
    env \
      CHANNEL_MPI_AUTOTUNE=report \
      CHANNEL_MPI_AUTOTUNE_REPEATS="${repeats}" \
      CHANNEL_EXIT_AFTER_MPI_AUTOTUNE=1 \
      CHANNEL_Y_PIPELINE_TIMING=1 \
      CHANNEL_DISABLE_RESTART_WRITE=1 \
      "$@" \
      nsys profile \
        --trace="${trace}" \
        --sample=none \
        --cpuctxsw=none \
        --stats=true \
        --force-overwrite=true \
        --export=sqlite \
        -o autotune_allranks \
        mpirun -np "${NP}" "${CHANNEL_EXE}" \
      > stdout.log 2> nsys.log
  )
}

# Same executable, same input, same autotune-only path. These are the direct
# profiling-overhead comparison cases.
run_channel_autotune_plain channel_plain_r1 "${AUTOTUNE_REPEATS}"
run_channel_autotune_nsys_per_rank channel_nsys_perrank_nvtx_r1 "${AUTOTUNE_REPEATS}" nvtx
run_channel_autotune_nsys_per_rank channel_nsys_perrank_nvtx_cuda_r1 "${AUTOTUNE_REPEATS}" nvtx,cuda

# This is useful on one node or if Nsight handles the two-node launch. It may
# fail on some systems, so keep the rest of the script alive.
run_channel_autotune_single_report channel_nsys_single_report_nvtx_r1 "${AUTOTUNE_REPEATS}" nvtx || true

# A repeat-only plain case estimates run-to-run jitter without profiler overhead.
if [[ "${AUTOTUNE_REPEATS}" == "1" ]]; then
  run_channel_autotune_plain channel_plain_r3 3
fi

if [[ "${RUN_STANDALONE_PROBES}" == "1" ]]; then
  cd "${SRC}"

  common_probe_env=(
    NVHPC_MODULE="${NVHPC_MODULE}"
    NP="${NP}"
    NPY="${NP}"
    NPXZ=1
    ACTIVE_N=0
    NLINES=0
    NX="${PROBE_NX}"
    NY="${PROBE_NY}"
    NZ="${PROBE_NZ}"
    ITERS=30
    WARMUP=5
  )

  run_probe() {
    local label="$1"
    shift

    echo "==> ${label}"
    env "${common_probe_env[@]}" OUT_DIR="${RUN_ROOT}/${label}" "$@" ./omp_y_lu_pipeline_repro.sh
  }

  # Standalone controls: same communication pattern, increasingly stripped down.
  run_probe lu8_real_sweep_workspace \
    SWEEP_BATCHES=1 KERNEL_MODE=real STORAGE_MODE=workspace
  run_probe lu8_dummy_sweep_workspace \
    SWEEP_BATCHES=1 KERNEL_MODE=dummy STORAGE_MODE=workspace
  run_probe lu8_none_sweep_workspace \
    SWEEP_BATCHES=1 KERNEL_MODE=none STORAGE_MODE=workspace

  # Profile the communication-only LU path. KERNEL_MODE=none keeps the same
  # solver-style workspace and device MPI calls, but removes LU/pack/halo
  # kernels so D/H behavior in Nsight is easier to attribute.
  run_probe lu8_none_profile_b1_nvtx_workspace \
    BATCHES=1 KERNEL_MODE=none STORAGE_MODE=workspace PROFILE=1 TRACE=nvtx ITERS=5 WARMUP=1
  run_probe lu8_none_profile_b1_nvtx_cuda_workspace \
    BATCHES=1 KERNEL_MODE=none STORAGE_MODE=workspace PROFILE=1 TRACE=nvtx,cuda ITERS=5 WARMUP=1
  run_probe lu8_none_profile_b2_nvtx_workspace \
    BATCHES=2 KERNEL_MODE=none STORAGE_MODE=workspace PROFILE=1 TRACE=nvtx ITERS=5 WARMUP=1
  run_probe lu8_none_profile_b2_nvtx_cuda_workspace \
    BATCHES=2 KERNEL_MODE=none STORAGE_MODE=workspace PROFILE=1 TRACE=nvtx,cuda ITERS=5 WARMUP=1

  # Detailed timings around the suspicious batch counts.
  run_probe lu8_real_detail_b1 \
    ITERS=5 WARMUP=1 BATCHES=1 DETAIL=1 KERNEL_MODE=real STORAGE_MODE=workspace
  run_probe lu8_real_detail_b2 \
    ITERS=5 WARMUP=1 BATCHES=2 DETAIL=1 KERNEL_MODE=real STORAGE_MODE=workspace
  run_probe lu8_real_detail_b4 \
    ITERS=5 WARMUP=1 BATCHES=4 DETAIL=1 KERNEL_MODE=real STORAGE_MODE=workspace

  # Two-rank placement controls: inter-node hop and same-node control.
  echo "==> lu2_two_nodes_real_sweep_workspace"
  env \
    NVHPC_MODULE="${NVHPC_MODULE}" \
    MPIRUN_ARGS="--map-by ppr:1:node" \
    NP=2 NPY=2 NPXZ=1 \
    ACTIVE_N=0 NLINES=0 \
    NX="${PROBE_NX}" NY="${PROBE_NY}" NZ="${PROBE_NZ}" \
    ITERS=30 WARMUP=5 \
    SWEEP_BATCHES=1 KERNEL_MODE=real STORAGE_MODE=workspace \
    OUT_DIR="${RUN_ROOT}/lu2_two_nodes_real_sweep_workspace" \
    ./omp_y_lu_pipeline_repro.sh

  echo "==> lu2_one_node_real_sweep_workspace"
  env \
    NVHPC_MODULE="${NVHPC_MODULE}" \
    MPIRUN_ARGS="--map-by ppr:2:node" \
    NP=2 NPY=2 NPXZ=1 \
    ACTIVE_N=0 NLINES=0 \
    NX="${PROBE_NX}" NY="${PROBE_NY}" NZ="${PROBE_NZ}" \
    ITERS=30 WARMUP=5 \
    SWEEP_BATCHES=1 KERNEL_MODE=real STORAGE_MODE=workspace \
    OUT_DIR="${RUN_ROOT}/lu2_one_node_real_sweep_workspace" \
    ./omp_y_lu_pipeline_repro.sh
fi

echo "Wrote ${RUN_ROOT}"
