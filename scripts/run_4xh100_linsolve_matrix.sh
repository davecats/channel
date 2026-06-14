#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_FILE="${1:-${ROOT_DIR}/dns.in}"
OUT_ROOT="${2:-${ROOT_DIR}/bench_runs_$(date +%Y%m%d_%H%M%S)_4xh100_linsolve}"

CURRENT_BUILD_DIR="${CURRENT_BUILD_DIR:-${ROOT_DIR}/build}"
CURRENT_EXE="${CURRENT_EXE:-${CURRENT_BUILD_DIR}/channel}"
PROFILE="${PROFILE:-nvtx}"
NSYS_TRACE="${NSYS_TRACE:-nvtx}"
CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0,1,2,3}"
CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER:-PCI_BUS_ID}"
NSYS_MULTINODE_MODE="${NSYS_MULTINODE_MODE:-auto}"
NSYS_DISABLE_CUDA_IPC="${NSYS_DISABLE_CUDA_IPC:-0}"
RUN_GPU2="${RUN_GPU2:-1}"

if [[ -f /etc/profile.d/lmod.sh ]]; then
  # shellcheck disable=SC1091
  . /etc/profile.d/lmod.sh
  module load toolkits/nvhpc/25.5
fi

detect_gpu2_root() {
  git -C "${ROOT_DIR}" worktree list --porcelain 2>/dev/null \
    | awk '/^worktree / { wt = $2 } /^branch refs\/heads\/gpu2$/ { print wt; exit }'
}

GPU2_ROOT="${GPU2_ROOT:-$(detect_gpu2_root)}"
GPU2_BUILD_DIR="${GPU2_BUILD_DIR:-${GPU2_ROOT:+${GPU2_ROOT}/build-nvhpc}}"
GPU2_EXE="${GPU2_EXE:-${GPU2_BUILD_DIR:+${GPU2_BUILD_DIR}/channel}}"

if [[ ! -x "${CURRENT_EXE}" ]]; then
  echo "Missing current executable: ${CURRENT_EXE}" >&2
  exit 1
fi
if [[ "${RUN_GPU2}" == "1" && ! -x "${GPU2_EXE}" ]]; then
  echo "Missing gpu2 executable. Set GPU2_EXE=/path/to/gpu2/build-nvhpc/channel or RUN_GPU2=0." >&2
  exit 1
fi

mkdir -p "${OUT_ROOT}"
nvidia-smi -L > "${OUT_ROOT}/gpu_inventory.txt" 2>/dev/null || true

write_input() {
  local npy="$1"
  local dst="$2"
  cp "${INPUT_FILE}" "${dst}"
  if grep -q '^\[parallel\]' "${dst}"; then
    awk -v npy="${npy}" '
      BEGIN { in_parallel = 0; wrote_npy = 0 }
      /^\[/ {
        if (in_parallel && !wrote_npy) {
          print "npy = " npy
          wrote_npy = 1
        }
        in_parallel = ($0 == "[parallel]")
        print
        next
      }
      in_parallel && /^[[:space:]]*npy[[:space:]]*=/ {
        if (!wrote_npy) {
          print "npy = " npy
          wrote_npy = 1
        }
        next
      }
      { print }
      END {
        if (in_parallel && !wrote_npy) print "npy = " npy
      }
    ' "${dst}" > "${dst}.tmp"
    mv "${dst}.tmp" "${dst}"
  else
    {
      echo
      echo "[parallel]"
      echo "npy = ${npy}"
    } >> "${dst}"
  fi
}

run_case() {
  local family="$1"
  local exe="$2"
  local np="$3"
  local npy="$4"
  local label="$5"
  shift 5
  local env_vars=("$@")
  local run_dir="${OUT_ROOT}/${family}_${label}"

  mkdir -p "${run_dir}"
  write_input "${npy}" "${run_dir}/dns.in"

  echo "==> ${family}_${label}"
  (
    cd "${run_dir}"
    if [[ "${PROFILE}" == "nsys" || "${PROFILE}" == "nvtx" ]]; then
      local use_rank_local_nsys=0
      if [[ "${NSYS_MULTINODE_MODE}" == "rank-local" ]]; then
        use_rank_local_nsys=1
      elif [[ "${NSYS_MULTINODE_MODE}" == "auto" && "${np}" -gt 4 ]]; then
        use_rank_local_nsys=1
      fi

      if [[ "${use_rank_local_nsys}" == "1" ]]; then
        local inner_cmd
        inner_cmd="export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES}; \
export CUDA_DEVICE_ORDER=${CUDA_DEVICE_ORDER}; \
if [[ \"${NSYS_DISABLE_CUDA_IPC}\" == \"1\" ]]; then export UCX_TLS=^cuda_ipc; fi; \
env CHANNEL_DISABLE_RESTART_WRITE=1"
        for kv in "${env_vars[@]}"; do
          inner_cmd+=" ${kv}"
        done
        inner_cmd+=" nsys profile --trace=${NSYS_TRACE} --sample=none --cpuctxsw=none"
        inner_cmd+=" --stats=true --force-overwrite=true --export=sqlite"
        inner_cmd+=" -o ${family}_${label}_r"'${OMPI_COMM_WORLD_RANK}'
        inner_cmd+=" ${exe}"

        mpirun -np "${np}" bash -lc "${inner_cmd}" \
          > stdout.log 2> nsys.log
      else
        env CHANNEL_DISABLE_RESTART_WRITE=1 "${env_vars[@]}" \
          nsys profile --trace="${NSYS_TRACE}" --sample=none --cpuctxsw=none \
            --stats=true --force-overwrite=true --export=sqlite \
            -o "${family}_${label}" mpirun -np "${np}" "${exe}" \
            > stdout.log 2> nsys.log
      fi
    else
      env CHANNEL_DISABLE_RESTART_WRITE=1 \
        CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES}" \
        CUDA_DEVICE_ORDER="${CUDA_DEVICE_ORDER}" \
        "${env_vars[@]}" \
        mpirun -np "${np}" "${exe}" \
          > stdout.log 2> stderr.log
    fi
  )
}

{
  echo "out_root=${OUT_ROOT}"
  echo "input=${INPUT_FILE}"
  echo "current_exe=${CURRENT_EXE}"
  echo "gpu2_exe=${GPU2_EXE:-}"
  echo "profile=${PROFILE}"
  echo "nsys_trace=${NSYS_TRACE}"
  echo "cuda_visible_devices=${CUDA_VISIBLE_DEVICES}"
  echo "cuda_device_order=${CUDA_DEVICE_ORDER}"
  echo "nsys_multinode_mode=${NSYS_MULTINODE_MODE}"
  echo "nsys_disable_cuda_ipc=${NSYS_DISABLE_CUDA_IPC}"
} > "${OUT_ROOT}/run_config.txt"

run_case head "${CURRENT_EXE}" 1 1 np1_npy1_local
run_case head "${CURRENT_EXE}" 2 1 np2_npy1_xz
run_case head "${CURRENT_EXE}" 2 2 np2_npy2_schur
run_case head "${CURRENT_EXE}" 4 1 np4_npy1_xz
run_case head "${CURRENT_EXE}" 4 2 np4_npy2_schur
run_case head "${CURRENT_EXE}" 4 4 np4_npy4_schur

  run_case head "${CURRENT_EXE}" 8 1 np8_npy1_xz
  run_case head "${CURRENT_EXE}" 8 4 np8_npy4_schur
  run_case head "${CURRENT_EXE}" 8 8 np8_npy8_schur

if [[ "${RUN_GPU2}" == "1" ]]; then
  run_case gpu2 "${GPU2_EXE}" 1 1 np1_npy1_xz
  run_case gpu2 "${GPU2_EXE}" 2 1 np2_npy1_xz
  run_case gpu2 "${GPU2_EXE}" 4 1 np4_npy1_xz
fi

echo "Wrote runs to ${OUT_ROOT}"
