#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT_FILE="${1:-${ROOT_DIR}/dns.in}"
OUT_ROOT="${2:-${ROOT_DIR}/bench_runs_$(date +%Y%m%d_%H%M%S)_4xh100_linsolve}"

CURRENT_BUILD_DIR="${CURRENT_BUILD_DIR:-${ROOT_DIR}/build-nvhpc}"
CURRENT_EXE="${CURRENT_EXE:-${CURRENT_BUILD_DIR}/channel}"
PROFILE="${PROFILE:-nvtx}"
NSYS_TRACE="${NSYS_TRACE:-nvtx}"
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
      env CHANNEL_DISABLE_RESTART_WRITE=1 "${env_vars[@]}" \
        nsys profile --trace="${NSYS_TRACE}" --sample=none --cpuctxsw=none \
          --stats=false --force-overwrite=true --export=sqlite \
          -o "${family}_${label}" mpirun -np "${np}" "${exe}" \
          > stdout.log 2> nsys.log
    else
      env CHANNEL_DISABLE_RESTART_WRITE=1 "${env_vars[@]}" \
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
} > "${OUT_ROOT}/run_config.txt"

run_case head "${CURRENT_EXE}" 1 1 np1_npy1_local
run_case head "${CURRENT_EXE}" 2 1 np2_npy1_xz
run_case head "${CURRENT_EXE}" 2 2 np2_npy2_schur
run_case head "${CURRENT_EXE}" 2 2 np2_npy2_yslab CHANNEL_USE_YSLAB_LINSOLVE=1
run_case head "${CURRENT_EXE}" 4 1 np4_npy1_xz
run_case head "${CURRENT_EXE}" 4 2 np4_npy2_schur
run_case head "${CURRENT_EXE}" 4 2 np4_npy2_yslab CHANNEL_USE_YSLAB_LINSOLVE=1
run_case head "${CURRENT_EXE}" 4 4 np4_npy4_yslab CHANNEL_USE_YSLAB_LINSOLVE=1

if [[ "${RUN_GPU2}" == "1" ]]; then
  run_case gpu2 "${GPU2_EXE}" 1 1 np1_npy1_xz
  run_case gpu2 "${GPU2_EXE}" 2 1 np2_npy1_xz
  run_case gpu2 "${GPU2_EXE}" 4 1 np4_npy1_xz
fi

echo "Wrote runs to ${OUT_ROOT}"
