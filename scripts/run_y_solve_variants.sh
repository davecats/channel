#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${ROOT_DIR}/build-nvhpc}"
CHANNEL_EXE="${CHANNEL_EXE:-${BUILD_DIR}/channel}"
INPUT_FILE="${1:-${ROOT_DIR}/dns.in}"
OUT_ROOT="${2:-${ROOT_DIR}/bench_runs_$(date +%Y%m%d_%H%M%S)_y_solve_variants}"

NP="${NP:-4}"
NPY_LIST="${NPY_LIST:-1 2 4}"
MODE_LIST="${MODE_LIST:-default}"
PROFILE="${PROFILE:-plain}"
NSYS_TRACE="${NSYS_TRACE:-nvtx}"

if [[ ! -x "${CHANNEL_EXE}" ]]; then
  echo "Missing executable: ${CHANNEL_EXE}" >&2
  exit 1
fi

if [[ -f /etc/profile.d/lmod.sh ]]; then
  # shellcheck disable=SC1091
  . /etc/profile.d/lmod.sh
  module load toolkits/nvhpc/25.5
fi

mkdir -p "${OUT_ROOT}"

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

mode_env() {
  local mode="$1"
  case "${mode}" in
    default)
      ;;
    *)
      echo "Unknown mode: ${mode}" >&2
      exit 1
      ;;
  esac
}

run_case() {
  local npy="$1"
  local mode="$2"
  local label="np${NP}_npy${npy}_${mode}"
  local run_dir="${OUT_ROOT}/${label}"
  local env_args=()

  mkdir -p "${run_dir}"
  write_input "${npy}" "${run_dir}/dns.in"
  mapfile -t env_args < <(mode_env "${mode}")

  echo "==> ${label}"
  (
    cd "${run_dir}"
    if [[ "${PROFILE}" == "nsys" || "${PROFILE}" == "nvtx" ]]; then
      env CHANNEL_DISABLE_RESTART_WRITE=1 "${env_args[@]}" \
        nsys profile --trace="${NSYS_TRACE}" --sample=none --cpuctxsw=none \
          --stats=false --force-overwrite=true --export=sqlite \
          -o "${label}" mpirun -np "${NP}" "${CHANNEL_EXE}" \
          > stdout.log 2> nsys.log
    else
      env CHANNEL_DISABLE_RESTART_WRITE=1 "${env_args[@]}" \
        mpirun -np "${NP}" "${CHANNEL_EXE}" \
          > stdout.log 2> stderr.log
    fi
  )
}

for npy in ${NPY_LIST}; do
  for mode in ${MODE_LIST}; do
    run_case "${npy}" "${mode}"
  done
done

echo "Wrote runs to ${OUT_ROOT}"
