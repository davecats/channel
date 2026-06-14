#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BASELINE_REF="${BASELINE_REF:-0e998519860d2413c99e0e3f264e951979e1a8df}"
INPUT_FILE="${1:-${ROOT_DIR}/dns.in}"
OUT_ROOT="${2:-${ROOT_DIR}/bench_runs_$(date +%Y%m%d_%H%M%S)_compare_since_baseline}"

NP="${NP:-2}"
NPY_LIST="${NPY_LIST:-1 2}"
MODE_LIST="${MODE_LIST:-default}"
PROFILE="${PROFILE:-plain}"
NSYS_TRACE="${NSYS_TRACE:-nvtx,cuda}"
MAX_REGRESSION_PCT="${MAX_REGRESSION_PCT:-}"

CURRENT_BUILD_DIR="${CURRENT_BUILD_DIR:-${ROOT_DIR}/build-nvhpc}"
CURRENT_EXE="${CURRENT_EXE:-${CURRENT_BUILD_DIR}/channel}"

BASELINE_TAG="$(git -C "${ROOT_DIR}" rev-parse --short "${BASELINE_REF}")"
BASELINE_WORKTREE="${BASELINE_WORKTREE:-/tmp/davecats-baseline-${BASELINE_TAG}}"
BASELINE_BUILD_DIR="${BASELINE_BUILD_DIR:-${BASELINE_WORKTREE}/build-nvhpc}"
BASELINE_EXE="${BASELINE_EXE:-${BASELINE_BUILD_DIR}/channel}"

ensure_module() {
  if [[ -f /etc/profile.d/lmod.sh ]]; then
    # shellcheck disable=SC1091
    . /etc/profile.d/lmod.sh
    module load toolkits/nvhpc/25.5
  fi
}

ensure_worktree() {
  if [[ ! -d "${BASELINE_WORKTREE}" ]]; then
    git -C "${ROOT_DIR}" worktree add --detach "${BASELINE_WORKTREE}" "${BASELINE_REF}"
    return
  fi

  local actual_ref
  actual_ref="$(git -C "${BASELINE_WORKTREE}" rev-parse HEAD)"
  if [[ "${actual_ref}" != "$(git -C "${ROOT_DIR}" rev-parse "${BASELINE_REF}")" ]]; then
    echo "Baseline worktree at ${BASELINE_WORKTREE} is not on ${BASELINE_REF}" >&2
    exit 1
  fi
}

ensure_build() {
  local source_dir="$1"
  local build_dir="$2"
  local exe="$3"

  if [[ -x "${exe}" ]]; then
    return
  fi

  cmake -S "${source_dir}" -B "${build_dir}"
  cmake --build "${build_dir}" -j 4
}

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
  local exe="$1"
  local run_root="$2"
  local series="$3"
  local npy="$4"
  local mode="$5"
  local label="np${NP}_npy${npy}_${mode}"
  local run_dir="${run_root}/${series}/${label}"
  local env_args=()

  mkdir -p "${run_dir}"
  write_input "${npy}" "${run_dir}/dns.in"
  mapfile -t env_args < <(mode_env "${mode}")

  echo "==> ${series}/${label}"
  (
    cd "${run_dir}"
    if [[ "${PROFILE}" == "nsys" || "${PROFILE}" == "nvtx" ]]; then
      echo "Running with NSYS profiling (trace=${NSYS_TRACE})..."
      env CHANNEL_DISABLE_RESTART_WRITE=1 "${env_args[@]}" \
        nsys profile --trace="${NSYS_TRACE}" --sample=none --cpuctxsw=none \
          --stats=true --force-overwrite=true --export=sqlite \
          -o "${label}" mpirun -np "${NP}" "${exe}" \
          > stdout.log 2> nsys.log
    else
      env CHANNEL_DISABLE_RESTART_WRITE=1 "${env_args[@]}" \
        mpirun -np "${NP}" "${exe}" \
          > stdout.log 2> stderr.log
    fi
  )
}

extract_post_warmup_avg() {
  local log_file="$1"
  awk '
    /TIME PER TIMESTEP/ {
      step += 1
      if (step > 1) {
        sum += $6
        n += 1
      }
    }
    END {
      if (n == 0) {
        print "nan"
      } else {
        printf "%.6f\n", sum / n
      }
    }
  ' "${log_file}"
}

write_summary() {
  local summary_tsv="$1"
  local summary_md="$2"

  {
    echo -e "case\tcurrent_s\tbaseline_s\tratio\tpercent_delta"
    for npy in ${NPY_LIST}; do
      for mode in ${MODE_LIST}; do
        local label="np${NP}_npy${npy}_${mode}"
        local current_log="${OUT_ROOT}/current/${label}/stdout.log"
        local baseline_log="${OUT_ROOT}/baseline/${label}/stdout.log"
        local current_avg baseline_avg ratio pct

        if [[ ! -f "${current_log}" || ! -f "${baseline_log}" ]]; then
          continue
        fi

        current_avg="$(extract_post_warmup_avg "${current_log}")"
        baseline_avg="$(extract_post_warmup_avg "${baseline_log}")"
        ratio="$(awk -v c="${current_avg}" -v b="${baseline_avg}" 'BEGIN { if (b == 0 || b == "nan" || c == "nan") print "nan"; else printf "%.6f", c / b }')"
        pct="$(awk -v c="${current_avg}" -v b="${baseline_avg}" 'BEGIN { if (b == 0 || b == "nan" || c == "nan") print "nan"; else printf "%.2f", 100.0 * (c - b) / b }')"
        echo -e "${label}\t${current_avg}\t${baseline_avg}\t${ratio}\t${pct}"
      done
    done
  } > "${summary_tsv}"

  {
    echo "| case | current (s) | baseline (s) | ratio | delta |"
    echo "| --- | ---: | ---: | ---: | ---: |"
    awk -F '\t' 'NR > 1 { printf "| %s | %s | %s | %s | %s%% |\n", $1, $2, $3, $4, $5 }' "${summary_tsv}"
  } > "${summary_md}"
}

check_regression_threshold() {
  local summary_tsv="$1"

  if [[ -z "${MAX_REGRESSION_PCT}" ]]; then
    return
  fi

  awk -F '\t' -v limit="${MAX_REGRESSION_PCT}" '
    NR == 1 { next }
    $5 == "nan" { next }
    ($5 + 0.0) > limit {
      printf "Performance regression exceeds %.2f%% for %s: %s%%\n", limit, $1, $5 > "/dev/stderr"
      failed = 1
    }
    END {
      exit failed
    }
  ' "${summary_tsv}"
}

main() {
  ensure_module
  ensure_worktree
  ensure_build "${ROOT_DIR}" "${CURRENT_BUILD_DIR}" "${CURRENT_EXE}"
  ensure_build "${BASELINE_WORKTREE}" "${BASELINE_BUILD_DIR}" "${BASELINE_EXE}"

  mkdir -p "${OUT_ROOT}/current" "${OUT_ROOT}/baseline"

  for npy in ${NPY_LIST}; do
    for mode in ${MODE_LIST}; do
      run_case "${CURRENT_EXE}" "${OUT_ROOT}" "current" "${npy}" "${mode}"
      run_case "${BASELINE_EXE}" "${OUT_ROOT}" "baseline" "${npy}" "${mode}"
    done
  done

  write_summary "${OUT_ROOT}/summary.tsv" "${OUT_ROOT}/summary.md"
  check_regression_threshold "${OUT_ROOT}/summary.tsv"
  echo "Wrote comparison to ${OUT_ROOT}"
  echo "Summary: ${OUT_ROOT}/summary.md"
}

main "$@"
