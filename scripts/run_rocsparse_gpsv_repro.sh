#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
src="${repo_root}/tools/rocsparse_gpsv_interleaved_batch_repro.cpp"
out="${ROCSparse_GPSV_REPRO_BIN:-${repo_root}/build/rocsparse_gpsv_interleaved_batch_repro}"
hipcc_bin="${HIPCC:-hipcc}"

mkdir -p "$(dirname "${out}")"

"${hipcc_bin}" -O2 -std=c++17 "${src}" -lhipsparse -o "${out}"

echo "built ${out}"
echo
echo "=== below suspected boundary: ny=510-ish, m=253, batch=527872 ==="
"${out}" --m 253 --batch-count 527872 "$@"

echo
echo "=== above suspected boundary: ny=520-ish, m=258, batch=527872 ==="
"${out}" --m 258 --batch-count 527872 "$@"
