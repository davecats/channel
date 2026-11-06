#!/bin/bash
# Launch channel on HLRS Hunter MI300A nodes with one APU per rank:
# rank i is bound to the cores of NUMA node i%4 and to GPU i%4, i.e.
#   mpiexec -np N --cpu-bind list:0-23:24-47:48-71:72-95 --gpu-bind list:0:1:2:3
# extended to N entries.  Without it the ranks share APU 0's cores, and with
# HSA_XNACK=1 (zero-copy) their arrays land in APU 0's memory, so every other
# GPU reads its data remotely: 3-4x slower steps on more than one rank.
# Usage: hunter_mpiexec.sh -n|-np N [mpiexec options] program [args]
# Set HUNTER_MPIEXEC_NOBIND=1 to launch without binding (for comparisons).
set -euo pipefail
launcher="$(command -v mpiexec || true)"
if [[ -z "$launcher" || "$launcher" == /opt/pbs/* ]]; then
    echo "hunter_mpiexec.sh: mpiexec is '${launcher:-missing}', not Cray PALS; source tools/hunter_env.sh" >&2
    exit 1
fi
args=("$@"); n=""
for ((i = 0; i < ${#args[@]} - 1; i++)); do
    case "${args[i]}" in -n | -np | --np) n="${args[i + 1]}" ;; esac
done
if [[ "${HUNTER_MPIEXEC_NOBIND:-0}" == 1 || -z "$n" ]]; then
    exec "$launcher" "$@"
fi
cpus=(0-23 24-47 48-71 72-95)
cpu_list=""; gpu_list=""
for ((r = 0; r < n; r++)); do
    cpu_list+="${cpu_list:+:}${cpus[r % 4]}"
    gpu_list+="${gpu_list:+:}$((r % 4))"
done
exec "$launcher" --cpu-bind "list:$cpu_list" --gpu-bind "list:$gpu_list" "$@"
