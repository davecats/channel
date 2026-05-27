#!/usr/bin/env bash
set -euo pipefail

module load toolkits/nvhpc/25.5

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$repo_root"

mkdir -p tests/convvelo
cp tests/data/dns_test_scalar.in tests/convvelo/dns.in
cp tests/data/start_field_scalar.out tests/convvelo/Dati.cart.out

pushd tests/convvelo >/dev/null
sed -i 's/^30\s\+-1\s\+7000\s\+\.TRUE\./1 -1 7000 .TRUE./' dns.in
sed -i 's/3\s\+! nstep/10   ! nstep/' dns.in
mpirun -n 2 ../../build/channel
mpirun -n 1 ../../build/post_pressure
for file in *.dat; do
  mv "$file" "${file%.dat}.fld"
done
../../.venv/bin/python ../../export_scalar_raw_statistics.py .
rm -f *.fld
popd >/dev/null

rm -f tests/convvelo/raw_statistics.bin tests/convvelo/raw_statistics.bin.fields
rm -f tests/convvelo/convvelo_runtime_minimal.bin tests/convvelo/convvelo_runtime_minimal.bin.fields
.venv/bin/python export_scalar_raw_statistics.py tests/convvelo --output raw_statistics.bin --prepend-means
.venv/bin/python export_scalar_raw_statistics.py tests/convvelo --output convvelo_runtime_minimal.bin --prepend-means --minimal
