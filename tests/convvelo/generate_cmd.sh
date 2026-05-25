# Generate the test file from the python demo
module load toolkits/nvhpc/25.5 
mkdir -p tests/convvelo 
cp tests/data/dns_test_scalar.in tests/convvelo/dns.in 
cp tests/data/start_field_scalar.out tests/convvelo/Dati.cart.out 
cd tests/convvelo 
sed -i 's/^30\s\+-1\s\+7000\s\+\.TRUE\./1 -1 7000 .TRUE./' dns.in 
sed -i 's/3\s\+! nstep/10   ! nstep/' dns.in
mpirun -n 2 ../../build/channel
mpirun -n 1 ../../build/post_pressure
# Rename output files from *.dat to *.fld
for file in *.dat; do
    mv "$file" "${file%.dat}.fld"
done
../../.venv/bin/python ../../export_scalar_raw_statistics.py .
rm *.fld
