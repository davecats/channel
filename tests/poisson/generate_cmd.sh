# Generate the test file from known good CPL code
module load toolkits/nvhpc/25.5 
mkdir -p tests/poisson 
cp tests/data/dns_test.in tests/poisson/dns.in 
cp tests/data/start_field.out tests/poisson/Dati.cart.out 
cd tests/poisson 
sed -i 's/^30\s\+-1\s\+7000\s\+\.TRUE\./1 -1 7000 .TRUE./' dns.in 
sed -i 's/^2\s*$/100/' dns.in && mpirun -n 2 ../../build/channel
sed -i 's/^1.5\s\+0.0\s\+2.0\s\+! a, ymin, ymax$/1.5    0.0    2.0                     ! a, ymin, ymax\n.FALSE.   1      0.161436                ! CPI, CPItype, gamma/' dns.in 
../../poisson/prepare-pressure 1 1 localhost 35 35 1 $(pwd)/ 
../../poisson/prepare-py 1 1 localhost 35 35 1 $(pwd)/
rm dns.in