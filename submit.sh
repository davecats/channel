#!/bin/bash
#SBATCH --time=0-00:15:00
#SBATCH --nodes=128
#SBATCH --ntasks-per-node=20
#SBATCH --job-name=couette-test
#SBATCH --output=log
#SBATCH --error=log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=davide.gatti@kit.edu
#SBATCH --partition=multinode


module load  numlib/mkl/2018
module load  compiler/intel/18.0
module load  numlib/fftw/3.3_serial
module load  mpi/openmpi/3.1

#export OMPI_MCA_btl_openib_allow_ib=1
#export OMPI_MCA_btl_openib_if_include="mlx4_0:1"

mpiexec --bind-to core --map-by core ./channel
