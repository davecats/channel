# Makefile for the GPU-accelerated DNS channel flow code
#
# Usage:
#   make           - Compiles the code
#   make clean     - Removes compiled files
#
# Requires NVIDIA HPC SDK (nvfortran compiler) and an MPI implementation.
#

# Compiler
FC = mpifort

# Compiler flags
# -mp=gpu         : Enable OpenMP GPU offloading
# -Minfo=mp       : Provide feedback on OpenMP parallelization
# -cuda           : Enable CUDA Fortran features (for cuFFT)
# -cudalib=cufft  : Link against the cuFFT library
FFLAGS = -mp=gpu -Minfo=mp -cuda -cudalib=cufft -cpp


# Executable name
TARGET = channel.x

# Source files and objects
SRCS = mpi_transpose.f90 ffts.f90 rbmat.f90  dnsdata.f90 channel.f90
OBJS = $(SRCS:.f90=.o)

# Default target
all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) $(MPI_LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) $(MPI_FFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f *.o *.mod $(TARGET)

.PHONY: all clean

