.NOTPARALLEL: 

# Fortran + MPI
F90 = mpifort

# Path to FFTW library
FFTW_INC = /usr/include/
FFTW_LIB = /usr/lib/



# COMPILER FLAGS

# CC = GCC (uncomment following line for GCC compiler) 
#FLAGS = -cpp -Ofast -malign-double -fall-intrinsics -ffree-line-length-none

# GCC, DEBUG (uncomment following line for GCC compiler)
FLAGS = -cpp -O0 -malign-double -fall-intrinsics -ffree-line-length-none -fbacktrace -g  -fcheck=all -Wall

OBJ = rbmat.o mpi_transpose.o ffts.o dnsdata.o
flags = -I$(FFTW_INC) -L$(FFTW_LIB) $(FLAGS)
libs = -lfftw3

channel: $(OBJ) channel.o
	$(F90) $(flags) -o  $@ $(OBJ) channel.o $(libs)
	make clean
%.o : %.f90
	$(F90) $(flags) -o $@ -c  $<
.PHONY: clean configure
clean:
	find . -type f -name '*.o' | xargs -t rm
	find . -type f -name '*.mod' |xargs -t rm
