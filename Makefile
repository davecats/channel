.NOTPARALLEL: 

####
F90 = mpifort
FFTW_INC = /path/to/fftw-3.3.8/include/
FFTW_LIB = /path/to/fftw-3.3.8/lib/
# INTEL (uncomment following line for INTEL compiler)
FLAGS = -cpp -Ofast -no-wrap-margin
# GCC (uncomment following line for GCC compiler) 
# FLAGS = -cpp -Ofast -malign-double -fall-intrinsics -ffree-line-length-none 
#### 

OBJ = typedef.o rbmat.o mpi_transpose.o ffts.o dnsdata.o
flags = -I$(FFTW_INC) -L$(FFTW_LIB) $(FLAGS)
libs = -lfftw3

channel: $(OBJ) channel.o
	$(F90) $(flags) -o  $@ $(OBJ) channel.o $(libs)
%.o : %.f90
	$(F90) $(flags) -c  $<
clean: 
	rm *.mod *.o

