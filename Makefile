.NOTPARALLEL: 

####
F90 = mpifort
FFTW_INC = /path/to/fftw-3.3.8/include/
FFTW_LIB = /path/to/fftw-3.3.8/lib/
# INTEL (uncomment following line for INTEL compiler)
FLAGS = -cpp -Ofast -no-wrap-margin
#FLAGS = -cpp -O0 -g -check all -fpe0 -warn -traceback -debug extended
# GCC (uncomment following line for GCC compiler) 
# FLAGS = -cpp -Ofast -malign-double -fall-intrinsics -ffree-line-length-none 
#### 

OBJ = typedef.o rbmat.o mpi_transpose.o ffts.o dnsdata.o
flags = -I$(FFTW_INC) -L$(FFTW_LIB) $(FLAGS)
libs = -lfftw3

channel: $(OBJ) channel.o
	$(F90) $(flags) -o  $@ $(OBJ) channel.o $(libs)
out2vtk: $(OBJ) out_exporter/out2vtk.o
	$(F90) $(flags) -o  out_exporter/$@ $(OBJ) out_exporter/out2vtk.o $(libs)
out2bin: $(OBJ) out_exporter/out2bin.o
	$(F90) $(flags) -o  out_exporter/$@ $(OBJ) out_exporter/out2bin.o $(libs)
postpro/tke/uiuj_largesmall: $(OBJ) postpro/tke/uiuj_largesmall.o
	$(F90) $(flags) -o $@ $(OBJ) postpro/tke/uiuj_largesmall.o $(libs)
%.o : %.f90
	$(F90) $(flags) -c  $<
clean: 
	rm *.mod *.o

