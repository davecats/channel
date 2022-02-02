.NOTPARALLEL: 

####
define CNFGSTRNG
# COMPILER
##########
# Fortran + MPI
F90 = mpifort


# IMPORTANT: PATH OF FFTW
#########################
# This might not be needed, if your fortran compiler can locate FFTW on its own
FFTW_INC = /path/to/fftw-3.3.8/include/
FFTW_LIB = /path/to/fftw-3.3.8/lib/



# COMPILER FLAGS
################

# INTEL (uncomment following line for INTEL compiler)
FLAGS = -cpp -Ofast -no-wrap-margin

# INTEL, DEBUG (uncomment following line for INTEL compiler)
#FLAGS = -cpp -O0 -g -check all -fpe0 -warn -traceback -debug extended

# GCC (uncomment following line for GCC compiler) 
# FLAGS = -cpp -Ofast -malign-double -fall-intrinsics -ffree-line-length-none

# GCC, DEBUG (uncomment following line for GCC compiler)
# FLAGS = -cpp -O0 -malign-double -fall-intrinsics -ffree-line-length-none -fbacktrace -g  -fcheck=all
endef
export CNFGSTRNG
#### 

config_file = compiler.settings
-include ${config_file}

OBJ = typedef.o rbmat.o mpi_transpose.o ffts.o dnsdata.o
flags = -I$(FFTW_INC) -L$(FFTW_LIB) $(FLAGS)
libs = -lfftw3

# Target-specific assignments
postpro/tke/uiuj_largesmall: flags += -DforceblockingY
postpro/tke/uiuj_spectra: flags += -DforceblockingY

channel: $(OBJ) channel.o
	$(F90) $(flags) -o  $@ $(OBJ) channel.o $(libs)
out_exporter/out2vtk: $(OBJ) out_exporter/out2vtk.o
	$(F90) $(flags) -o  $@ $(OBJ) out_exporter/out2vtk.o $(libs)
out_exporter/out2bin: $(OBJ) out_exporter/out2bin.o
	$(F90) $(flags) -o  $@ $(OBJ) out_exporter/out2bin.o $(libs)
postpro/tke/uiuj_largesmall: $(OBJ) postpro/tke/uiuj_largesmall.o
	$(F90) $(flags) -o $@ $(OBJ) postpro/tke/uiuj_largesmall.o $(libs)
postpro/tke/uiuj_spectra: $(OBJ) postpro/tke/uiuj_spectra.o
	$(F90) $(flags) -o $@ $(OBJ) postpro/tke/uiuj_spectra.o $(libs)
postpro/am/camstar: $(OBJ) postpro/am/camstar.o
	$(F90) $(flags) -o $@ $(OBJ) postpro/am/camstar.o $(libs)
%.o : %.f90
	$(F90) $(flags) -o $@ -c  $<
.PHONY: clean configure
clean:
	find . -type f -name '*.o' | xargs -t rm
	find . -type f -name '*.mod' |xargs -t rm
configure:
	if [ -e ${config_file} ]; then echo "Configuration file already exists."; else echo "$$CNFGSTRNG" > ${config_file}; fi
