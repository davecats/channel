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

# CC = GCC (uncomment following line for GCC compiler) 
# FLAGS = -cpp -Ofast -malign-double -fall-intrinsics -ffree-line-length-none

# GCC, DEBUG (uncomment following line for GCC compiler)
# FLAGS = -cpp -O0 -malign-double -fall-intrinsics -ffree-line-length-none -fbacktrace -g  -fcheck=all
endef
export CNFGSTRNG
#### 

config_file = compiler.settings
-include ${config_file}

OBJ = typedef.o rbmat.o mpi_transpose.o ffts.o dnsdata.o runtime.o
flags = -I$(FFTW_INC) -L$(FFTW_LIB) $(FLAGS)
libs = -lfftw3

# Target-specific assignments
#############################

# for uiuj (various versions)
postpro/tke/uiuj_largesmall: flags += -DforceblockingY
postpro/tke/uiuj_spectra: flags += -DforceblockingY

# for runtime statistics
runtime: flags += -Druntimestats

postpro/conditional/Velocity_cut: flags += -DforceblockingY
postpro/conditional/zero_crossings: flags += -DforceblockingY

# Actual recipes
################

channel: $(OBJ) channel.o
	$(F90) $(flags) -o  $@ $(OBJ) channel.o $(libs)
	make clean
runtime: $(OBJ) channel.o
	$(F90) $(flags) -o  $@ $(OBJ) channel.o $(libs)
	make clean
out_exporter/out2vtk: $(OBJ) out_exporter/out2vtk.o
	$(F90) $(flags) -o  $@ $(OBJ) out_exporter/out2vtk.o $(libs)
	make clean
out_exporter/out2bin: $(OBJ) out_exporter/out2bin.o
	$(F90) $(flags) -o  $@ $(OBJ) out_exporter/out2bin.o $(libs)
	make clean
postpro/tke/uiuj_largesmall: $(OBJ) postpro/tke/uiuj_largesmall.o
	$(F90) $(flags) -o $@ $(OBJ) postpro/tke/uiuj_largesmall.o $(libs)
	make clean
postpro/tke/uiuj_spectra: $(OBJ) postpro/tke/uiuj_spectra.o
	$(F90) $(flags) -o $@ $(OBJ) postpro/tke/uiuj_spectra.o $(libs)
	make clean
postpro/conditional/Velocity_cut: $(OBJ) postpro/conditional/Velocity_cut.o
	$(CC) $(flags) -o $@ $(OBJ) postpro/conditional/Velocity_cut.o $(libs)
	make clean
postpro/conditional/zero_crossings: $(OBJ) postpro/conditional/zero_crossings.o
	$(CC) $(flags) -o $@ $(OBJ) postpro/conditional/zero_crossings.o $(libs)
	make clean
postpro/am/camstar: $(OBJ) postpro/am/camstar.o
	$(F90) $(flags) -o $@ $(OBJ) postpro/am/camstar.o $(libs)
	make clean
%.o : %.f90
	$(F90) $(flags) -o $@ -c  $<
.PHONY: clean configure
clean:
	find . -type f -name '*.o' | xargs -t rm
	find . -type f -name '*.mod' |xargs -t rm
configure:
	if [ -e ${config_file} ]; then echo "Configuration file already exists."; else echo "$$CNFGSTRNG" > ${config_file}; fi
