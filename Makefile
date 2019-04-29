#F90 = /home/fh1-project-reeffect/xt8786/openmpi-4.0.1/bin/mpifort
#INC = /home/fh1-project-reeffect/xt8786/openmpi-4.0.1/include/
#LIB = /home/fh1-project-reeffect/xt8786/openmpi-4.0.1/lib/

F90 = ${MPI_ROOT}/bin/mpifort
INC = ${MPI_ROOT}/lib 
LIB = ${MPI_ROOT}/lib

OBJ = typedef.o rbmat.o mpi_transpose.o ffts.o dnsdata.o
#flags = -malign-double -fall-intrinsics -ffree-line-length-none -I$(INC) -cpp -L$(LIB)  -Ofast 
#flags = -malign-double -fall-intrinsics -ffree-line-length-none -I$(INC) -cpp -L$(LIB)  -Og -g -Wall -Wextra -fcheck=all -fbacktrace -ffpe-trap=invalid
flags = -I$(INC) -L$(LIB) -cpp  -Ofast -no-wrap-margin   
libs = -lfftw3 

channel: $(OBJ) channel.o
	$(F90) $(flags) -o  $@ $(OBJ) channel.o $(libs)
%.o : %.f90
	$(F90) $(flags) -c  $<  
clean: 
	rm *.mod *.o

