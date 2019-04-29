F90 = mpifort
INC = /usr/include/
OPT = optimize
OBJ = typedef.o mpi_transpose.o rbmat.o ffts.o dnsdata.o channel.o

ifeq ($(OPT),check)
  flags =  -malign-double -fall-intrinsics -ffree-line-length-none -I$(INC) -cpp -g -fcheck=all -Wall -Wextra -fbacktrace -ffpe-trap=invalid,zero,overflow -lfftw3 #-fopenmp
else ifeq ($(OPT),profile)
  flags =  -malign-double -fall-intrinsics -ffree-line-length-none -I$(INC) -cpp -finline-functions -g -pg -lfftw3 #-fopenmp
else ifeq ($(OPT),googleprofile)
  flags =  -malign-double -fall-intrinsics -ffree-line-length-none -I$(INC) -cpp -Ofast -lfftw3 -lprofiler #-fopenmp
else
  flags =  -malign-double -fall-intrinsics -ffree-line-length-none  -I$(INC) -cpp -lfftw3  -Ofast #-fopenmp
endif
					 

channel: $(OBJ)
	$(F90) $(flags) -o $@ $(OBJ)
%.o : %.f90
	$(F90) $(flags) -c $<
clean: 
	rm *.mod *.o *.gcda *.gcno *.gcov gmon*
