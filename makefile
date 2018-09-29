HOM = $(PWD)
FC  = /usr/local/mpi.intel/openmpi1.6.4/bin/mpif90
OBJ = tokRZ_mpi.o
INC = /usr/local/intel/Compiler/11.1/ifort/mkl/include/fftw
LIB = /usr/local/intel/Compiler/11.1/ifort/mkl/lib/em64t

FFLAGS = -c -w -O3

a.out : $(OBJ)
	$(FC) -i_dynamic -mkl -o a.out -I$(INC) -L$(LIB) -liomp5 -lpthread tokRZ_mpi.o  

tokRZ_mpi.o : $(HOM)/tokRZ_mpi.f90
	$(FC) $(FFLAGS) -mkl -I$(INC) $(HOM)/tokRZ_mpi.f90
clean:
	rm $(OBJ)
