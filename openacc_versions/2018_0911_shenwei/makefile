export PGI_ACC_TIME=1
# fc = /usr/sw-mpp/bin/mpif90 -OPT:IEEE_arith=1
# fc = /usr/sw-mpp/bin/swafort -OPT:IEEE_arith=1 -arrayAnalyse -ldmAnalyse -priv
# fc = /usr/sw-cluster/mpi2/bin/mpiifort
# fc = /usr/local/mpi3-dynamic/bin/mpif90 -diag-disable remark -O2 -parallel
fc = /opt/pgi/linux86-64/2018/mpi/openmpi/bin/mpif90
# fc = /usr/sw-mpp/bin/swafort -OPT:IEEE_arith=1 -dumpcommand mk
# fc = /opt/intel/impi/5.0.3.048/intel64/bin/mpiifort -diag-disable remark

objects = clt_01_modules.o clt_02_program.o clt_03_initialization.o clt_04_calculation.o clt_05_boundary.o clt_06_diagnose.o clt_07_record.o clt_08_mpitransf.o clt_09_interpolation.o clt_10_rmp.o clt_11_else.o

a.out : $(objects)
	$(fc) -o a.out $(objects) -acc 

clt_01_modules.o : clt_01_modules.f90
	$(fc) -c clt_01_modules.f90
clt_02_program.o : clt_02_program.f90
	$(fc) -c clt_02_program.f90
clt_03_initialization.o : clt_03_initialization.f90
	$(fc) -c clt_03_initialization.f90
clt_04_calculation.o : clt_04_calculation.f90 
	$(fc) -c clt_04_calculation.f90 -acc -Minfo
clt_05_boundary.o : clt_05_boundary.f90
	$(fc) -c clt_05_boundary.f90
clt_06_diagnose.o : clt_06_diagnose.f90
	$(fc) -c clt_06_diagnose.f90
clt_07_record.o : clt_07_record.f90
	$(fc) -c clt_07_record.f90
clt_08_mpitransf.o : clt_08_mpitransf.f90
	$(fc) -c clt_08_mpitransf.f90
clt_09_interpolation.o : clt_09_interpolation.f90
	$(fc) -c clt_09_interpolation.f90
clt_10_rmp.o : clt_10_rmp.f90
	$(fc) -c clt_10_rmp.f90
clt_11_else.o : clt_11_else.f90
	$(fc) -c clt_11_else.f90

clean :
	rm a.out $(objects)
