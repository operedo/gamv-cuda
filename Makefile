#all: gamv-omp_instr.exe
#all: gamv-omp.exe
#all: gamv.exe

CC=gcc
FC=gfortran
NVCC=nvcc
#FFLAGS= -cpp -Wall -pedantic -O3
FFLAGS= -g -cpp -O3
CFLAGS= -g -O3 
#NVFLAGS= -ccbin=/Soft/gcc/4.7.2/bin/gcc -m64 -c -O3 -arch=sm_35
NVFLAGS= -g -m64 -c -O3 -arch=sm_35

#OPENMP=-fopenmp
#INCS= -I.
#LIBS= -lgomp -lm
#else
#OPENMP=
#INCS= -I.
#LIBS= -lm
#endif

#ifdef CUDA
OPENMP= -fopenmp
INCS= -I. 
NVINCS= -I/usr/local/cuda-8.0/include
#NVINCS= -I/Soft/cuda/6.5.14/include/
OPENMPLIBS= -lgomp 
LIBS= -lm 
#NVLIBS= -L/usr/local/cuda-6.5/lib64 -lcuda -lcudart /usr/local/cuda-6.5/lib64/libcudadevrt.a
NVLIBS= -L/usr/local/cuda-8.0/lib64 -lcuda -lcudart 
#NVLIBS= -L/Soft/cuda/6.5.14/lib64 -lcuda -lcudart 
#NVLIBS= /Soft/cuda/6.5.14/lib64/libcudart_static.a -L/Soft/cuda/6.5.14/lib64 -lcuda 
#endif

#FC=ifort
#FFLAGS=-fpp -O3
#OPENMP=-openmp

F_OBJECTS= extractStatisticsFortran.o
C_OBJECTS= extractStatisticsCwrapper.o 
CU_OBJECTS= extractStatisticsCUDAwrapper.o  
CUOMP_OBJECTS= extractStatisticsCUDAOMPwrapper.o  

GSLIB=../gslib-alges/gslib/gslib.a 

all: seq-for seq-c cuda cudaomp

#seq-for: gslib $(F_OBJECTS)
seq-for: gslib
	$(FC) -DFORTRAN -c $(FFLAGS) gamv.for
	$(FC) -c $(FFLAGS) extractStatisticsFortran.for
	$(FC) -DFORTRAN $(FFLAGS) $(INCS) gamv.o $(F_OBJECTS) -o gamvSeqFor.exe $(GSLIB) $(LIBS) 

#seq-c: gslib $(F_OBJECTS) $(C_OBJECTS)
seq-c: gslib
	$(FC) -DANSIC -c $(FFLAGS) gamv.for
	$(CC) -c $(CFLAGS) extractStatisticsCwrapper.c
	$(FC) -DANSIC $(FFLAGS) $(INCS) gamv.o $(C_OBJECTS) -o gamvSeqC.exe $(GSLIB) $(LIBS) 

#par-for: gslib $(F_OBJECTS)
par-for: gslib
	$(FC) -DFORTRAN -c $(FFLAGS) $(OPENMP) gamv.for
	$(FC) -c $(FFLAGS) $(OPENMP) extractStatisticsFortran.for
	$(FC) -DFORTRAN $(FFLAGS) $(INCS) gamv.o $(F_OBJECTS) -o gamvParFor.exe $(GSLIB) $(OPENMPLIBS) $(LIBS) 

#cuda: gslib $(F_OBJECTS) $(C_OBJECTS) $(CU_OBJECTS)
cuda: gslib
	$(FC) -DCUDA -c $(FFLAGS) gamv.for
	$(NVCC)  $(NVFLAGS)  extractStatisticsCUDAwrapper.cu
#	$(NVCC)  $(CU_OBJECTS) -arch=sm_35 -lcudadevrt -dlink -o extractStatisticsCUDAwrapper.o
#	$(FC) -DCUDA $(FFLAGS) $(INCS) $(NVINCS) gamv.o dlink.o -o gamv.exe $(GSLIB) $(NVLIBS) $(LIBS)
	$(FC) -DCUDA $(FFLAGS) $(INCS) $(NVINCS) gamv.o $(CU_OBJECTS) -o gamvCUDA.exe $(GSLIB) $(NVLIBS) $(LIBS) 

cudaomp: gslib
	$(FC) -DCUDAOMP -c $(FFLAGS) $(OPENMP) gamv.for
	$(NVCC)  $(NVFLAGS)  -Xcompiler -fopenmp extractStatisticsCUDAOMPwrapper.cu
#	$(NVCC)  $(CU_OBJECTS) -arch=sm_35 -lcudadevrt -dlink -o extractStatisticsCUDAwrapper.o
#	$(FC) -DCUDA $(FFLAGS) $(INCS) $(NVINCS) gamv.o dlink.o -o gamv.exe $(GSLIB) $(NVLIBS) $(LIBS)
	$(FC) -DCUDAOMP $(FFLAGS) $(INCS) $(NVINCS) gamv.o $(CUOMP_OBJECTS) -o gamvCUDAOMP.exe $(GSLIB) $(OPENMPLIBS) $(NVLIBS) $(LIBS) 


gslib:
	cd ../gslib-alges/gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP=" "; cd -

#	cd ../gslib90/gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP=" "; cd ../../pargamv
#	cd ../gslib90/gslib; make clean; make gslib.a COMPILER="$(FC)" FLAGS="$(FFLAGS)" OMP="$(OPENMP)"; cd ../../pargamv

clean:
	rm *.exe *.o 

#.SUFFIXES: .o .for .c .cu
#
#ifdef openmp
#.for.o : ; $(FC) -c $(FFLAGS) $(OPENMP) $*.for
#else
#.for.o : ; $(FC) -c $(FFLAGS) $*.for
#endif
#
#.c.o : ; $(CC) -c $(CFLAGS) $*.c
#.cu.o : ; $(NVCC) -c $(NVFLAGS) $*.cu 
