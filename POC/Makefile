CUDA_PATH=/opt/cuda

all: clean POC 
POC:
	gcc -O2 -I$(CUDA_PATH)/include -c -DCUBLAS_GFORTRAN $(CUDA_PATH)/src/fortran.c
	fypp mfi_blas.fpp mfi_blas.f90
	fypp main.fpp main.f90
	gfortran -c mfi_blas.f90 -lblas
	gfortran -o POC main.f90 fortran.o mfi_blas.o -lblas -L$(CUDA_PATH)/lib64 -lcublas
	./POC
	MFI_USE_CUDA=0 ./POC
	MFI_USE_CUDA=1 ./POC
clean:;@rm -f *.mod *.f90 *.o POC 
