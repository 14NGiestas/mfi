FC = gfortran
FPP = fypp
MODULES := blas
SRCS := $(MODULES:%=%.fpp)
SRCS += $(MODULES:%=%.f90)
OBJS := $(MODULES:%=%.o)

%.f90: %.fpp
	$(FPP) $< $@

%.o: %.f90
	$(FC) -c $< -lblas -llapack -lpthread

test_%: test_blas.f90 %.o
	$(FC) $^ -o $@ -lblas -llapack -lpthread
