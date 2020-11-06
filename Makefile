FC = gfortran -O2
FPP = fypp
MODULES := mfi_blas
SRCS := $(MODULES:%=%.fpp)
SRCS += $(MODULES:%=%.f90)
OBJS := $(MODULES:%=%.o)
OBJS += $(MODULES:%=%.mod)

%.f90: %.fpp
	$(FPP) $< $@

%.o: %.f90
	$(FC) -c $< -lblas -llapack -lpthread

%.o: %.mod

test_%: test_mfi_blas.f90 %.o
	$(FC) $^ -o $@ -lblas -llapack -lpthread
