FC = gfortran -O2
FPP = fypp
MODULES := f77_blas mfi_blas f77_lapack mfi_lapack
FPPS = $(MODULES:%=%.fpp)
SRCS = $(MODULES:%=%.f90)
OBJS = $(MODULES:%=%.o)
MODS = $(MODULES:%=%.mod)
.PRECIOUS: %.f90 %.mod %.o

all: $(MODULES)

clean:
	$(RM) *.o *.mod *.f90 *.a *.so

%: %.fpp
	$(FPP) $< $@.f90
	$(FC) -c $@.f90 -lblas -llapack

test_%: all
	$(FPP) $@.fpp $@.f90
	$(FC) $@.f90 *.o -o $@ -lblas -llapack

static-lib: $(OBJS)
	@echo "Creating static library."
	@ar rs libmfi.a $^

shared-lib: $(OBJS)
	@echo "Creating shared library."
	@ld -shared -o libmfi.so $^
