FPP=fypp
FYPPFLAGS=-DMFI_EXTENSIONS -DMFI_USE_CUBLAS
fpp_files=$(shell find test src -name "*.fpp" ! -name "blas.fpp")
f90_files=$(patsubst %.fpp,%.f90,$(fpp_files))
mod_files=$(patsubst %.fpp,%.mod,$(fpp_files))
all: $(f90_files) src/mfi/blas_*.f90 src/f77/blas_*.f90 src/mfi/blas.f90 src/f77/blas.f90

.PHONY: split_mfi split_f77
split_mfi: src/mfi/blas.fpp
	$(FPP) -m os $(FYPPFLAGS) -I. $< | awk \
		'/^module mfi_blas_/{name=$$2; gsub(/ /,"",name)} \
		 name{print > "src/mfi/"name".f90"} \
		 /^end module/{if(name)close("src/mfi/"name".f90")}'

split_f77: src/f77/blas.fpp
	$(FPP) -m os $(FYPPFLAGS) -I. $< | awk \
		'/^module f77_blas_/{name=$$2; gsub(/ /,"",name)} \
		 name{print > "src/f77/"name".f90"} \
		 /^end module/{if(name)close("src/f77/"name".f90")}'

%.f90: %.fpp
	$(FPP) $(FYPPFLAGS) -I. $< $@

clean:
	-rm -f $(f90_files) $(mod_files) src/mfi/blas_*.f90 src/f77/blas_*.f90 src/mfi/blas.f90 src/f77/blas.f90
