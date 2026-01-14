FPP=fypp
FYPPFLAGS=-DMFI_EXTENSIONS -DMFI_USE_CUBLAS
fpp_files=$(shell find test src -name "*.fpp")
f90_files=$(patsubst %.fpp,%.f90,$(fpp_files))
mod_files=$(patsubst %.fpp,%.mod,$(fpp_files))
all: $(f90_files)
%.f90: %.fpp; $(FPP) $(FYPPFLAGS) -I. $< $@
clean:; -rm -f $(f90_files) $(mod_files)
