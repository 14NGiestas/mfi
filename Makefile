FPP=fypp
fpp_files=$(shell find test src -name "*.fpp")
f90_files=$(patsubst %.fpp,%.f90,$(fpp_files))
all: $(f90_files)
%.f90: %.fpp; $(FPP) -I. $< $@
clean:; -rm $(f90_files)
