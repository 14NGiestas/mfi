FPP=fypp
fpp_files=$(shell find test src -name "*.fpp" ! -name "blas.fpp")
f90_files=$(patsubst %.fpp,%.f90,$(fpp_files))
mod_files=$(patsubst %.fpp,%.mod,$(fpp_files))
all: $(f90_files) split

# Umbrella → split into per-module .f90 files with simple names
split: src/mfi/blas.fpp src/f77/blas.fpp
	$(FPP) -m os -I. src/mfi/blas.fpp > .blas_mfi.tmp
	csplit -s -f .mfi_ -z .blas_mfi.tmp '/^module /' '{*}'
	@for f in .mfi_*; do \
		mod=$$(head -1 $$f | sed 's/module \([^ ]*\).*/\1/'); \
		if [ "$$mod" = "mfi_blas" ]; then \
			mv "$$f" "src/mfi/blas.f90"; \
		else \
			case "$$mod" in mfi_blas_*) name=$${mod#mfi_blas_} ;; *) name=$$mod ;; esac; \
			mv "$$f" "src/mfi/blas/$$name.f90"; \
		fi \
	done
	$(FPP) -m os -I. src/f77/blas.fpp > .blas_f77.tmp
	csplit -s -f .f77_ -z .blas_f77.tmp '/^module /' '{*}'
	@for f in .f77_*; do \
		mod=$$(head -1 $$f | sed 's/module \([^ ]*\).*/\1/'); \
		if [ "$$mod" = "f77_blas" ]; then \
			mv "$$f" "src/f77/blas.f90"; \
		else \
			case "$$mod" in f77_blas_*) name=$${mod#f77_blas_} ;; *) name=$$mod ;; esac; \
			mv "$$f" "src/f77/blas/$$name.f90"; \
		fi \
	done
	@rm -f .blas_mfi.tmp .blas_f77.tmp

%.f90: %.fpp
	$(FPP) -I. $< $@

clean:
	-rm -f $(f90_files) $(mod_files) src/mfi/blas/*.f90 src/f77/blas/*.f90 src/mfi/blas.f90 src/f77/blas.f90 .mfi_* .f77_* .blas_*.tmp
