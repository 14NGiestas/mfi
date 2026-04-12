FPP=fypp

# Source targets
src_fpp_files=$(shell find src -name "*.fpp")
src_f90_files=$(patsubst %.fpp,%.f90,$(src_fpp_files))

# Test targets (umbrella files only)
test_umbrella_files=test/blas.fpp test/lapack.fpp

all: src_split lapack_split test_blas_split test_lapack_split

# Source: split umbrella blas.fpp into per-module .f90 files
src_split: src/mfi/blas.fpp src/f77/blas.fpp
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

# Source: split umbrella lapack.fpp into per-module .f90 files
lapack_split: src/mfi/lapack.fpp src/f77/lapack.fpp
	$(FPP) -m os -I. src/mfi/lapack.fpp > src/mfi/lapack.f90
	$(FPP) -m os -I. src/f77/lapack.fpp > src/f77/lapack.f90

# Test umbrella → split into per-test .f90 files
test_blas_split: test/blas.fpp
	$(FPP) -m os -I. $< > .test_blas.tmp
	csplit -s -z -f test/blas/xx .test_blas.tmp '/^ *program /' '{*}'
	@for f in test/blas/xx*; do \
		name=$$(head -1 $$f | sed 's/.*program \([^ ]*\).*/\1/'); \
		mv "$$f" "test/blas/$$name.f90"; \
	done
	@rm -f .test_blas.tmp test/blas.f90

test_lapack_split: test/lapack.fpp
	$(FPP) -m os -I. $< > .test_lapack.tmp
	csplit -s -z -f test/lapack/xx .test_lapack.tmp '/^ *program /' '{*}'
	@for f in test/lapack/xx*; do \
		name=$$(head -1 $$f | sed 's/.*program \([^ ]*\).*/\1/'); \
		mv "$$f" "test/lapack/$$name.f90"; \
	done
	@rm -f .test_lapack.tmp test/lapack.f90

%.f90: %.fpp
	$(FPP) -I. $< $@

clean:
	-rm -f $(src_f90_files) test/blas/*.f90 test/lapack/*.f90 src/mfi/blas/*.f90 src/f77/blas/*.f90 src/mfi/lapack/*.f90 src/f77/lapack/*.f90 src/mfi/blas.f90 src/f77/blas.f90 src/mfi/lapack.f90 src/f77/lapack.f90 .mfi_* .f77_* .mfi_lapack_* .f77_lapack_* .blas_*.tmp .lapack_*.tmp .test_blas.tmp .test_lapack.tmp

.PHONY: regenerate
regenerate:
	@echo "Forcing regeneration of all .f90 files from .fpp sources..."
	-rm -f $(src_f90_files) test/blas/*.f90 test/lapack/*.f90 src/mfi/blas/*.f90 src/f77/blas/*.f90 src/mfi/lapack/*.f90 src/f77/lapack/*.f90 .mfi_* .f77_* .mfi_lapack_* .f77_lapack_* .blas_*.tmp .lapack_*.tmp .test_blas.tmp .test_lapack.tmp
	$(MAKE) all
