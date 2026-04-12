#:include "common.fpp"
#:include "test/blas/macros/axpy.fypp"
#:include "test/blas/macros/asum_nrm2.fypp"
#:include "test/blas/macros/copy.fypp"
#:include "test/blas/macros/dot_product.fypp"
#:include "test/blas/macros/gemv.fypp"
#:include "test/blas/macros/gemm.fypp"
#:include "test/blas/macros/iamin_iamax.fypp"
#:include "test/blas/macros/lamch.fypp"
#:include "test/blas/macros/rot.fypp"
#:include "test/blas/macros/rotg.fypp"
#:include "test/blas/macros/scal.fypp"
#:include "test/blas/macros/swap.fypp"
#:include "test/blas/macros/trsm.fypp"

#:set _COLLECT = [                                                 &
     ('?axpy', DEFAULT_TYPES,                    axpy),             &
     ('?copy', DEFAULT_TYPES,                    copy),             &
     ('?swap', DEFAULT_TYPES,                    swap),             &
     ('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL, scal),             &
     ('?rot',  DEFAULT_TYPES + MIX_COMPLEX_REAL, rot),              &
     ('?rotg', DEFAULT_TYPES,                    rotg),             &
     ('?dot',  REAL_TYPES,                       dot_product),      &
     ('?dotc', COMPLEX_TYPES,                    dot_product),      &
     ('?dotu', COMPLEX_TYPES,                    dot_product),      &
     ('?asum', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),        &
     ('?nrm2', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),        &
     ('?gemv', DEFAULT_TYPES,                    gemv),             &
     ('?gemm', DEFAULT_TYPES,                    gemm),             &
     ('?trsm', DEFAULT_TYPES,                    trsm),             &
     ('?lamch',REAL_TYPES,                       lamch),            &
]
#:set _modname = lambda name: name.replace('?','')

 #:for name, supported_types, code in _COLLECT
  program ${_modname(name)}$
  use iso_fortran_env
  use mfi_blas
  implicit none
  $:test_run(name, supported_types)
  contains
  $:fmt_time_fn()
  $:test_get_size_fn()
  $:test_implement(name, supported_types, code)

  #:include "test/assert.inc"

  end program

#:endfor

#:for name, supported_types, code in _COLLECT
  program ${_modname(name)}$_gpu
  use iso_fortran_env
  use mfi_blas
  implicit none
  $:test_run_gpu(name, supported_types)
  contains
  $:fmt_time_fn()
  $:test_get_size_fn()
  $:test_implement_gpu(name, supported_types, code)

  #:include "test/assert.inc"

  end program

#:endfor
