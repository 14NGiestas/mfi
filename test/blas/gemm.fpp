
#:include "common.fpp"
#:include "test/blas/macros/gemm.fypp"

program test_gemm
use iso_fortran_env
implicit none
$:test_run('?gemm', DEFAULT_TYPES)
contains
$:test_implement('?gemm', DEFAULT_TYPES, gemm)

#:include "test/assert.inc"

end program

