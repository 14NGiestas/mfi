
#:include "common.fpp"
#:include "test/blas/macros/gemv.fypp"

program test_gemv
use iso_fortran_env
implicit none
$:test_run('?gemv', DEFAULT_TYPES)
contains
$:test_implement('?gemv', DEFAULT_TYPES, gemv)

#:include "test/assert.inc"

end program

