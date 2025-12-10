
#:include "common.fpp"
#:include "test/blas/macros/dot_product.fypp"

program test_dot
use iso_fortran_env
implicit none
$:test_run('?dot',  REAL_TYPES)
contains
$:test_implement('?dot',  REAL_TYPES, dot_product)

#:include "test/assert.inc"

end program

