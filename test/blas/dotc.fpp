
#:include "common.fpp"
#:include "test/blas/macros/dot_product.fypp"

program test_dotc
use iso_fortran_env
implicit none
$:test_run('?dotc', COMPLEX_TYPES)
contains
$:test_implement('?dotc', COMPLEX_TYPES, dot_product)

#:include "test/assert.inc"

end program

