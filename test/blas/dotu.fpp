
#:include "common.fpp"
#:include "test/blas/macros/dot_product.fypp"

program test_dotu
use iso_fortran_env
implicit none
$:test_run('?dotu', COMPLEX_TYPES)
contains
$:test_implement('?dotu', COMPLEX_TYPES, dot_product)

#:include "test/assert.inc"

end program

