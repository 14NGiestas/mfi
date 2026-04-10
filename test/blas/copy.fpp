
#:include "common.fpp"
#:include "test/blas/macros/copy.fypp"

program test_copy
use iso_fortran_env
implicit none
$:test_run('?copy', DEFAULT_TYPES)
contains
$:test_implement('?copy', DEFAULT_TYPES, copy)

#:include "test/assert.inc"

end program

