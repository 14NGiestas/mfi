
#:include "common.fpp"
#:include "test/blas/macros/axpy.fypp"

program test_axpy
use iso_fortran_env
implicit none
$:test_run('?axpy', DEFAULT_TYPES)
contains
$:test_implement('?axpy', DEFAULT_TYPES, axpy)

#:include "test/assert.inc"

end program

