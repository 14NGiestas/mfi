
#:include "common.fpp"
#:include "test/blas/macros/rotg.fypp"

program test_rotg
use iso_fortran_env
implicit none
$:test_run('?rotg', DEFAULT_TYPES)
contains
$:test_implement('?rotg', DEFAULT_TYPES, rotg)

#:include "test/assert.inc"

end program

