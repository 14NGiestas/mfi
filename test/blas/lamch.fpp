
#:include "common.fpp"
#:include "test/blas/macros/lamch.fypp"

program test_lamch
use iso_fortran_env
implicit none
$:test_run('?lamch', REAL_TYPES)
contains
$:test_implement('?lamch', REAL_TYPES, lamch)

#:include "test/assert.inc"

end program

