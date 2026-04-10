
#:include "common.fpp"
#:include "test/blas/macros/trsm.fypp"

program test_trsm
use iso_fortran_env
implicit none
$:test_run('?trsm', DEFAULT_TYPES)
contains
$:test_implement('?trsm', DEFAULT_TYPES, trsm)

#:include "test/assert.inc"

end program
