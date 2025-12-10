
#:include "common.fpp"
#:include "test/blas/macros/scal.fypp"

program test_scal
use iso_fortran_env
implicit none
$:test_run('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL)
contains
$:test_implement('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL, scal)

#:include "test/assert.inc"

end program

