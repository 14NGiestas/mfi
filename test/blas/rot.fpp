
#:include "common.fpp"
#:include "test/blas/macros/rot.fypp"

program test_rot
use iso_fortran_env
implicit none
$:test_run('?rot', DEFAULT_TYPES + MIX_COMPLEX_REAL)
contains
$:test_implement('?rot', DEFAULT_TYPES + MIX_COMPLEX_REAL, rot)

#:include "test/assert.inc"

end program

