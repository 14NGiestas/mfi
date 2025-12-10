
#:include "common.fpp"
#:include "test/blas/macros/asum_nrm2.fypp"

program test_nrm2
use iso_fortran_env
implicit none
$:test_run('?nrm2', REAL_TYPES + MIX_REAL_COMPLEX)
contains
$:test_implement('?nrm2', REAL_TYPES + MIX_REAL_COMPLEX, asum_nrm2)

#:include "test/assert.inc"

end program

