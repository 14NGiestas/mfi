#:include "common.fpp"
#:include "test/blas/macros/asum_nrm2.fypp"

program test_asum
use iso_fortran_env
implicit none
$:test_run('?asum', REAL_TYPES + MIX_REAL_COMPLEX)
contains
$:test_implement('?asum', REAL_TYPES + MIX_REAL_COMPLEX, asum_nrm2)

#:include "test/assert.inc"

end program
