
#:include "common.fpp"
#:include "test/blas/macros/swap.fypp"

program test_swap
use iso_fortran_env
implicit none
$:test_run('?swap', DEFAULT_TYPES)
contains
$:test_implement('?swap', DEFAULT_TYPES, swap)

#:include "test/assert.inc"

end program

