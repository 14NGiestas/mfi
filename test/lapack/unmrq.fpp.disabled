#:include "common.fpp"
#:include "test/lapack/macros/ormqr.fypp"

program test_unmrq
    use iso_fortran_env
    implicit none
    $:test_run('?unmrq', COMPLEX_TYPES)
contains

$:test_implement('?unmrq', COMPLEX_TYPES, ormqr)

#:include "test/assert.inc"

end program