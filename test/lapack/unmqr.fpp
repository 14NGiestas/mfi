#:include "common.fpp"
#:include "test/lapack/macros/ormqr.fypp"

program test_unmqr
    use iso_fortran_env
    implicit none
    $:test_run('?unmqr', COMPLEX_TYPES)
contains

$:test_implement('?unmqr', COMPLEX_TYPES, ormqr)

#:include "test/assert.inc"

end program