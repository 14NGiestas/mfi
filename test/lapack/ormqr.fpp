#:include "common.fpp"
#:include "test/lapack/macros/ormqr.fypp"

program test_ormqr
    use iso_fortran_env
    implicit none
    $:test_run('?ormqr', REAL_TYPES)
contains

$:test_implement('?ormqr', REAL_TYPES, ormqr)

#:include "test/assert.inc"

end program