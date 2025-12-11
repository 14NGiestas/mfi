#:include "common.fpp"
#:include "test/lapack/macros/ormqr.fypp"

program test_ormrq
    use iso_fortran_env
    implicit none
    $:test_run('?ormrq', REAL_TYPES)
contains

$:test_implement('?ormrq', REAL_TYPES, ormqr)

#:include "test/assert.inc"

end program