#:include "common.fpp"
#:include "test/lapack/macros/ormr2.fypp"

program test_ormr2
    use iso_fortran_env
    implicit none
    $:test_run('?ormr2', REAL_TYPES)
contains

$:test_implement('?ormr2', REAL_TYPES, ormr2)

#:include "test/assert.inc"

end program