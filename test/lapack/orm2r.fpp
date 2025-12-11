#:include "common.fpp"
#:include "test/lapack/macros/orm2r.fypp"

program test_orm2r
    use iso_fortran_env
    implicit none
    $:test_run('?orm2r', REAL_TYPES)
contains

$:test_implement('?orm2r', REAL_TYPES, orm2r)

#:include "test/assert.inc"

end program