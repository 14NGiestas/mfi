#:include "common.fpp"
#:include "test/lapack/macros/orgqr.fypp"

program test_ungqr
    use iso_fortran_env
    implicit none
    $:test_run('?ungqr', COMPLEX_TYPES)
contains

$:test_implement('?ungqr', COMPLEX_TYPES, orgqr)

#:include "test/assert.inc"

end program