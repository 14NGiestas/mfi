#:include "common.fpp"
#:include "test/lapack/macros/orgqr.fypp"

program test_orgqr
    use iso_fortran_env
    implicit none
    $:test_run('?orgqr', REAL_TYPES)
contains

$:test_implement('?orgqr', REAL_TYPES, orgqr)

#:include "test/assert.inc"

end program
