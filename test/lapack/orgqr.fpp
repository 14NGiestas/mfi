#:include "common.fpp"
#:include "test/lapack/macros/orgqr.fypp"

program test_orgqr
    use iso_fortran_env
    implicit none
    $:test_run('?orgqr', DEFAULT_TYPES)
contains

$:test_implement('?orgqr', DEFAULT_TYPES, orgqr)

#:include "test/assert.inc"

end program