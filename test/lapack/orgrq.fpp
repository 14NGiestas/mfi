#:include "common.fpp"
#:include "test/lapack/macros/orgqr.fypp"

program test_orgrq
    use iso_fortran_env
    implicit none
    $:test_run('?orgrq', REAL_TYPES)
contains

$:test_implement('?orgrq', REAL_TYPES, orgqr)

#:include "test/assert.inc"

end program