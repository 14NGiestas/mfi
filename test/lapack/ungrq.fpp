#:include "common.fpp"
#:include "test/lapack/macros/orgqr.fypp"

program test_ungrq
    use iso_fortran_env
    implicit none
    $:test_run('?ungrq', COMPLEX_TYPES)
contains

$:test_implement('?ungrq', COMPLEX_TYPES, orgqr)

#:include "test/assert.inc"

end program