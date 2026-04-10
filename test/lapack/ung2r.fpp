#:include "common.fpp"
#:include "test/lapack/macros/ung2r.fypp"

program test_ung2r
    use iso_fortran_env
    implicit none
    $:test_run('?ung2r', COMPLEX_TYPES)
contains

$:test_implement('?ung2r', COMPLEX_TYPES, ung2r)

#:include "test/assert.inc"

end program