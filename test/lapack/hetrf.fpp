#:include "common.fpp"
#:include "test/lapack/macros/hetrf.fypp"

program test_hetrf
    use iso_fortran_env
    implicit none
    $:test_run('?hetrf', COMPLEX_TYPES)
contains
$:test_implement('?hetrf', COMPLEX_TYPES, hetrf)

#:include "test/assert.inc"

end program
