#:include "common.fpp"
#:include "test/lapack/macros/unm2r.fypp"

program test_unm2r
    use iso_fortran_env
    implicit none
    $:test_run('?unm2r', COMPLEX_TYPES)
contains

$:test_implement('?unm2r', COMPLEX_TYPES, unm2r)

#:include "test/assert.inc"

end program