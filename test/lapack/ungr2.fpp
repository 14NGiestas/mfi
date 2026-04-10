#:include "common.fpp"
#:include "test/lapack/macros/ungr2.fypp"

program test_ungr2
    use iso_fortran_env
    implicit none
    $:test_run('?ungr2', COMPLEX_TYPES)
contains

$:test_implement('?ungr2', COMPLEX_TYPES, ungr2)

#:include "test/assert.inc"

end program