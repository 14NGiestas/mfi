#:include "common.fpp"
#:include "test/lapack/macros/heevr.fypp"

program test_heevr
    use iso_fortran_env
    implicit none
    $:test_run('?heevr', COMPLEX_TYPES)
contains

$:test_implement('?heevr', COMPLEX_TYPES, heevr)

#:include "test/assert.inc"

end program