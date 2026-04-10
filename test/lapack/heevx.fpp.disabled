#:include "common.fpp"
#:include "test/lapack/macros/heevx.fypp"

program test_heevx
    use iso_fortran_env
    implicit none
    $:test_run('?heevx', COMPLEX_TYPES)
contains

$:test_implement('?heevx', COMPLEX_TYPES, heevx)

#:include "test/assert.inc"

end program