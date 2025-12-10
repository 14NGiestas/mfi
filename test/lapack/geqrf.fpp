#:include "common.fpp"
#:include "test/lapack/macros/geqrf_gerqf.fypp"

program test_geqrf
    use iso_fortran_env
    implicit none
    $:test_run('?geqrf', DEFAULT_TYPES)
contains

$:test_implement('?geqrf', DEFAULT_TYPES, geqrf_gerqf)

#:include "test/assert.inc"

end program
