#:include "common.fpp"
#:include "test/lapack/macros/geqrf_gerqf.fypp"

program test_gerqf
    use iso_fortran_env
    implicit none
    $:test_run('?gerqf', DEFAULT_TYPES)
contains

$:test_implement('?gerqf', DEFAULT_TYPES, geqrf_gerqf)

#:include "test/assert.inc"

end program