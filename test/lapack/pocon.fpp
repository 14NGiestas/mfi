#:include "common.fpp"
#:include "test/lapack/macros/pocon.fypp"

program test_pocon
    use iso_fortran_env
    implicit none
    $:test_run('?pocon', DEFAULT_TYPES)
contains

$:test_implement('?pocon', DEFAULT_TYPES, pocon)

#:include "test/assert.inc"

end program
