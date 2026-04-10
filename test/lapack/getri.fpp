#:include "common.fpp"
#:include "test/lapack/macros/getri.fypp"

program test_getri
    use iso_fortran_env
    implicit none
    $:test_run('?getri', DEFAULT_TYPES)
contains

$:test_implement('?getri', DEFAULT_TYPES, getri)

#:include "test/assert.inc"

end program
