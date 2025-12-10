#:include "common.fpp"
#:include "test/lapack/macros/getrf.fypp"

program test_getrf
    use iso_fortran_env
    implicit none
    $:test_run('?getrf', DEFAULT_TYPES)
contains

$:test_implement('?getrf', DEFAULT_TYPES, getrf)

#:include "test/lapack/test_common.inc"

end program
