#:include "common.fpp"
#:include "test/lapack/macros/potrf_potri.fypp"

program test_potrf
    use iso_fortran_env
    implicit none
    $:test_run('?potrf', DEFAULT_TYPES)
contains

$:test_implement('?potrf', DEFAULT_TYPES, potrf_potri)

#:include "test/lapack/test_common.inc"

end program
