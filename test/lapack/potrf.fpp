#:include "common.fpp"
#:include "test/lapack/macros/potrf_potri.fypp"

program test_potrf
    use iso_fortran_env
    implicit none

    call test_spotrf
    call test_dpotrf
    call test_cpotrf
    call test_zpotrf

contains

$:test_implement('?potrf', DEFAULT_TYPES, potrf_potri)

#:include "test/lapack/test_common.inc"

end program
