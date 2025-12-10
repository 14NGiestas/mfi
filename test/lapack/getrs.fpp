#:include "common.fpp"
#:include "test/lapack/macros/getrs.fypp"

program test_getrs
    use iso_fortran_env
    implicit none

    call test_sgetrs
    call test_dgetrs
    call test_cgetrs
    call test_zgetrs

contains

$:test_implement('?getrs', DEFAULT_TYPES, getrs)

#:include "test/lapack/test_common.inc"

end program
