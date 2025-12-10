#:include "common.fpp"
#:include "test/lapack/macros/trtrs.fypp"

program test_trtrs
    use iso_fortran_env
    implicit none

    write(*,*) 'Starting trtrs tests...'
    call test_strtrs
    call test_dtrtrs
    call test_ctrtrs
    call test_ztrtrs
    write(*,*) 'All trtrs tests completed successfully.'

contains

$:test_implement('?trtrs', DEFAULT_TYPES, trtrs)

#:include "test/lapack/test_common.inc"

end program
