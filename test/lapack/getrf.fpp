#:include "common.fpp"
#:include "test/lapack/macros/getrf.fypp"

program test_getrf
    use iso_fortran_env
    implicit none

    write(*,'(A)') 'Starting getrf tests...'
    call test_sgetrf
    call test_dgetrf
    call test_cgetrf
    call test_zgetrf
    write(*,'(A)') 'All getrf tests completed successfully.'

contains

$:test_implement('?getrf', DEFAULT_TYPES, getrf)

#:include "test/lapack/test_common.inc"

end program
