#:include "common.fpp"
#:include "test/lapack/macros/getri.fypp"

program test_getri
    use iso_fortran_env
    implicit none

    write(*,'(A)') 'Starting getri tests...'
    call test_sgetri
    call test_dgetri
    call test_cgetri
    call test_zgetri
    write(*,'(A)') 'All getri tests completed successfully.'

contains

$:test_implement('?getri', DEFAULT_TYPES, getri)

#:include "test/lapack/test_common.inc"

end program
