#:include "common.fpp"
#:include "test/lapack/macros/gesvd.fypp"

program test_gesvd
    use iso_fortran_env
    implicit none

    write(*,'(A)') 'Starting gesvd tests...'
    call test_sgesvd
    call test_dgesvd
    call test_cgesvd
    call test_zgesvd
    write(*,'(A)') 'All gesvd tests completed successfully.'

contains

$:test_implement('?gesvd', DEFAULT_TYPES, gesvd)

#:include "test/lapack/test_common.inc"

end program
