#:include "common.fpp"
#:include "test/lapack/macros/pocon.fypp"

program test_pocon
    use iso_fortran_env
    implicit none

    write(*,'(A)') 'Starting pocon tests...'
    call test_spocon
    call test_dpocon
    call test_cpocon
    call test_zpocon
    write(*,'(A)') 'All pocon tests completed successfully.'

contains

$:test_implement('?pocon', DEFAULT_TYPES, pocon)

#:include "test/lapack/test_common.inc"

end program
