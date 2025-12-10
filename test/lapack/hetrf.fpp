#:include "common.fpp"
#:include "test/lapack/macros/hetrf.fypp"

program test_hetrf
    use iso_fortran_env
    implicit none

    write(*,'(A)') 'Starting hetrf tests...'
    call test_chetrf
    call test_zhetrf
    write(*,'(A)') 'All hetrf tests completed successfully.'

contains

$:test_implement('?hetrf', COMPLEX_TYPES, hetrf)

#:include "test/lapack/test_common.inc"

end program
