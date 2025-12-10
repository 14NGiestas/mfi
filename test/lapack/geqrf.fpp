#:include "common.fpp"
#:include "test/lapack/macros/geqrf_gerqf.fypp"

program test_geqrf
    use iso_fortran_env
    implicit none

    write(*,*) 'Starting geqrf tests...'
    call test_sgeqrf
    call test_dgeqrf
    call test_cgeqrf
    call test_zgeqrf
    write(*,*) 'All geqrf tests completed successfully.'

contains

$:test_implement('?geqrf', DEFAULT_TYPES, geqrf_gerqf)

#:include "test/lapack/test_common.inc"

end program
