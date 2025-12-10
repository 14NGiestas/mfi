#:include "common.fpp"
#:include "test/lapack/macros/getrs.fypp"

program test_getrs
    use iso_fortran_env
    implicit none
    $:test_run('?getrs', DEFAULT_TYPES)
contains

$:test_implement('?getrs', DEFAULT_TYPES, getrs)

#:include "test/lapack/test_common.inc"

end program
