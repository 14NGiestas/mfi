#:include "common.fpp"
#:include "test/lapack/macros/trtrs.fypp"

program test_trtrs
    use iso_fortran_env
    implicit none
    $:test_run('?trtrs', DEFAULT_TYPES)
contains

$:test_implement('?trtrs', DEFAULT_TYPES, trtrs)

#:include "test/lapack/test_common.inc"

end program
