#:include "common.fpp"
#:include "test/lapack/macros/sytrf.fypp"

program test_sytrf
    use iso_fortran_env
    implicit none
    $:test_run('?sytrf', REAL_TYPES)
contains

$:test_implement('?sytrf', REAL_TYPES, sytrf)

#:include "test/lapack/test_common.inc"

end program
