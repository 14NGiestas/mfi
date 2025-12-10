#:include "common.fpp"
#:include "test/lapack/macros/gesvd.fypp"

program test_gesvd
    use iso_fortran_env
    implicit none
    $:test_run('?gesvd', DEFAULT_TYPES)
contains

$:test_implement('?gesvd', DEFAULT_TYPES, gesvd)

#:include "test/lapack/test_common.inc"

end program
