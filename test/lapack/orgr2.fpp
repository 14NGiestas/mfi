#:include "common.fpp"
#:include "test/lapack/macros/orgr2.fypp"

program test_orgr2
    use iso_fortran_env
    implicit none
    $:test_run('?orgr2', REAL_TYPES)
contains

$:test_implement('?orgr2', REAL_TYPES, orgr2)

#:include "test/assert.inc"

end program