#:include "common.fpp"
#:include "test/lapack/macros/org2r.fypp"

program test_org2r
    use iso_fortran_env
    implicit none
    $:test_run('?org2r', REAL_TYPES)
contains

$:test_implement('?org2r', REAL_TYPES, org2r)

#:include "test/assert.inc"

end program