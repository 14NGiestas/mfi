#:include "common.fpp"
#:include "test/lapack/macros/unmr2.fypp"

program test_unmr2
    use iso_fortran_env
    implicit none
    $:test_run('?unmr2', COMPLEX_TYPES)
contains

$:test_implement('?unmr2', COMPLEX_TYPES, unmr2)

#:include "test/assert.inc"

end program