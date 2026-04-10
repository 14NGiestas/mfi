
#:include "common.fpp"

program test_iamin
use iso_fortran_env
implicit none
#:if defined('MFI_EXTENSIONS')
#:include "test/blas/macros/iamin_iamax.fypp"
$:test_run('i?amin', DEFAULT_TYPES)
contains
$:test_implement('i?amin', DEFAULT_TYPES, iamin_iamax)

#:include "test/assert.inc"

end program
#:else
    write(*,*) 'i?amin tests skipped: extensions not enabled'
end program
#:endif

