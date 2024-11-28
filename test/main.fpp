#:include "common.fpp"
#:include "test/blas/asum.fypp"
program main
use iso_fortran_env
implicit none
$:test_run('?asum', DEFAULT_TYPES, & 
    f=lambda pfx: 'sc' if pfx == 'c' else &
                  'dz' if pfx == 'z' else pfx)
contains
$:test_implement('?asum', DEFAULT_TYPES, asum, &
    f=lambda pfx: 'sc' if pfx == 'c' else &
                  'dz' if pfx == 'z' else pfx)

    pure subroutine assert(test, msg)
        logical, intent(in) :: test
        character(*), intent(in) :: msg
        if (.not. test) then
            error stop msg 
        end if
    end subroutine
end program
