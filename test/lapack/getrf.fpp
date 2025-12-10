#:include "common.fpp"
#:include "test/lapack/macros/getrf.fypp"

program test_getrf
    use iso_fortran_env
    implicit none

    call test_sgetrf
    call test_dgetrf
    call test_cgetrf
    call test_zgetrf

contains

$:test_implement('?getrf', DEFAULT_TYPES, getrf)

    pure subroutine assert(test, msg, info)
        logical, intent(in) :: test
        character(*), intent(in) :: msg
        integer, intent(in), optional :: info
        character(1024) :: buffer
        if (.not. test) then
            if (present(info)) then
                write(buffer, *) 'Error ', info, ': ', msg
            else
                write(buffer, *) 'Error: ', msg
            end if
            error stop buffer
        end if
    end subroutine

end program