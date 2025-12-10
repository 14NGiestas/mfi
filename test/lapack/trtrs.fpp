#:include "common.fpp"
#:include "test/lapack/macros/trtrs.fypp"

program test_trtrs
    use iso_fortran_env
    implicit none

    call test_strtrs
    call test_dtrtrs
    call test_ctrtrs
    call test_ztrtrs

contains

$:test_implement('?trtrs', DEFAULT_TYPES, trtrs)

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