#:include "common.fpp"
#:include "test/lapack/macros/getrs.fypp"

program test_getrs
    use iso_fortran_env
    implicit none

    call test_sgetrs
    call test_dgetrs
    call test_cgetrs
    call test_zgetrs

contains

$:test_implement('?getrs', DEFAULT_TYPES, getrs)

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