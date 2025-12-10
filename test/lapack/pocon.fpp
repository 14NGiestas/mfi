#:include "common.fpp"
#:include "test/lapack/macros/pocon.fypp"

program test_pocon
    use iso_fortran_env
    implicit none

    call test_spocon
    call test_dpocon
    call test_cpocon
    call test_zpocon

contains

$:test_implement('?pocon', DEFAULT_TYPES, pocon)

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