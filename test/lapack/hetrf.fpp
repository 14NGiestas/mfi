#:include "common.fpp"
#:include "test/lapack/macros/hetrf.fypp"

program test_hetrf
    use iso_fortran_env
    implicit none

    call test_chetrf
    call test_zhetrf

contains

$:test_implement('?hetrf', COMPLEX_TYPES, hetrf)

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