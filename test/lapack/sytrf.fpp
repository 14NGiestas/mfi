#:include "common.fpp"
#:include "test/lapack/macros/sytrf.fypp"

program test_sytrf
    use iso_fortran_env
    implicit none

    call test_ssytrf
    call test_dsytrf

contains

$:test_implement('?sytrf', REAL_TYPES, sytrf)

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