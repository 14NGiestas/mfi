#:include "common.fpp"
#:include "test/lapack/macros/potrf_potri.fypp"

program test_potrf
    use iso_fortran_env
    implicit none

    call test_spotrf
    call test_dpotrf
    call test_cpotrf
    call test_zpotrf

contains

$:test_implement('?potrf', DEFAULT_TYPES, potrf_potri)

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