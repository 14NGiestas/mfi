#:include "common.fpp"
#:include "test/lapack/macros/geqrf_gerqf.fypp"

program test_geqrf
    use iso_fortran_env
    implicit none

    call test_sgeqrf
    call test_dgeqrf
    call test_cgeqrf
    call test_zgeqrf

contains

$:test_implement('?geqrf', DEFAULT_TYPES, geqrf_gerqf)

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