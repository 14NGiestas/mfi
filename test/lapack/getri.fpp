#:include "common.fpp"
#:include "test/lapack/macros/getri.fypp"

program test_getri
    use iso_fortran_env
    implicit none

    call test_sgetri
    call test_dgetri
    call test_cgetri
    call test_zgetri

contains

$:test_implement('?getri', DEFAULT_TYPES, getri)

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