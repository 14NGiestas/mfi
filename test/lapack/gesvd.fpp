#:include "common.fpp"
#:include "test/lapack/macros/gesvd.fypp"

program test_gesvd
    use iso_fortran_env
    implicit none

    call test_sgesvd
    call test_dgesvd
    call test_cgesvd
    call test_zgesvd

contains

$:test_implement('?gesvd', DEFAULT_TYPES, gesvd)

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