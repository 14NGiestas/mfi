#:mute
#:include "common.fpp"
#:include "test/lapack/geqrf_gerqf.fypp"
#:include "test/lapack/getrf.fypp"
#:include "test/lapack/getri.fypp"
#:include "test/lapack/getrs.fypp"
#:include "test/lapack/hetrf.fypp"
#:include "test/lapack/gesvd.fypp"
#:include "test/lapack/hegv.fypp"
#:include "test/lapack/heevd.fypp"
#:include "test/lapack/potrf_potri.fypp"
#:include "test/lapack/potrs.fypp"
#:include "test/lapack/pocon.fypp"
#:include "test/lapack/trtrs.fypp"
#:include "test/lapack/sytrf.fypp"
#:set COLLECT = [                                  &
    ('?geqrf',  DEFAULT_TYPES, geqrf_gerqf),       &
    ('?gerqf',  DEFAULT_TYPES, geqrf_gerqf),       &
    ('?getrf',  DEFAULT_TYPES, getrf),             &
    ('?getri',  DEFAULT_TYPES, getri),             &
    ('?getrs',  DEFAULT_TYPES, getrs),             &
    ('?hetrf',  COMPLEX_TYPES, hetrf),             &
    ('?hegv',   COMPLEX_TYPES, hegv),              &
    ('?heevd',  COMPLEX_TYPES, heevd),             &
    ('?gesvd',  DEFAULT_TYPES, gesvd),             &
    ('?potrf',  DEFAULT_TYPES, potrf_potri),       &
    ('?potri',  DEFAULT_TYPES, potrf_potri),       &
    ('?potrs',  DEFAULT_TYPES, potrs),             &
    ('?pocon',  DEFAULT_TYPES, pocon),             &
    ('?trtrs',  DEFAULT_TYPES, trtrs),             &
    ('?sytrf',  REAL_TYPES,    sytrf),             &
]
#:endmute

program test_mfi_lapack
    use iso_fortran_env
    implicit none

#:for name, supported_types, code in COLLECT
$:test_run(name, supported_types)
#:endfor

contains

#:for name, supported_types, code in COLLECT
$:test_implement(name, supported_types, code)
#:endfor

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
