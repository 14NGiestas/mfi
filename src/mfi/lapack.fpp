#:mute
#:include "common.fpp"
#:include "src/mfi/lapack/geqrf_gerqf.fypp"
#:include "src/mfi/lapack/getrf.fypp"
#:include "src/mfi/lapack/getri.fypp"
#:include "src/mfi/lapack/getrs.fypp"
#:include "src/mfi/lapack/hetrf.fypp"
#:include "src/mfi/lapack/gesvd.fypp"
#:include "src/mfi/lapack/hegv.fypp"
#:include "src/mfi/lapack/heevd.fypp"
#:include "src/mfi/lapack/potrf_potri.fypp"
#:include "src/mfi/lapack/potrs.fypp"
#:include "src/mfi/lapack/pocon.fypp"
#:set COLLECT = [                            &
    ('?geqrf',  DEFAULT_TYPES, geqrf_gerqf), &
    ('?gerqf',  DEFAULT_TYPES, geqrf_gerqf), &
    ('?getrf',  DEFAULT_TYPES, getrf),       &
    ('?getri',  DEFAULT_TYPES, getri),       &
    ('?getrs',  DEFAULT_TYPES, getrs),       &
    ('?hetrf',  COMPLEX_TYPES, hetrf),       &
    ('?hegv',   COMPLEX_TYPES, hegv),        &
    ('?heevd',  COMPLEX_TYPES, heevd),       &
    ('?gesvd',  DEFAULT_TYPES, gesvd),       &
    ('?potrf',  DEFAULT_TYPES, potrf_potri), &
    ('?potri',  DEFAULT_TYPES, potrf_potri), &
    ('?potrs',  DEFAULT_TYPES, potrs),       &
    ('?pocon',  DEFAULT_TYPES, pocon),       &
]
#:endmute
!> Modern fortran interfaces for LAPACK
module mfi_lapack
use iso_fortran_env
use f77_lapack
use f77_lapack, only: mfi_lartg => f77_lartg
implicit none

#:for name, supported_types, code in COLLECT
$:mfi_interface(name, supported_types)
#:endfor

contains

#:for name, supported_types, code in COLLECT
$:mfi_implement(name, supported_types, code)
#:endfor

    pure subroutine mfi_error(name, info)
        character(*), intent(in) :: name
        integer, intent(in) :: info
        call f77_xerbla(name, info)
    end subroutine

end module
