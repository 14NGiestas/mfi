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
#:include "src/mfi/lapack/trtrs.fypp"
#:include "src/mfi/lapack/sytrf.fypp"
#:include "src/mfi/lapack/orgqr.fypp"
#:include "src/mfi/lapack/ormqr.fypp"
#:include "src/mfi/lapack/org2r.fypp"
#:include "src/mfi/lapack/orm2r.fypp"
#:include "src/mfi/lapack/orgr2.fypp"
#:include "src/mfi/lapack/ormr2.fypp"
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
    ('?orgqr',  REAL_TYPES,    orgqr),       &
    ('?orgrq',  REAL_TYPES,    orgqr),       &
    ('?ungqr',  COMPLEX_TYPES, orgqr),       &
    ('?ungrq',  COMPLEX_TYPES, orgqr),       &
    ('?ormqr',  REAL_TYPES,    ormqr),       &
    ('?ormrq',  REAL_TYPES,    ormqr),       &
    ('?unmqr',  COMPLEX_TYPES, ormqr),       &
    ('?unmrq',  COMPLEX_TYPES, ormqr),       &
    ('?org2r',  REAL_TYPES,    org2r),       &
    ('?ung2r',  COMPLEX_TYPES, org2r),       &
    ('?orm2r',  REAL_TYPES,    orm2r),       &
    ('?unm2r',  COMPLEX_TYPES, orm2r),       &
    ('?orgr2',  REAL_TYPES,    orgr2),       &
    ('?ungr2',  COMPLEX_TYPES, orgr2),       &
    ('?ormr2',  REAL_TYPES,    ormr2),       &
    ('?unmr2',  COMPLEX_TYPES, ormr2),       &
    ('?potrf',  DEFAULT_TYPES, potrf_potri), &
    ('?potri',  DEFAULT_TYPES, potrf_potri), &
    ('?potrs',  DEFAULT_TYPES, potrs),       &
    ('?pocon',  DEFAULT_TYPES, pocon),       &
    ('?trtrs',  DEFAULT_TYPES, trtrs),       &
    ('?sytrf',  REAL_TYPES,    sytrf),       &
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
