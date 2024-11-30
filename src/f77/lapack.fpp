#:mute
#:include "common.fpp"
#:include "src/f77/lapack/lartg.fypp"
#:include "src/f77/lapack/geqrf_gerqf.fypp"
#:include "src/f77/lapack/getrf.fypp"
#:include "src/f77/lapack/getri.fypp"
#:include "src/f77/lapack/getrs.fypp"
#:include "src/f77/lapack/hetrf.fypp"
#:include "src/f77/lapack/gesvd.fypp"
#:include "src/f77/lapack/hegv.fypp"
#:include "src/f77/lapack/heevx.fypp"
#:include "src/f77/lapack/heevr.fypp"
#:include "src/f77/lapack/heevd.fypp"
#:include "src/f77/lapack/potrf_potri.fypp"
#:include "src/f77/lapack/potrs.fypp"
#:include "src/f77/lapack/pocon.fypp"
#:include "src/f77/lapack/gels_gelst_getsls.fypp"
#:include "src/f77/lapack/gelsd.fypp"
#:include "src/f77/lapack/gelss.fypp"
#:include "src/f77/lapack/gelsy.fypp"
#:include "src/f77/lapack/gglse.fypp"
#:include "src/f77/lapack/gglsm.fypp"
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
    ('?heevx',  COMPLEX_TYPES, heevx),             &
    ('?heevr',  COMPLEX_TYPES, heevr),             &
    ('?gels',   DEFAULT_TYPES, gels_gelst_getsls), &
    ('?gelst',  DEFAULT_TYPES, gels_gelst_getsls), &
    ('?getsls', DEFAULT_TYPES, gels_gelst_getsls), &
    ('?gelsd',  DEFAULT_TYPES, gelsd),             &
    ('?gelss',  DEFAULT_TYPES, gelss),             &
    ('?gelsy',  DEFAULT_TYPES, gelsy),             &
    ('?gglse',  DEFAULT_TYPES, gglse),             &
    ('?gglsm',  DEFAULT_TYPES, gglsm),             &
    ('?lartg',  DEFAULT_TYPES, lartg),             &
]
#:endmute                                          
!> Improved and original F77 interfaces for LAPACK
module f77_lapack
use iso_fortran_env
implicit none


#:for name, supported_types, code in COLLECT
$:f77_original(name, supported_types, code)
#:endfor

#:for name, supported_types, code in COLLECT
$:f77_improved(name, supported_types)
#:endfor

    interface f77_xerbla
        pure subroutine xerbla(name,info)
            character(*), intent(in) :: name
            integer, intent(in) :: info
        end subroutine
    end interface f77_xerbla

end module

