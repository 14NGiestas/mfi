#:mute
#:include "common.fpp"
#:include "src/f77/lapack/aux_lartg.fypp"
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
#:endmute
!> Improved and original F77 interfaces for LAPACK
module f77_lapack
use iso_fortran_env
implicit none

$:f77_interface('?geqrf',  DEFAULT_TYPES, geqrf_gerqf)
$:f77_interface('?gerqf',  DEFAULT_TYPES, geqrf_gerqf)
$:f77_interface('?getrf',  DEFAULT_TYPES, getrf)
$:f77_interface('?getri',  DEFAULT_TYPES, getri)
$:f77_interface('?getrs',  DEFAULT_TYPES, getrs)
$:f77_interface('?hetrf',  COMPLEX_TYPES, hetrf)
$:f77_interface('?hegv',   COMPLEX_TYPES, hegv)
$:f77_interface('?heevd',  COMPLEX_TYPES, heevd)
$:f77_interface('?heevx',  COMPLEX_TYPES, heevx)
$:f77_interface('?heevr',  COMPLEX_TYPES, heevr)
$:f77_interface('?gesvd',  DEFAULT_TYPES, gesvd)
$:f77_interface('?potrf',  DEFAULT_TYPES, potrf_potri)
$:f77_interface('?potri',  DEFAULT_TYPES, potrf_potri)
$:f77_interface('?potrs',  DEFAULT_TYPES, potrs)
$:f77_interface('?pocon',  DEFAULT_TYPES, pocon)
$:f77_interface('?gels',   DEFAULT_TYPES, gels_gelst_getsls)
$:f77_interface('?gelst',  DEFAULT_TYPES, gels_gelst_getsls)
$:f77_interface('?getsls', DEFAULT_TYPES, gels_gelst_getsls)
$:f77_interface('?gelsd',  DEFAULT_TYPES, gelsd)
$:f77_interface('?gelss',  DEFAULT_TYPES, gelss)
$:f77_interface('?gelsy',  DEFAULT_TYPES, gelsy)
$:f77_interface('?gglse',  DEFAULT_TYPES, gglse)
$:f77_interface('?gglsm',  DEFAULT_TYPES, gglsm)

! Other Auxiliary Routines
$:f77_interface('?lartg',  DEFAULT_TYPES, aux_lartg)

    interface f77_xerbla
        pure subroutine xerbla(name,info)
            character(*), intent(in) :: name
            integer, intent(in) :: info
        end subroutine
    end interface f77_xerbla

end module

