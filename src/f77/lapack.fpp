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
#:endmute
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

! Other Auxiliary Routines
$:f77_interface('?lartg',  DEFAULT_TYPES, aux_lartg)

    interface f77_xerbla
        pure subroutine xerbla(name,info)
            character(*), intent(in) :: name
            integer, intent(in) :: info
        end subroutine
    end interface f77_xerbla

end module

