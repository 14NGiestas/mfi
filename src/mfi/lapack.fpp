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
#:endmute
module mfi_lapack
use iso_fortran_env
use f77_lapack
use f77_lapack, only: mfi_lartg => f77_lartg
implicit none

$:mfi_interface('?geqrf',  DEFAULT_TYPES)
$:mfi_interface('?gerqf',  DEFAULT_TYPES)
$:mfi_interface('?getrf',  DEFAULT_TYPES)
$:mfi_interface('?getri',  DEFAULT_TYPES)
$:mfi_interface('?getrs',  DEFAULT_TYPES)
$:mfi_interface('?hetrf',  COMPLEX_TYPES)
$:mfi_interface('?hegv',   COMPLEX_TYPES)
$:mfi_interface('?heevd',  COMPLEX_TYPES)
$:mfi_interface('?gesvd',  DEFAULT_TYPES)
$:mfi_interface('?potrf',  DEFAULT_TYPES)
$:mfi_interface('?potri',  DEFAULT_TYPES)
$:mfi_interface('?potrs',  DEFAULT_TYPES)

contains

$:mfi_implement('?geqrf',  DEFAULT_TYPES, geqrf_gerqf)
$:mfi_implement('?gerqf',  DEFAULT_TYPES, geqrf_gerqf)
$:mfi_implement('?getrf',  DEFAULT_TYPES, getrf)
$:mfi_implement('?getri',  DEFAULT_TYPES, getri)
$:mfi_implement('?getrs',  DEFAULT_TYPES, getrs)
$:mfi_implement('?hetrf',  COMPLEX_TYPES, hetrf)
$:mfi_implement('?hegv',   COMPLEX_TYPES, hegv)
$:mfi_implement('?heevd',  COMPLEX_TYPES, heevd)
$:mfi_implement('?gesvd',  DEFAULT_TYPES, gesvd)
$:mfi_implement('?potrf',  DEFAULT_TYPES, potrf_potri)
$:mfi_implement('?potri',  DEFAULT_TYPES, potrf_potri)
$:mfi_implement('?potrs',  DEFAULT_TYPES, potrs)

    pure subroutine mfi_error(name, info)
        character(*), intent(in) :: name
        integer, intent(in) :: info
        call f77_xerbla(name, info)
    end subroutine

end module
