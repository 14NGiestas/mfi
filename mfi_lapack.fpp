#:include "common.fpp"

#:def potrf_potri(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, uplo, info)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, a(:,:))
@:optional(character, in, uplo)
@:optional(integer, out, info)
@:localvars(integer, n, lda)
@:defaults(uplo='U')
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(uplo_,n,a,lda,info_)
    if (present(info)) then
        info = info_
    else if (info_ /= 0) then
        !call mfi_error('name', info_)
        error stop info_
    end if
end subroutine
#:enddef

module mfi_lapack
use iso_fortran_env
use f77_lapack
implicit none

$:mfi_interface('?potrf',  DEFAULT_TYPES)
$:mfi_interface('?potri',  DEFAULT_TYPES)

contains

$:mfi_implement('?potrf',  DEFAULT_TYPES, potrf_potri)
$:mfi_implement('?potri',  DEFAULT_TYPES, potrf_potri)

end module
