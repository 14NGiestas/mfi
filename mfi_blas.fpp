#:include "common.fpp"

#:def dotc(NAME,TYPE,KIND)
function mfi_${NAME}$(x, y, incx, incy)
    integer, external :: ${NAME}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in, x(:), y(:))
@:localvars(integer, n, mfi_${NAME}$)
@:optional(integer, incx, incy)
@:defaults(incx=1, incy=1)
    N = size(X)
    mfi_${NAME}$ = ${NAME}$(n,x,incx,y,incy)
end function
#:enddef

#:def axpy(NAME,TYPE,KIND)
subroutine mfi_${NAME}$(x, y, a, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:))
@:args(${TYPE}$, inout, y(:))
@:optional(integer, a, incx, incy)
@:localvars(integer, n)
@:defaults(a=1, incx=1, incy=1)
    N = size(X)
    call ${NAME}$(n,local_a,x,incx,y,incy)
end subroutine
#:enddef

#:def iamin_iamax(NAME,TYPE,KIND)
function mfi_${NAME}$(x, incx)
    integer, external :: ${NAME}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in, x(:))
@:optional(integer, incx)
@:localvars(integer, n, mfi_${NAME}$)
@:defaults(incx=1)
    n = size(x)
    mfi_${NAME}$ = ${NAME}$(n,x,local_incx)
end function
#:enddef

#:def gemv(NAME,TYPE,KIND)
subroutine mfi_${NAME}$(a, x, y, trans, alpha, beta, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in, a(:,:), x(:))
@:args(${TYPE}$, inout, y(:))
@:optional(character, trans)
@:optional(${TYPE}$,  alpha, beta)
@:optional(integer, incx,  incy)
@:localvars(integer, m, n, lda)
@:defaults(trans='N', alpha=1, beta=0, incx=1, incy=1)
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call ${NAME}$(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef

#:def gemm(NAME,TYPE,KIND)
subroutine mfi_${NAME}$(a, b, c, transa, transb, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in, a(:,:), b(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, transa, transb)
@:optional(${TYPE}$,  alpha, beta)
@:localvars(integer, m, n, k, lda, ldb, ldc)
@:defaults(transa='N', transb='N', alpha=1, beta=0)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    call ${NAME}$(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:enddef

#:def herk(NAME,TYPE,KIND)
subroutine mfi_${NAME}$(a, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, trans, uplo)
@:optional(${TYPE}$,  alpha, beta)
@:localvars(integer, n, k, lda, ldc)
@:defaults(trans='N', uplo='U', alpha=1, beta=0)
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call ${NAME}$(local_uplo, local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
#:enddef

module mfi_blas
use iso_fortran_env
implicit none

$:mfi_interface('?dotc',  COMPLEX_TYPES)
$:mfi_interface('?axpy',  DEFAULT_TYPES)
$:mfi_interface('i?amin', DEFAULT_TYPES)
$:mfi_interface('i?amax', DEFAULT_TYPES)
$:mfi_interface('?gemm',  DEFAULT_TYPES)
$:mfi_interface('?gemv',  DEFAULT_TYPES)
$:mfi_interface('?herk',  COMPLEX_TYPES)

contains

$:mfi_implement('?dotc',  COMPLEX_TYPES, dotc)
$:mfi_implement('?axpy',  DEFAULT_TYPES, axpy)
$:mfi_implement('i?amin', DEFAULT_TYPES, iamin_iamax)
$:mfi_implement('i?amax', DEFAULT_TYPES, iamin_iamax)
$:mfi_implement('?gemm',  DEFAULT_TYPES, gemm)
$:mfi_implement('?gemv',  DEFAULT_TYPES, gemv)
$:mfi_implement('?herk',  COMPLEX_TYPES, herk)

end module
