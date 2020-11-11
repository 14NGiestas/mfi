#:include "common.fpp"

#:def axpy(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(x, y, a, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:))
@:args(${TYPE}$, inout, y(:))
@:optional(${TYPE}$, a)
@:optional(integer, incx, incy)
@:localvars(integer, n)
@:defaults(a=1, incx=1, incy=1)
    N = size(X)
    call ${F77_NAME}$(n,a_,x,incx_,y,incy_)
end subroutine
#:enddef

#:def copy_swap(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(x, y, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:))
@:args(${TYPE}$, inout, y(:))
@:optional(integer, incx, incy)
@:localvars(integer, n)
@:defaults(incx=1, incy=1)
    N = size(X)
    call ${F77_NAME}$(n,x,incx_,y,incy_)
end subroutine
#:enddef

#:def dot_product(MFI_NAME,F77_NAME,TYPE,KIND)
pure function ${MFI_NAME}$(x, y, incx, incy)
@:parameter(integer, wp=${KIND}$)
    ${TYPE}$ :: ${MFI_NAME}$
@:args(${TYPE}$, in, x(:), y(:))
@:localvars(integer, n)
@:optional(integer, incx, incy)
@:defaults(incx=1, incy=1)
    N = size(X)
    ${MFI_NAME}$ = ${F77_NAME}$(n,x,incx_,y,incy_)
end function
#:enddef

#:def rotm(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(x, y, param, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, x(:), y(:))
@:args(${TYPE}$, in, param(5))
@:optional(integer, incx, incy)
@:localvars(integer, n)
@:defaults(incx=1, incy=1)
    N = size(X)
    call ${F77_NAME}$(n,x,incx_,y,incy_,param)
end subroutine
#:enddef

! FIXME MFI and F77 interfaces are the same (this subroutine is not needed)
#:def rotmg(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(d1, d2, x1, y1, param)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    y1)
@:args(${TYPE}$, out,   param(5))
@:args(${TYPE}$, inout, d1, d2, x1)
    call ${F77_NAME}$(d1, d2, x1, y1, param)
end subroutine
#:enddef

#:def iamin_iamax(MFI_NAME,F77_NAME,TYPE,KIND)
pure function ${MFI_NAME}$(x, incx)
@:parameter(integer, wp=${KIND}$)
    integer :: ${MFI_NAME}$
@:args(${TYPE}$, in, x(:))
@:optional(integer, incx)
@:localvars(integer, n)
@:defaults(incx=1)
    n = size(x)
    ${MFI_NAME}$ = ${F77_NAME}$(n,x,incx_)
end function
#:enddef

#:def gbmv(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, x, y, kl, m, alpha, beta, trans, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in, a(:,:), x(:))
@:args(${TYPE}$, inout, y(:))
@:optional(character, trans)
@:optional(${TYPE}$,  alpha, beta)
@:optional(integer, kl, m, incx,  incy)
@:localvars(integer, n, ku, lda)
    n = size(a,2)
    lda = max(1,size(a,1))
@:defaults(kl=(lda-1)/2, m=n, trans='N', alpha=1, beta=0, incx=1, incy=1)
    ku = lda-kl_-1
    call ${F77_NAME}$(trans_,m_,n,kl_,ku,alpha_,a,lda,x,incx_,beta_,y,incy_)
end subroutine
#:enddef

#:def gemv(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, x, y, trans, alpha, beta, incx, incy)
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
    call ${F77_NAME}$(trans_,m,n,alpha_,a,lda,x,incx_,beta_,y,incy_)
end subroutine
#:enddef

#:def gemm(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, b, c, transa, transb, alpha, beta)
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
    if (transa_ == 'N' .or. transa_ == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    call ${F77_NAME}$(transa_,transb_,m,n,k,alpha_,a,lda,b,ldb,beta_,c,ldc)
end subroutine
#:enddef

#:def herk(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, trans, uplo)
@:optional(${TYPE}$,  alpha, beta)
@:localvars(integer, n, k, lda, ldc)
@:defaults(trans='N', uplo='U', alpha=1, beta=0)
    n = size(c,2)
    if (trans_ == 'N' .or. trans_ == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call ${F77_NAME}$(uplo_,trans_,n,k,alpha_,a,lda,beta_,c,ldc)
end subroutine
#:enddef

module mfi_blas
use iso_fortran_env
use f77_blas
implicit none

! BLAS level 1
!$:mfi_interface('?asum',  DEFAULT_TYPES)
$:mfi_interface('?axpy',  DEFAULT_TYPES)
$:mfi_interface('?copy',  DEFAULT_TYPES)
!$:mfi_interface('?dot',   REAL_TYPES)
!$:mfi_interface('sdsdot', REAL_TYPES)
$:mfi_interface('?dotu',  COMPLEX_TYPES)
$:mfi_interface('?dotc',  COMPLEX_TYPES)
!$:mfi_interface('?nrm2',  DEFAULT_TYPES)
!$:mfi_interface('?rot',   DEFAULT_TYPES)
!$:mfi_interface('?rotg',  DEFAULT_TYPES)
$:mfi_interface('?rotm',  REAL_TYPES)
$:mfi_interface('?rotmg', REAL_TYPES)
!$:f77_interface('?scal')
$:mfi_interface('?swap',  DEFAULT_TYPES)
$:mfi_interface('i?amin', DEFAULT_TYPES)
$:mfi_interface('i?amax', DEFAULT_TYPES)

! BLAS level 2
$:mfi_interface('?gbmv',  DEFAULT_TYPES)
$:mfi_interface('?gemv',  DEFAULT_TYPES)

! BLAS level 3
$:mfi_interface('?gemm',  DEFAULT_TYPES)
$:mfi_interface('?herk',  COMPLEX_TYPES)

contains

! BLAS level 1
!$:mfi_interface('?asum',  DEFAULT_TYPES)
$:mfi_implement('?axpy',  DEFAULT_TYPES, axpy)
$:mfi_implement('?copy',  DEFAULT_TYPES, copy_swap)
!$:mfi_interface('?dot',   REAL_TYPES)
!$:mfi_interface('sdsdot', REAL_TYPES)
$:mfi_implement('?dotu',  COMPLEX_TYPES, dot_product)
$:mfi_implement('?dotc',  COMPLEX_TYPES, dot_product)
!$:mfi_interface('?nrm2',  DEFAULT_TYPES)
!$:mfi_interface('?rot',   DEFAULT_TYPES)
!$:mfi_interface('?rotg',  DEFAULT_TYPES)
$:mfi_implement('?rotm',  REAL_TYPES, rotm)
$:mfi_implement('?rotmg', REAL_TYPES, rotmg)
!$:f77_interface('?scal')
$:mfi_implement('?swap',  DEFAULT_TYPES, copy_swap)
$:mfi_implement('i?amin', DEFAULT_TYPES, iamin_iamax)
$:mfi_implement('i?amax', DEFAULT_TYPES, iamin_iamax)

! BLAS level 2
$:mfi_implement('?gbmv',  DEFAULT_TYPES, gbmv)
$:mfi_implement('?gemv',  DEFAULT_TYPES, gemv)

! BLAS level 3
$:mfi_implement('?gemm',  DEFAULT_TYPES, gemm)
$:mfi_implement('?herk',  COMPLEX_TYPES, herk)

end module
