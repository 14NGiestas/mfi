#:mute
#:include "common.fpp"

#:def lamch(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure function ${MFI_NAME}$(cmach, kind) result(res)
@:parameter(integer, wp=${KIND}$)
@:args(character, in, cmach)
@:args(${TYPE}$, in, kind)
    !! Just a kind placeholder
    ${TYPE}$ :: res
    res = ${ORIGINAL_NAME}$(cmach)
end function
#:enddef

#:def axpy(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(x, y, a, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:))
@:args(${TYPE}$, inout, y(:))
@:optional(${TYPE}$, in, a)
@:optional(integer, in, incx, incy)
    integer :: n
@:defaults(a=1.0_wp, incx=1, incy=1)
    N = size(X)
    call ${F77_NAME}$(n,local_a,x,local_incx,y,local_incy)
end subroutine
#:enddef

#:def copy_swap(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(x, y, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:))
@:args(${TYPE}$, inout, y(:))
@:optional(integer, in, incx, incy)
    integer :: n
@:defaults(incx=1, incy=1)
    N = size(X)
    call ${F77_NAME}$(n,x,local_incx,y,local_incy)
end subroutine
#:enddef

#:def dot_product(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure function ${MFI_NAME}$(x, y, incx, incy)
@:parameter(integer, wp=${KIND}$)
    ${TYPE}$ :: ${MFI_NAME}$
@:args(${TYPE}$, in, x(:), y(:))
    integer :: n
@:optional(integer, in, incx, incy)
@:defaults(incx=1, incy=1)
    N = size(X)
    ${MFI_NAME}$ = ${F77_NAME}$(n,x,local_incx,y,local_incy)
end function
#:enddef

#:def rotm(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(x, y, param, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, x(:), y(:))
@:args(${TYPE}$, in, param(5))
@:optional(integer, in, incx, incy)
    integer :: n
@:defaults(incx=1, incy=1)
    N = size(X)
    call ${F77_NAME}$(n,x,local_incx,y,local_incy,param)
end subroutine
#:enddef

#:def iamin_iamax(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure function ${MFI_NAME}$(x, incx)
@:parameter(integer, wp=${KIND}$)
    integer :: ${MFI_NAME}$
@:args(${TYPE}$, in, x(:))
@:optional(integer, in, incx)
    integer :: n
@:defaults(incx=1)
    n = size(x)
    ${MFI_NAME}$ = ${F77_NAME}$(n,x,local_incx)
end function
#:enddef

#:def gbmv(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, y, kl, m, alpha, beta, trans, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:), x(:))
@:args(${TYPE}$, inout, y(:))
@:optional(character, in, trans)
@:optional(${TYPE}$,  in, alpha, beta)
@:optional(integer,   in, kl, m, incx,  incy)
    integer :: n, ku, lda
    n = size(a,2)
    lda = max(1,size(a,1))
@:defaults(kl=(lda-1)/2, m=n, trans='N', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    ku = lda-local_kl-1
    call ${F77_NAME}$(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef

#:def gemv(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, y, trans, alpha, beta, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:), x(:))
@:args(${TYPE}$, inout, y(:))
@:optional(character, in, trans)
@:optional(${TYPE}$,  in, alpha, beta)
@:optional(integer,   in, incx,  incy)
    integer :: m, n, lda
@:defaults(trans='N', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call ${F77_NAME}$(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef

#:def ger_gerc_geru(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, y, alpha, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:), y(:))
@:args(${TYPE}$, inout, a(:,:))
@:optional(${TYPE}$,  in, alpha)
@:optional(integer,   in, incx, incy)
    integer :: m, n, lda
@:defaults(alpha=1.0_wp, incx=1, incy=1)
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call ${F77_NAME}$(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
#:enddef

#:def hbmv_sbmv(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, y, uplo, alpha, beta, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in, x(:), a(:,:))
@:args(${TYPE}$, inout, y(:))
@:optional(character, in, uplo)
@:optional(${TYPE}$,  in, alpha, beta)
@:optional(integer,   in, incx, incy)
    integer :: n, k, lda
@:defaults(uplo='U', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef

#:def hemv_symv(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, y, uplo, alpha, beta, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in, x(:), a(:,:))
@:args(${TYPE}$, inout, y(:))
@:optional(character, in, uplo)
@:optional(${TYPE}$,  in, alpha, beta)
@:optional(integer,   in, incx, incy)
    integer :: n, lda
@:defaults(uplo='U', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef

#:def her(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, uplo, alpha, incx)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:))
@:args(${TYPE}$, inout, a(:,:))
@:optional(character, in, uplo)
@:optional(real(wp),  in, alpha)
@:optional(integer,   in, incx)
    integer :: n, lda
@:defaults(uplo='U', alpha=1.0_wp, incx=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
#:enddef

#:def syr(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, uplo, alpha, incx)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:))
@:args(${TYPE}$, inout, a(:,:))
@:optional(character, in, uplo)
@:optional(${TYPE}$,  in, alpha)
@:optional(integer,   in, incx)
    integer :: n, lda
@:defaults(uplo='U', alpha=1.0_wp, incx=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
#:enddef

#:def her_syr2(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, y, uplo, alpha, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:), y(:))
@:args(${TYPE}$, inout, a(:,:))
@:optional(character, in, uplo)
@:optional(${TYPE}$,  in, alpha)
@:optional(integer,   in, incx, incy)
    integer :: n, lda
@:defaults(uplo='U', alpha=1.0_wp, incx=1, incy=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
#:enddef

#:def hpmv_spmv(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(ap, x, y, uplo, alpha, beta, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:), ap(:))
@:args(${TYPE}$, inout, y(:))
@:optional(character, in, uplo)
@:optional(${TYPE}$,  in, alpha, beta)
@:optional(integer,   in, incx, incy)
    integer :: n
@:defaults(uplo='U', alpha=1.0_wp, beta=0.0_wp, incx=1, incy=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
#:enddef

#:def hpr(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(ap, x, uplo, alpha, incx)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:))
@:args(${TYPE}$, inout, ap(:))
@:optional(character, in, uplo)
@:optional(real(wp),  in, alpha)
@:optional(integer,   in, incx)
    integer :: n
@:defaults(uplo='U', alpha=1.0_wp, incx=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
#:enddef

#:def spr(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(ap, x, uplo, alpha, incx)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:))
@:args(${TYPE}$, inout, ap(:))
@:optional(character, in, uplo)
@:optional(${TYPE}$,  in, alpha)
@:optional(integer,   in, incx)
    integer :: n
@:defaults(uplo='U', alpha=1.0_wp, incx=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
#:enddef

#:def hpr_spr2(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(ap, x, y, uplo, alpha, incx, incy)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(:), y(:))
@:args(${TYPE}$, inout, ap(:))
@:optional(character, in, uplo)
@:optional(${TYPE}$,  in, alpha)
@:optional(integer,   in, incx, incy)
    integer :: n
@:defaults(uplo='U', alpha=1.0_wp, incx=1, incy=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
#:enddef

#:def tbmv_tbsv(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, uplo, trans, diag, incx)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in, a(:,:))
@:args(${TYPE}$, inout, x(:))
@:optional(character, in, uplo, trans, diag)
@:optional(integer,   in, incx)
    integer :: n, k, lda
@:defaults(uplo='U', trans='N', diag='N', incx=1)
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
#:enddef

#:def tpmv_tpsv(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(ap, x, uplo, trans, diag, incx)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    ap(:))
@:args(${TYPE}$, inout, x(:))
@:optional(character, in, uplo, trans, diag)
@:optional(integer,   in, incx)
    integer :: n
@:defaults(uplo='U', trans='N', diag='N', incx=1)
    n = size(x)
    call ${F77_NAME}$(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
#:enddef

#:def trmv_trsv(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, x, uplo, trans, diag, incx)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in, a(:,:))
@:args(${TYPE}$, inout, x(:))
@:optional(character, in, uplo, trans, diag)
@:optional(integer,   in, incx)
    integer :: n, lda
@:defaults(uplo='U', trans='N', diag='N', incx=1)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
#:enddef

#:def gemm(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, b, c, transa, transb, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:), b(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, in, transa, transb)
@:optional(${TYPE}$,  in, alpha, beta)
    integer :: m, n, k, lda, ldb, ldc
@:defaults(transa='N', transb='N', alpha=1.0_wp, beta=0.0_wp)
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
    call ${F77_NAME}$(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:enddef

#:def hemm_symm(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, b, c, side, uplo, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:), b(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, in, side,  uplo)
@:optional(${TYPE}$,  in, alpha, beta)
    integer :: m, n, lda, ldb, ldc
@:defaults(side='L', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call ${F77_NAME}$(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:enddef

#:def herk(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, in, trans, uplo)
@:optional(real(wp),  in, alpha, beta)
    integer :: n, k, lda, ldc
@:defaults(trans='N', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call ${F77_NAME}$(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
#:enddef

#:def syrk(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, in, trans, uplo)
@:optional(${TYPE}$,  in, alpha, beta)
    integer :: n, k, lda, ldc
@:defaults(trans='N', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call ${F77_NAME}$(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
#:enddef

#:def her2k(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, b, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:))
@:args(${TYPE}$, in,    b(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, in, trans, uplo)
@:optional(${TYPE}$,  in, alpha)
@:optional(real(wp),  in, beta)
    integer :: n, k, lda, ldb, ldc
@:defaults(trans='N', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call ${F77_NAME}$(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:enddef

#:def syr2k(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, b, c, uplo, trans, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:))
@:args(${TYPE}$, in,    b(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, in, trans, uplo)
@:optional(${TYPE}$,  in, alpha, beta)
    integer :: n, k, lda, ldb, ldc
@:defaults(trans='N', uplo='U', alpha=1.0_wp, beta=0.0_wp)
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call ${F77_NAME}$(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:enddef

#:def trmm_trsm(MFI_NAME,F77_NAME,TYPE,KIND,ORIGINAL_NAME)
pure subroutine ${MFI_NAME}$(a, b, side, uplo, transa, diag, alpha)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:))
@:args(${TYPE}$, inout, b(:,:))
@:optional(character, in, side, uplo, transa, diag)
@:optional(${TYPE}$,  in, alpha)
    integer :: m, n, lda, ldb
@:defaults(side='L', uplo='U', transa='N', diag='N', alpha=1.0_wp)
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call ${F77_NAME}$(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
#:enddef
#:endmute
module mfi_blas
use iso_fortran_env
use f77_blas
use f77_blas, only: mfi_rotmg => f77_rotmg
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
!$:mfi_interface('?rotmg', REAL_TYPES)
!$:f77_interface('?scal')
$:mfi_interface('?swap',  DEFAULT_TYPES)

! BLAS level 2
$:mfi_interface('?gbmv',  DEFAULT_TYPES)
$:mfi_interface('?gemv',  DEFAULT_TYPES)
$:mfi_interface('?ger',   REAL_TYPES)
$:mfi_interface('?gerc',  COMPLEX_TYPES)
$:mfi_interface('?geru',  COMPLEX_TYPES)
$:mfi_interface('?hbmv',  COMPLEX_TYPES)
$:mfi_interface('?hemv',  COMPLEX_TYPES)
$:mfi_interface('?her',   COMPLEX_TYPES)
$:mfi_interface('?her2',  COMPLEX_TYPES)
$:mfi_interface('?hpmv',  COMPLEX_TYPES)
$:mfi_interface('?hpr',   COMPLEX_TYPES)
$:mfi_interface('?hpr2',  COMPLEX_TYPES)
$:mfi_interface('?sbmv',  REAL_TYPES)
$:mfi_interface('?spmv',  REAL_TYPES)
$:mfi_interface('?spr',   REAL_TYPES)
$:mfi_interface('?spr2',  REAL_TYPES)
$:mfi_interface('?symv',  REAL_TYPES)
$:mfi_interface('?syr',   REAL_TYPES)
$:mfi_interface('?syr2',  REAL_TYPES)
$:mfi_interface('?tbmv',  DEFAULT_TYPES)
$:mfi_interface('?tbsv',  DEFAULT_TYPES)
$:mfi_interface('?tpmv',  DEFAULT_TYPES)
$:mfi_interface('?tpsv',  DEFAULT_TYPES)
$:mfi_interface('?trmv',  DEFAULT_TYPES)
$:mfi_interface('?trsv',  DEFAULT_TYPES)

! BLAS level 3
$:mfi_interface('?gemm',  DEFAULT_TYPES)
$:mfi_interface('?hemm',  COMPLEX_TYPES)
$:mfi_interface('?herk',  COMPLEX_TYPES)
$:mfi_interface('?her2k', COMPLEX_TYPES)
$:mfi_interface('?symm',  REAL_TYPES)
$:mfi_interface('?syrk',  REAL_TYPES)
$:mfi_interface('?syr2k', REAL_TYPES)
$:mfi_interface('?trmm',  DEFAULT_TYPES)
$:mfi_interface('?trsm',  DEFAULT_TYPES)

! Extensions
! BLAS level 1 - Utils / Extensions
$:mfi_interface('i?amax', DEFAULT_TYPES)
#:if defined('MFI_EXTENSIONS')
$:mfi_interface('i?amin', DEFAULT_TYPES)
#:endif
$:mfi_interface('?lamch', REAL_TYPES)

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
!$:mfi_implement('?rotmg', REAL_TYPES, rotmg)
!$:f77_interface('?scal')
$:mfi_implement('?swap',  DEFAULT_TYPES, copy_swap)

! BLAS level 2
$:mfi_implement('?gbmv',  DEFAULT_TYPES, gbmv)
$:mfi_implement('?gemv',  DEFAULT_TYPES, gemv)
$:mfi_implement('?ger',   REAL_TYPES,    ger_gerc_geru)
$:mfi_implement('?gerc',  COMPLEX_TYPES, ger_gerc_geru)
$:mfi_implement('?geru',  COMPLEX_TYPES, ger_gerc_geru)
$:mfi_implement('?hbmv',  COMPLEX_TYPES, hbmv_sbmv)
$:mfi_implement('?hemv',  COMPLEX_TYPES, hemv_symv)
$:mfi_implement('?her',   COMPLEX_TYPES, her)
$:mfi_implement('?her2',  COMPLEX_TYPES, her_syr2)
$:mfi_implement('?hpmv',  COMPLEX_TYPES, hpmv_spmv)
$:mfi_implement('?hpr',   COMPLEX_TYPES, hpr)
$:mfi_implement('?hpr2',  COMPLEX_TYPES, hpr_spr2)
$:mfi_implement('?sbmv',  REAL_TYPES,    hbmv_sbmv)
$:mfi_implement('?spmv',  REAL_TYPES,    hpmv_spmv)
$:mfi_implement('?spr',   REAL_TYPES,    spr)
$:mfi_implement('?spr2',  REAL_TYPES,    hpr_spr2)
$:mfi_implement('?symv',  REAL_TYPES,    hemv_symv)
$:mfi_implement('?syr',   REAL_TYPES,    syr)
$:mfi_implement('?syr2',  REAL_TYPES,    her_syr2)
$:mfi_implement('?tbmv',  DEFAULT_TYPES, tbmv_tbsv)
$:mfi_implement('?tbsv',  DEFAULT_TYPES, tbmv_tbsv)
$:mfi_implement('?tpmv',  DEFAULT_TYPES, tpmv_tpsv)
$:mfi_implement('?tpsv',  DEFAULT_TYPES, tpmv_tpsv)
$:mfi_implement('?trmv',  DEFAULT_TYPES, trmv_trsv)
$:mfi_implement('?trsv',  DEFAULT_TYPES, trmv_trsv)

! BLAS level 3
$:mfi_implement('?gemm',  DEFAULT_TYPES, gemm)
$:mfi_implement('?hemm',  COMPLEX_TYPES, hemm_symm)
$:mfi_implement('?herk',  COMPLEX_TYPES, herk)
$:mfi_implement('?her2k', COMPLEX_TYPES, her2k)
$:mfi_implement('?symm',  REAL_TYPES,    hemm_symm)
$:mfi_implement('?syrk',  REAL_TYPES,    syrk)
$:mfi_implement('?syr2k', REAL_TYPES,    syr2k)
$:mfi_implement('?trmm',  DEFAULT_TYPES, trmm_trsm)
$:mfi_implement('?trsm',  DEFAULT_TYPES, trmm_trsm)

! Extensions
! BLAS level 1 - Utils / Extensions
$:mfi_implement('i?amax', DEFAULT_TYPES, iamin_iamax)
#:if defined('MFI_EXTENSIONS')
$:mfi_implement('i?amin', DEFAULT_TYPES, iamin_iamax)
#:endif
$:mfi_implement('?lamch', REAL_TYPES, lamch)

end module
