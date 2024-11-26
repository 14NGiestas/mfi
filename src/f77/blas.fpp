#:mute
#:include "common.fpp"

#:def asum(NAME,TYPE,KIND)
pure function ${NAME}$(n, x, incx)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
    ${TYPE}$ :: ${NAME}$
@:args(${TYPE}$, in, x(*))
@:args(integer,  in, n, incx)
end function
#:enddef

#:def axpy(NAME,TYPE,KIND)
pure subroutine ${NAME}$(n, a, x, incx, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(*), a)
@:args(${TYPE}$, inout, y(*))
@:args(integer,  in,    n, incx, incy)
end subroutine
#:enddef

#:def copy_swap(NAME,TYPE,KIND)
pure subroutine ${NAME}$(n, x, incx, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    x(*))
@:args(${TYPE}$, inout, y(*))
@:args(integer,  in,    n, incx, incy)
end subroutine
#:enddef

#:def dot_product(NAME,TYPE,KIND)
pure function ${NAME}$(n, x, incx, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
    ${TYPE}$ :: ${NAME}$
@:args(${TYPE}$, in, x(*), y(*))
@:args(integer,  in, n, incx, incy)
end function
#:enddef

#:def rotm(NAME,TYPE,KIND)
pure subroutine ${NAME}$(n, x, incx, y, incy, param)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, x(*), y(*))
@:args(${TYPE}$, in, param(5))
@:args(integer,  in, n, incx, incy)
end subroutine
#:enddef

#:def rotmg(NAME,TYPE,KIND)
pure subroutine ${NAME}$(d1, d2, x1, y1, param)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    y1)
@:args(${TYPE}$, out,   param(5))
@:args(${TYPE}$, inout, d1, d2, x1)
end subroutine
#:enddef

#:def gbmv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*), x(*))
@:args(${TYPE}$,  inout, y(*))
@:args(character, in,    trans)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    m, n, kl, ku, lda, incx, incy)
end subroutine
#:enddef

#:def gemv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*), x(*))
@:args(${TYPE}$,  inout, y(*))
@:args(character, in,    trans)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    m, n, lda, incx, incy)
end subroutine
#:enddef

#:def ger_gerc_geru(NAME,TYPE,KIND)
pure subroutine ${NAME}$(m, n, alpha, x, incx, y, incy, a, lda)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    x(*), y(*))
@:args(${TYPE}$,  inout, a(lda,*))
@:args(${TYPE}$,  in,    alpha)
@:args(integer,   in,    m, n, lda, incx, incy)
end subroutine
#:enddef

#:def hbmv_sbmv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*), x(*))
@:args(${TYPE}$,  inout, y(*))
@:args(character, in,    uplo)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    n, k, lda, incx, incy)
end subroutine
#:enddef

#:def hemv_symv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*), x(*))
@:args(${TYPE}$,  inout, y(*))
@:args(character, in,    uplo)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    n, lda, incx, incy)
end subroutine
#:enddef

#:def her(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, a, lda)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    x(*))
@:args(${TYPE}$,  inout, a(lda,*))
@:args(character, in,    uplo)
@:args(real(wp),  in,    alpha)
@:args(integer,   in,    n, lda, incx)
end subroutine
#:enddef

#:def syr(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, a, lda)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    x(*))
@:args(${TYPE}$,  inout, a(lda,*))
@:args(character, in,    uplo)
@:args(${TYPE}$,  in,    alpha)
@:args(integer,   in,    n, lda, incx)
end subroutine
#:enddef

#:def her_syr2(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    x(*), y(*))
@:args(${TYPE}$,  inout, a(lda,*))
@:args(character, in,    uplo)
@:args(${TYPE}$,  in,    alpha)
@:args(integer,   in,    n, lda, incx, incy)
end subroutine
#:enddef

#:def hpmv_spmv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    ap(*), x(*))
@:args(${TYPE}$,  inout, y(*))
@:args(character, in,    uplo)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    n, incx, incy)
end subroutine
#:enddef

#:def hpr(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, ap)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    x(*))
@:args(${TYPE}$,  inout, ap(*))
@:args(character, in,    uplo)
@:args(real(wp),  in,    alpha)
@:args(integer,   in,    n, incx)
end subroutine
#:enddef

#:def spr(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, ap)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    x(*))
@:args(${TYPE}$,  inout, ap(*))
@:args(character, in,    uplo)
@:args(${TYPE}$,  in,    alpha)
@:args(integer,   in,    n, incx)
end subroutine
#:enddef

#:def hpr_spr2(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, alpha, x, incx, y, incy, ap)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    x(*), y(*))
@:args(${TYPE}$,  inout, ap(*))
@:args(character, in,    uplo)
@:args(${TYPE}$,  in,    alpha)
@:args(integer,   in,    n, incx, incy)
end subroutine
#:enddef

#:def tbmv_tbsv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*))
@:args(${TYPE}$,  inout, x(*))
@:args(character, in,    uplo, trans, diag)
@:args(integer,   in,    n, k, lda, incx)
end subroutine
#:enddef

#:def tpmv_tpsv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, trans, diag, n, ap, x, incx)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    ap(*))
@:args(${TYPE}$,  inout, x(*))
@:args(character, in,    uplo, trans, diag)
@:args(integer,   in,    n, incx)
end subroutine
#:enddef

#:def trmv_trsv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, trans, diag, n, a, lda, x, incx)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*))
@:args(${TYPE}$,  inout, x(*))
@:args(character, in,    uplo, trans, diag)
@:args(integer,   in,    n, lda, incx)
end subroutine
#:enddef

#:def gemm(NAME,TYPE,KIND)
pure subroutine ${NAME}$(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*), b(ldb,*))
@:args(${TYPE}$,  inout, c(ldc,*))
@:args(character, in,    transa, transb)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    m, n, k, lda, ldb, ldc)
end subroutine
#:enddef

#:def hemm_symm(NAME,TYPE,KIND)
pure subroutine ${NAME}$(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*), b(ldb,*))
@:args(${TYPE}$,  inout, c(ldc,*))
@:args(character, in,    side, uplo)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    m, n, lda, ldb, ldc)
end subroutine
#:enddef

#:def herk(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*))
@:args(${TYPE}$,  inout, c(ldc,*))
@:args(character, in,    trans, uplo)
@:args(real(wp),  in,    alpha, beta)
@:args(integer,   in,    n, k, lda, ldc)
end subroutine
#:enddef

#:def syrk(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*))
@:args(${TYPE}$,  inout, c(ldc,*))
@:args(character, in,    trans, uplo)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    n, k, lda, ldc)
end subroutine
#:enddef

#:def her2k(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*))
@:args(${TYPE}$,  in,    b(ldb,*))
@:args(${TYPE}$,  inout, c(ldc,*))
@:args(character, in,    trans, uplo)
@:args(${TYPE}$,  in,    alpha)
@:args(real(wp),  in,    beta)
@:args(integer,   in,    n, k, lda, ldb, ldc)
end subroutine
#:enddef

#:def syr2k(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*))
@:args(${TYPE}$,  in,    b(ldb,*))
@:args(${TYPE}$,  inout, c(ldc,*))
@:args(character, in,    trans, uplo)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    n, k, lda, ldb, ldc)
end subroutine
#:enddef

#:def trmm_trsm(NAME,TYPE,KIND)
pure subroutine ${NAME}$(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*))
@:args(${TYPE}$,  inout, b(ldb,*))
@:args(character, in,    side, uplo, transa, diag)
@:args(${TYPE}$,  in,    alpha)
@:args(integer,   in,    m, n, lda, ldb)
end subroutine
#:enddef

! BLAS Level 1 - Extensions
#:def iamax_iamin(NAME,TYPE,KIND)
pure function ${NAME}$(n, x, incx)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
    integer :: ${NAME}$
@:args(${TYPE}$, in, x(*))
@:args(integer,  in, n, incx)
end function
#:enddef

#:def iamin_stub(NAME,TYPE,KIND)
pure function ${NAME}$(n, x, incx)
@:parameter(integer, wp=${KIND}$)
    integer :: ${NAME}$
@:args(${TYPE}$, in, x(*))
@:args(integer,  in, n, incx)
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        ${NAME}$ = 0
        return
    end if
#:if TYPE is COMPLEX_TYPE
    ${NAME}$ = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
#:else
    ${NAME}$ = minloc(x(1:n:incx),dim=1)
#:endif
end function
#:enddef

#:endmute
module f77_blas
use iso_fortran_env
implicit none

!FIXME rot, dot, rotg, nrm2: problem with functions that have TYPE /= TYPE_result
!https://spec.oneapi.com/versions/latest/elements/oneMKL/source/domains/blas/asum.html#onemkl-blas-asum
!FIXME sdsdot: Weird specific interface: computes a vec-vect dot product but perform a sum
!FIXME scal: problem with functions that have TYPE /= TYPE_scalar
!https://spec.oneapi.com/versions/latest/elements/oneMKL/source/domains/blas/scal.html#onemkl-blas-scal

! BLAS level 1
!!$:f77_interface('?asum',  DEFAULT_TYPES, asum, result=REAL_TYPES)
$:f77_interface('?axpy',  DEFAULT_TYPES, axpy)
$:f77_interface('?copy',  DEFAULT_TYPES, copy_swap)
!$:f77_interface('?dot',  REAL_TYPES, dot_product, result=REAL_TYPES)
!$:f77_interface('sdsdot')
$:f77_interface('?dotu',  COMPLEX_TYPES, dot_product)
$:f77_interface('?dotc',  COMPLEX_TYPES, dot_product)
!$:f77_interface('?nrm2', DEFAULT_TYPES, nrm2, result=REAL_TYPES)
!$:f77_interface('?rot',  DEFAULT_TYPES, rot,  result=REAL_TYPES)
!$:f77_interface('?rotg', DEFAULT_TYPES, rotg, result=REAL_TYPES)
$:f77_interface('?rotm',  REAL_TYPES,    rotm)
$:f77_interface('?rotmg', REAL_TYPES,    rotmg)
!$:f77_interface('?scal')
$:f77_interface('?swap',  DEFAULT_TYPES, copy_swap)

! BLAS level 2
$:f77_interface('?gbmv',  DEFAULT_TYPES, gbmv)
$:f77_interface('?gemv',  DEFAULT_TYPES, gemv)
$:f77_interface('?ger',   REAL_TYPES,    ger_gerc_geru)
$:f77_interface('?gerc',  COMPLEX_TYPES, ger_gerc_geru)
$:f77_interface('?geru',  COMPLEX_TYPES, ger_gerc_geru)
$:f77_interface('?hbmv',  COMPLEX_TYPES, hbmv_sbmv)
$:f77_interface('?hemv',  COMPLEX_TYPES, hemv_symv)
$:f77_interface('?her',   COMPLEX_TYPES, her)
$:f77_interface('?her2',  COMPLEX_TYPES, her_syr2)
$:f77_interface('?hpmv',  COMPLEX_TYPES, hpmv_spmv)
$:f77_interface('?hpr',   COMPLEX_TYPES, hpr)
$:f77_interface('?hpr2',  COMPLEX_TYPES, hpr_spr2)
$:f77_interface('?sbmv',  REAL_TYPES,    hbmv_sbmv)
$:f77_interface('?spmv',  REAL_TYPES,    hpmv_spmv)
$:f77_interface('?spr',   REAL_TYPES,    spr)
$:f77_interface('?spr2',  REAL_TYPES,    hpr_spr2)
$:f77_interface('?symv',  REAL_TYPES,    hemv_symv)
$:f77_interface('?syr',   REAL_TYPES,    syr)
$:f77_interface('?syr2',  REAL_TYPES,    her_syr2)
$:f77_interface('?tbmv',  DEFAULT_TYPES, tbmv_tbsv)
$:f77_interface('?tbsv',  DEFAULT_TYPES, tbmv_tbsv)
$:f77_interface('?tpmv',  DEFAULT_TYPES, tpmv_tpsv)
$:f77_interface('?tpsv',  DEFAULT_TYPES, tpmv_tpsv)
$:f77_interface('?trmv',  DEFAULT_TYPES, trmv_trsv)
$:f77_interface('?trsv',  DEFAULT_TYPES, trmv_trsv)

! BLAS level 3
$:f77_interface('?gemm',  DEFAULT_TYPES, gemm)
$:f77_interface('?hemm',  COMPLEX_TYPES, hemm_symm)
$:f77_interface('?herk',  COMPLEX_TYPES, herk)
$:f77_interface('?her2k', COMPLEX_TYPES, her2k)
$:f77_interface('?symm',  REAL_TYPES,    hemm_symm)
$:f77_interface('?syrk',  REAL_TYPES,    syrk)
$:f77_interface('?syr2k', REAL_TYPES,    syr2k)
$:f77_interface('?trmm',  DEFAULT_TYPES, trmm_trsm)
$:f77_interface('?trsm',  DEFAULT_TYPES, trmm_trsm)


! Specific interfaces for slamch and dlamch 
! as fortran can't differentiate them by the return kind
interface
    pure real(REAL32) function slamch(cmach)
        import :: REAL32
        character, intent(in) :: cmach
    end function

    pure real(REAL64) function dlamch(cmach)
        import :: REAL64
        character, intent(in) :: cmach
    end function
end interface

! Extensions
! BLAS Level 1 - Utils / Extensions
$:f77_interface('i?amax', DEFAULT_TYPES, iamax_iamin)
#:if defined('MFI_EXTENSIONS')
  #:if defined('MFI_LINK_EXTERNAL')
! Link with a external source
$:f77_interface('i?amin', DEFAULT_TYPES, iamax_iamin)
  #:else
! Implement the blas extensions in
$:f77_interface_internal('i?amin', DEFAULT_TYPES)
contains
$:f77_implement('i?amin', DEFAULT_TYPES, iamin_stub)
  #:endif
#:endif

end module

