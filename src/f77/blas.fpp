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

#:def iamin_iamax(NAME,TYPE,KIND)
pure function ${NAME}$(n, x, incx)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
    integer :: ${NAME}$
@:args(${TYPE}$, in, x(*))
@:args(integer,  in, n, incx)
end function
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

#:def herk(NAME,TYPE,KIND)
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
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    n, k, lda, ldb, ldc)
end subroutine
#:enddef

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
$:f77_interface('i?amin', DEFAULT_TYPES, iamin_iamax)
$:f77_interface('i?amax', DEFAULT_TYPES, iamin_iamax)

! BLAS level 2
$:f77_interface('?gbmv',  DEFAULT_TYPES, gbmv)
$:f77_interface('?gemv',  DEFAULT_TYPES, gemv)

! BLAS level 3
$:f77_interface('?gemm',  DEFAULT_TYPES, gemm)
$:f77_interface('?herk',  COMPLEX_TYPES, herk)
$:f77_interface('?her2k', COMPLEX_TYPES, her2k)

end module

