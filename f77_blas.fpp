#:include "common.fpp"

#:def dotc(NAME,TYPE,KIND)
pure function ${NAME}$(n, x, incx, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
    ${TYPE}$ :: ${NAME}$
@:args(${TYPE}$, in, x(*), y(*))
@:args(integer,  in, n, incx, incy)
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

#:def iamin_iamax(NAME,TYPE,KIND)
pure function ${NAME}$(n, x, incx)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
    integer :: ${NAME}$
@:args(${TYPE}$, in, x(*))
@:args(integer,  in, n, incx)
end function
#:enddef

#:def gemv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in,    a(lda,*), x(*))
@:args(${TYPE}$,  inout, y(*))
@:args(character, in,    trans)
@:args(${TYPE}$,  in,    alpha, beta)
@:args(integer,   in,    incx,  incy)
@:args(integer,   in,    m, n, lda)
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

module f77_blas
use iso_fortran_env
implicit none

$:f77_interface('?dotc',  COMPLEX_TYPES, dotc)
$:f77_interface('?axpy',  DEFAULT_TYPES, axpy)
$:f77_interface('i?amin', DEFAULT_TYPES, iamin_iamax)
$:f77_interface('i?amax', DEFAULT_TYPES, iamin_iamax)
$:f77_interface('?gemm',  DEFAULT_TYPES, gemm)
$:f77_interface('?gemv',  DEFAULT_TYPES, gemv)
$:f77_interface('?herk',  COMPLEX_TYPES, herk)

end module

