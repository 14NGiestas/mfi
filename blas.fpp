#:include "common.fpp"
module f77_blas
#:for name, types in BLAS_ROUTINES.items()
$:f77_interface(name, types)
#:endfor
end module

module mfi_blas
use f77_blas
use iso_fortran_env
implicit none

#:for name, types in BLAS_ROUTINES.items()
$:mfi_interface(name, types)
#:endfor

contains

#:block mfi_implement('?gemv')
pure subroutine mfi_SNAME(a, b, c, trans, alpha, beta, incx, incy)
@:args(DTYPE, in, a(:,:), x(:))
@:args(DTYPE, inout, y(:))
@:optargs(character, trans)
@:optargs(DTYPE,     alpha, beta)
@:optargs(integer,   incx,  incy)
@:localvars(integer, m, n, lda)
@:defaults(trans='N', alpha=1, beta=0, incx=1, incy=1)
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,incx,local_beta,y,incy)
end subroutine
#:endblock

#:block mfi_implement('?gemm')
pure subroutine mfi_SNAME(a, b, c, transa, transb, alpha, beta)
@:args(DTYPE, in, a(:,:), b(:,:))
@:args(DTYPE, inout, c(:,:))
@:optargs(character, transa, transb)
@:optargs(DTYPE,     alpha, beta)
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
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
#:endblock

end module
