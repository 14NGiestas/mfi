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

#:block mfi_implement('?rotg')
#:endblock

#:block mfi_implement('?rotmg')
#:endblock

#:block mfi_implement('?rot')
#:endblock

#:block mfi_implement('?rotm')
#:endblock

#:block mfi_implement('?swap')
#:endblock

#:block mfi_implement('?scal')
#:endblock

#:block mfi_implement('?copy')
#:endblock

#:block mfi_implement('?axpy')
pure subroutine mfi_SNAME(X, Y, a, incx, incy)
@:args(DTYPE, in,    X(:))
@:args(DTYPE, inout, Y(:))
@:optargs(integer, a, incx, incy)
@:localvars(integer, N)
@:defaults(a=1, incx=1, incy=1)
N = size(X)
call f77_axpy(N,local_a,X,incx,Y,incy)
end subroutine
#:endblock

#:block mfi_implement('?dot')
#:endblock

#:block mfi_implement('?dotu')
#:endblock

#:block mfi_implement('?dotc')
pure function mfi_SNAME(X, Y, incx, incy)
@:args(DTYPE, in, X(:), Y(:))
@:localvars(integer, N, mfi_SNAME)
@:optargs(integer, incx, incy)
@:defaults(incx=1, incy=1)
N = size(X)
mfi_SNAME = f77_dotc(N,X,incx,Y,incy)
end function
#:endblock

#:block mfi_implement('?nrm2')
#:endblock

#:block mfi_implement('?asum')
#:endblock

#:block mfi_implement('i?amax')
pure function mfi_SNAME(X, incx)
@:args(DTYPE, in, X(:))
@:optargs(integer, incx)
@:localvars(integer, N, mfi_SNAME)
@:defaults(incx=1)
N = size(X)
mfi_SNAME = f77_iamin(N,X,incx)
end function
#:endblock

#:block mfi_implement('i?amin')
pure function mfi_SNAME(X, incx)
@:args(DTYPE, in, X(:))
@:optargs(integer, incx)
@:localvars(integer, N, mfi_SNAME)
@:defaults(incx=1)
N = size(X)
mfi_SNAME = f77_iamax(N,X,incx)
end function
#:endblock

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

#:block mfi_implement('?gbmv')
#:endblock

#:block mfi_implement('?hemv')
#:endblock

#:block mfi_implement('?hbmv')
#:endblock

#:block mfi_implement('?hpmv')
#:endblock

#:block mfi_implement('?symv')
#:endblock

#:block mfi_implement('?sbmv')
#:endblock

#:block mfi_implement('?spmv')
#:endblock

#:block mfi_implement('?trmv')
#:endblock

#:block mfi_implement('?tbmv')
#:endblock

#:block mfi_implement('?tpmv')
#:endblock

#:block mfi_implement('?trsv')
#:endblock

#:block mfi_implement('?tbsv')
#:endblock

#:block mfi_implement('?tpsv')
#:endblock

#:block mfi_implement('?ger')
#:endblock

#:block mfi_implement('?geru')
#:endblock

#:block mfi_implement('?gerc')
#:endblock

#:block mfi_implement('?her')
#:endblock

#:block mfi_implement('?hpr')
#:endblock

#:block mfi_implement('?her2')
#:endblock

#:block mfi_implement('?hpr2')
#:endblock

#:block mfi_implement('?syr')
#:endblock

#:block mfi_implement('?spr')
#:endblock

#:block mfi_implement('?syr2')
#:endblock

#:block mfi_implement('?spr2')
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

#:block mfi_implement('?symm')
#:endblock

#:block mfi_implement('?hemm')
#:endblock

#:block mfi_implement('?syrk')
#:endblock

#:block mfi_implement('?herk')
pure subroutine mfi_SNAME(a, c, uplo, trans, alpha, beta)
@:args(DTYPE, in,    a(:,:))
@:args(DTYPE, inout, c(:,:))
@:optargs(character, trans, uplo)
@:optargs(DTYPE,     alpha, beta)
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
    call f77_herk(local_uplo, local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
#:endblock

#:block mfi_implement('?syr2k')
#:endblock

#:block mfi_implement('?her2k')
#:endblock

#:block mfi_implement('?trmm')
#:endblock

#:block mfi_implement('?trsm')
#:endblock

#:block mfi_implement('?axpyi')
#:endblock

#:block mfi_implement('?doti')
#:endblock

#:block mfi_implement('?dotci')
#:endblock

#:block mfi_implement('?dotui')
#:endblock

#:block mfi_implement('?gthr')
#:endblock

#:block mfi_implement('?gthrz')
#:endblock

#:block mfi_implement('?sctr')
#:endblock

#:block mfi_implement('?roti')
#:endblock

end module
