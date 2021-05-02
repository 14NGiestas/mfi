#:mute
#:include "../common.fpp"
#:include "cuda.fpp"

#:def mfi_gemm(MFI_NAME,F77_NAME,TYPE,KIND)
#:if defined('CUDA_SUPPORT')
subroutine ${MFI_NAME}$(a, b, c, transa, transb, alpha, beta)
#:else
pure subroutine ${MFI_NAME}$(a, b, c, transa, transb, alpha, beta)
#:endif
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:), b(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, in, transa, transb)
@:optional(${TYPE}$,  in, alpha, beta)
    integer :: m, n, k, lda, ldb, ldc
#:if defined('CUDA_SUPPORT')
    integer(INT64) :: device_a, device_b, device_c
#:endif
@:defaults(transa='N', transb='N', alpha=1_wp, beta=0_wp)
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
#:if defined('CUDA_SUPPORT')
    @:copyin(a,b,c)
    call ${F77_NAME}$(local_transa,local_transb,m,n,k,local_alpha,&
        device_a,lda,device_b,ldb,local_beta,device_c,ldc)
    @:copyout(c)
    @:delete(a,b)
#:else
    call ${F77_NAME}$(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
#:endif
end subroutine
#:enddef

#:def f77_gemm(NAME,TYPE,KIND)
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
#:endmute

module mfi_blas
use iso_fortran_env
implicit none

$:f77_interface('?gemm',  DEFAULT_TYPES, f77_gemm)
$:mfi_interface('?gemm',  DEFAULT_TYPES)

contains

$:mfi_implement('?gemm',  DEFAULT_TYPES, mfi_gemm)

end module
