#:mute
#:include "../common.fpp"
#:include "cuda.fpp"

#:def mfi_gemm(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, b, c, transa, transb, alpha, beta)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, in,    a(:,:), b(:,:))
@:args(${TYPE}$, inout, c(:,:))
@:optional(character, in, transa, transb)
@:optional(${TYPE}$,  in, alpha, beta)
    integer :: m, n, k, lda, ldb, ldc
#:if defined('CUDA_SUPPORT')
    type(c_ptr) :: device_a, device_b, device_c
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
    call ${F77_NAME.replace('f77','gpu')}$(local_transa,local_transb,m,n,k,local_alpha,&
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

#:def gpu_gemm(NAME,TYPE,KIND)
pure subroutine ${NAME}$(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) &
    bind(C,name="${NAME}$")
    use iso_c_binding
    use iso_fortran_env
@:parameter(integer, wp=${KIND}$)
@:args({type(c_ptr), value}, in, a, b)
    type(c_ptr), value :: c ! void *
@:args({character(c_char),value}, in,  transa, transb)
@:args({${TYPE}$,value}, in,  alpha, beta)
@:args({integer(c_int),value}, in, m, n, k, lda, ldb, ldc)
end subroutine
#:enddef

#:endmute
module mfi_blas
use, intrinsic :: iso_c_binding
use iso_fortran_env

interface
    pure subroutine cublas_alloc(n, elemSize, devicePtr) &
        bind(c,name="cublasAlloc")
        use, intrinsic :: iso_c_binding
        use iso_fortran_env
        integer(c_int), intent(in), value :: n, elemSize
        type(c_ptr), target, intent(out) :: devicePtr ! void **devicePtr
    end subroutine

    pure subroutine cublas_free(devicePtr) &
        bind(c,name="cublasFree")
        use, intrinsic :: iso_c_binding
        use iso_fortran_env
        type(c_ptr), value, intent(in) :: devicePtr
    end subroutine

    pure subroutine cublas_set_vector(n, elemSize, x, incx, y, incy) &
        bind(c,name="cublasSetVector")
        use, intrinsic :: iso_c_binding
        use iso_fortran_env
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(in), value :: elemSize
        integer(c_int), intent(in), value :: incx, incy
        type(*),        intent(in) :: x(*)
        type(c_ptr),    value :: y
    end subroutine

    pure subroutine cublas_get_vector(n, elemSize, x, incx, y, incy) &
        bind(c,name="cublasGetVector")
        use, intrinsic :: iso_c_binding
        use iso_fortran_env
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(in), value :: elemSize
        integer(c_int), intent(in), value :: incx, incy
        type(c_ptr),    intent(in), value :: x
        type(*),        intent(inout) :: y(*)
    end subroutine

    pure subroutine cublas_set_matrix(rows, cols, elemSize, a, lda, b, ldb) &
        bind(c,name="cublasSetMatrix")
        use, intrinsic :: iso_c_binding
        use iso_fortran_env
        integer(c_int), intent(in), value :: rows, cols
        integer(c_int), intent(in), value :: elemSize
        integer(c_int), intent(in), value :: lda, ldb
        type(*),        intent(in)    :: a(lda,*)
        type(c_ptr),    value :: b
    end subroutine

    pure subroutine cublas_get_matrix(rows, cols, elemSize, a, lda, b, ldb) &
        bind(c,name="cublasGetMatrix")
        use, intrinsic :: iso_c_binding
        use iso_fortran_env
        integer(c_int), intent(in), value :: rows, cols
        integer(c_int), intent(in), value :: elemSize
        integer(c_int), intent(in), value :: lda, ldb
        type(c_ptr),    intent(in), value :: a
        type(*),        intent(inout) :: b(ldb,*)
    end subroutine

end interface

$:f77_interface('?gemm',  DEFAULT_TYPES, f77_gemm)
$:gpu_interface('?gemm',  DEFAULT_TYPES, gpu_gemm)
$:mfi_interface('?gemm',  DEFAULT_TYPES)

contains

$:mfi_implement('?gemm',  DEFAULT_TYPES, mfi_gemm)

end module
