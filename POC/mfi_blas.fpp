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
    if (MFI_USE_CUDA == 1) then
        block
        @:allocate(a,b,c)
        @:set_matrix(a,b)
        call ${F77_NAME}$(local_transa,local_transb,m,n,k, &
                 local_alpha,device_a,lda,device_b,ldb,local_beta,device_c,ldc)
        @:get_matrix(c)
        @:deallocate(a,b,c)
        end block
    else
        call ${F77_NAME}$(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
    end if
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
implicit none

interface
    pure subroutine cublas_alloc(n, elemSize, devicePtr) &
        bind(c,name="cublasAlloc")
        use, intrinsic :: iso_c_binding
        integer(c_int), intent(in), value :: n, elemSize
        type(c_ptr), target, intent(out) :: devicePtr ! void **devicePtr
    end subroutine

    pure subroutine cublas_free(devicePtr) &
        bind(c,name="cublasFree")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: devicePtr ! void *devicePtr
    end subroutine

    pure subroutine cublas_set_vector(n, elemSize, x, incx, y, incy) &
        bind(c,name="cublasSetVector")
        use, intrinsic :: iso_c_binding
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(in), value :: elemSize
        integer(c_int), intent(in), value :: incx, incy
        type(*),        intent(in) :: x(*)
        type(c_ptr),    value :: y
    end subroutine

    pure subroutine cublas_get_vector(n, elemSize, x, incx, y, incy) &
        bind(c,name="cublasGetVector")
        use, intrinsic :: iso_c_binding
        integer(c_int), intent(in), value :: n
        integer(c_int), intent(in), value :: elemSize
        integer(c_int), intent(in), value :: incx, incy
        type(c_ptr),    intent(in), value :: x
        type(*),        intent(inout) :: y(*)
    end subroutine

    pure subroutine cublas_set_matrix(rows, cols, elemSize, a, lda, b, ldb) &
        bind(c,name="cublasSetMatrix")
        use, intrinsic :: iso_c_binding
        integer(c_int), intent(in), value :: rows, cols
        integer(c_int), intent(in), value :: elemSize
        integer(c_int), intent(in), value :: lda, ldb
        type(*),        intent(in)    :: a(lda,*)
        type(c_ptr),    value :: b
    end subroutine

    pure subroutine cublas_get_matrix(rows, cols, elemSize, a, lda, b, ldb) &
        bind(c,name="cublasGetMatrix")
        use, intrinsic :: iso_c_binding
        integer(c_int), intent(in), value :: rows, cols
        integer(c_int), intent(in), value :: elemSize
        integer(c_int), intent(in), value :: lda, ldb
        type(c_ptr),    intent(in), value :: a
        type(*),        intent(inout)     :: b(ldb,*)
    end subroutine
end interface

$:f77_interface('?gemm',  DEFAULT_TYPES, f77_gemm)
$:gpu_interface('?gemm',  DEFAULT_TYPES, gpu_gemm)
$:mfi_interface('?gemm',  DEFAULT_TYPES)

integer(INT32), private :: MFI_USE_CUDA = 0

contains

$:mfi_implement('?gemm',  DEFAULT_TYPES, mfi_gemm)

subroutine mfi_init()
    integer :: info
    character :: env_mfi_use_cuda
    call get_environment_variable("MFI_USE_CUDA",env_mfi_use_cuda)
    read(env_mfi_use_cuda,*,iostat=info) MFI_USE_CUDA
end subroutine

end module
