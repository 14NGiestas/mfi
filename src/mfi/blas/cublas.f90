module mfi_blas_cublas
    use iso_c_binding
    implicit none

interface
    pure subroutine cuda_malloc(devPtr, size, stat) bind(c,name="mfi_cuda_malloc")
        import
        type(c_ptr), intent(out) :: devPtr
        integer(c_size_t), value, intent(in) :: size
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine cuda_free(devPtr, stat) bind(c,name="mfi_cuda_free")
        import
        type(c_ptr), value, intent(in) :: devPtr
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine cudaMemcpy(dst, src, count, kind) bind(c,name="cudaMemcpy")
        import
        type(c_ptr), value, intent(in) :: dst
        type(c_ptr), value, intent(in) :: src
        integer(c_size_t), value, intent(in) :: count
        integer(c_int), value, intent(in) :: kind
    end subroutine

    pure subroutine cublasCreate(handle, stat) bind(c,name="mfi_cublas_create")
        import
        type(c_ptr), intent(out) :: handle
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine cublasDestroy(handle, stat) bind(c,name="mfi_cublas_destroy")
        import
        type(c_ptr), value, intent(in) :: handle
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine cublasSetPointerMode(handle, mode, stat) bind(c,name="mfi_cublas_set_pointer_mode")
        import
        type(c_ptr), value, intent(in) :: handle
        integer(c_int), value, intent(in) :: mode
        integer(c_int), intent(out) :: stat
    end subroutine
end interface

!> cuBLAS operation constants
integer(c_int), parameter :: CUBLAS_OP_N = 0
integer(c_int), parameter :: CUBLAS_OP_T = 1
integer(c_int), parameter :: CUBLAS_OP_C = 2

!> cuBLAS fill mode constants (standard: UPPER=0, LOWER=1)
integer(c_int), parameter :: CUBLAS_FILL_MODE_UPPER = 0
integer(c_int), parameter :: CUBLAS_FILL_MODE_LOWER = 1

!> cuBLAS fill mode for TRSM/TRMM — empirically inverted vs standard
integer(c_int), parameter :: CUBLAS_TRSM_FILL_UPPER = 1
integer(c_int), parameter :: CUBLAS_TRSM_FILL_LOWER = 0

!> cuBLAS diagonal constants
integer(c_int), parameter :: CUBLAS_DIAG_NON_UNIT = 0
integer(c_int), parameter :: CUBLAS_DIAG_UNIT = 1

!> cuBLAS side constants
integer(c_int), parameter :: CUBLAS_SIDE_LEFT = 0
integer(c_int), parameter :: CUBLAS_SIDE_RIGHT = 1

!> CUDA memory copy direction constants
integer(c_int), parameter :: cudaMemcpyHostToDevice = 1
integer(c_int), parameter :: cudaMemcpyDeviceToHost = 2

!> cuBLAS pointer mode
integer(c_int), parameter :: CUBLAS_POINTER_MODE_HOST = 0
    type(c_ptr), save :: mfi_cublas_handle = c_null_ptr
end module

!> Extensions module
