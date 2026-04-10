#:mute

#! Create boilerplate to declare variables for GPU allocation (must be at top of block)
#:def declare(*varlist)
    integer(c_int) :: cuda_allocation_status
    #:for var in varlist
        type(c_ptr) :: device_${var}$
    #:endfor
#:enddef

#! Create boilerplate to allocate memory in GPU device using CUDA Runtime API
#:def allocate(*varlist)
    #:for var in varlist
        cuda_allocation_status = cuda_malloc(device_${var}$, &
                              int(size(${var}$) * storage_size(${var}$)/8, c_size_t))
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to copy a matrix to GPU device using cudaMemcpy
#:def set_matrix(*varlist)
    #:for var in varlist
        call cudaMemcpy(device_${var}$, c_loc(${var}$), &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
    #:endfor
#:enddef

#! Create boilerplate to get a matrix from GPU device using cudaMemcpy
#:def get_matrix(*varlist)
    #:for var in varlist
        call cudaMemcpy(c_loc(${var}$), device_${var}$, &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
    #:endfor
#:enddef

#! Create boilerplate to copy a vector to GPU device using cudaMemcpy
#:def set_vector(*varlist)
    #:for var in varlist
        call cudaMemcpy(device_${var}$, c_loc(${var}$), &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
    #:endfor
#:enddef

#! Create boilerplate to get a vector from GPU device using cudaMemcpy
#:def get_vector(*varlist)
    #:for var in varlist
        call cudaMemcpy(c_loc(${var}$), device_${var}$, &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
    #:endfor
#:enddef

#! Create boilerplate to free memory in GPU device using cudaFree
#:def deallocate(*varlist)
    #:for var in varlist
        cuda_allocation_status = cuda_free(device_${var}$)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating ${var}$'
    #:endfor
#:enddef

#:def cublas_interfaces()
interface
    pure function cuda_malloc(devPtr, size) bind(c,name="cudaMalloc") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(out) :: devPtr
        integer(c_size_t), value, intent(in) :: size
        integer(c_int) :: stat
    end function

    pure function cuda_free(devPtr) bind(c,name="cudaFree") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: devPtr
        integer(c_int) :: stat
    end function

    pure subroutine cudaMemcpy(dst, src, count, kind) bind(c,name="cudaMemcpy")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: dst
        type(c_ptr), value, intent(in) :: src
        integer(c_size_t), value, intent(in) :: count
        integer(c_int), value, intent(in) :: kind
    end subroutine

    pure function cublasCreate(handle) bind(c,name="cublasCreate_v2") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(out) :: handle
        integer(c_int) :: stat
    end function

    pure function cublasDestroy(handle) bind(c,name="cublasDestroy_v2") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: handle
        integer(c_int) :: stat
    end function

    pure function cublasSetPointerMode(handle, mode) bind(c,name="cublasSetPointerMode_v2") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: handle
        integer(c_int), value, intent(in) :: mode
        integer(c_int) :: stat
    end function
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

interface
    pure subroutine mfi_cublas_error(stat, name)
        integer(c_int), value, intent(in) :: stat
        character(*), intent(in) :: name
    end subroutine
end interface
#:enddef
#:endmute
