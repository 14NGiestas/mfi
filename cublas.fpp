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
        call cuda_malloc(device_${var}$, &
                              int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to copy a matrix to GPU device using cudaMemcpy
#:def set_matrix(*varlist)
    #:for var in varlist
        call cudaMemcpy(device_${var}$, c_loc(${var}$), &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        cudaMemcpyHostToDevice, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMemcpy (H->D) failed for ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to get a matrix from GPU device using cudaMemcpy
#:def get_matrix(*varlist)
    #:for var in varlist
        call cudaMemcpy(c_loc(${var}$), device_${var}$, &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        cudaMemcpyDeviceToHost, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMemcpy (D->H) failed for ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to copy a vector to GPU device using cudaMemcpy
#:def set_vector(*varlist)
    #:for var in varlist
        call cudaMemcpy(device_${var}$, c_loc(${var}$), &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        cudaMemcpyHostToDevice, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMemcpy (H->D) failed for ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to get a vector from GPU device using cudaMemcpy
#:def get_vector(*varlist)
    #:for var in varlist
        call cudaMemcpy(c_loc(${var}$), device_${var}$, &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        cudaMemcpyDeviceToHost, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMemcpy (D->H) failed for ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to free memory in GPU device using cudaFree
#:def deallocate(*varlist)
    #:for var in varlist
        call cuda_free(device_${var}$, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating ${var}$'
    #:endfor
#:enddef

#:def cublas_interfaces()
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

    pure subroutine cudaMemcpy(dst, src, count, kind, stat) bind(c,name="mfi_cuda_memcpy")
        import
        type(c_ptr), value, intent(in) :: dst
        type(c_ptr), value, intent(in) :: src
        integer(c_size_t), value, intent(in) :: count
        integer(c_int), value, intent(in) :: kind
        integer(c_int), intent(out) :: stat
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

    pure subroutine mfi_cublas_set_count(count) bind(c,name="mfi_cublas_set_count")
        import
        integer(c_int), value, intent(in) :: count
    end subroutine

    pure subroutine mfi_cublas_init_handles(stat) bind(c,name="mfi_cublas_init_handles")
        import
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine mfi_cublas_get_thread_handle(out_handle, stat) bind(c,name="mfi_cublas_get_thread_handle")
        import
        type(c_ptr), intent(out) :: out_handle
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine mfi_cublas_finalize_all(stat) bind(c,name="mfi_cublas_finalize_all")
        import
        integer(c_int), intent(out) :: stat
    end subroutine

    pure function mfi_cublas_read_env() bind(c,name="mfi_cublas_read_env")
        import
        integer(c_int) :: mfi_cublas_read_env
    end function

    pure function mfi_cublas_read_omp_threads() bind(c,name="mfi_cublas_read_omp_threads")
        import
        integer(c_int) :: mfi_cublas_read_omp_threads
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

!> cuBLAS status constants
integer(c_int), parameter :: CUBLAS_STATUS_SUCCESS = 0
integer(c_int), parameter :: CUBLAS_STATUS_NOT_INITIALIZED = 1
integer(c_int), parameter :: CUBLAS_STATUS_ALLOC_FAILED = 3
integer(c_int), parameter :: CUBLAS_STATUS_INVALID_VALUE = 7
integer(c_int), parameter :: CUBLAS_STATUS_ARCH_MISMATCH = 8
integer(c_int), parameter :: CUBLAS_STATUS_MAPPING_ERROR = 11
integer(c_int), parameter :: CUBLAS_STATUS_EXECUTION_FAILED = 13
integer(c_int), parameter :: CUBLAS_STATUS_INTERNAL_ERROR = 14
integer(c_int), parameter :: CUBLAS_STATUS_NOT_SUPPORTED = 15
integer(c_int), parameter :: CUBLAS_STATUS_LICENSE_ERROR = 16
#:enddef
#:endmute
