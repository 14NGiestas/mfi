#:mute

#! Native HIP and rocBLAS support macros. This file mirrors cublas.fpp.
#! Fypp builds one global macro table. A second macro named declare would
#! overwrite the cuBLAS one, because both files pass through fypp. So every
#! macro here carries the hip_ prefix.

#! Create boilerplate to declare variables for GPU allocation (must be at top of block)
#:def hip_declare(*varlist)
    integer(c_int) :: hip_allocation_status
    #:for var in varlist
        type(c_ptr) :: device_${var}$
    #:endfor
#:enddef

#! Create boilerplate to allocate memory in the GPU device using the HIP runtime
#:def hip_allocate(*varlist)
    #:for var in varlist
        call hip_malloc(device_${var}$, &
                              int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                              hip_allocation_status)
        if (hip_allocation_status /= 0) error stop 'hipMalloc failed allocating ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to copy a matrix to the GPU device
#:def hip_set_matrix(*varlist)
    #:for var in varlist
        call hip_memcpy(device_${var}$, c_loc(${var}$), &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        HIP_MEMCPY_HOST_TO_DEVICE, hip_allocation_status)
        if (hip_allocation_status /= 0) error stop 'hipMemcpy (H->D) failed for ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to get a matrix from the GPU device
#:def hip_get_matrix(*varlist)
    #:for var in varlist
        call hip_memcpy(c_loc(${var}$), device_${var}$, &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        HIP_MEMCPY_DEVICE_TO_HOST, hip_allocation_status)
        if (hip_allocation_status /= 0) error stop 'hipMemcpy (D->H) failed for ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to copy a vector to the GPU device
#:def hip_set_vector(*varlist)
    #:for var in varlist
        call hip_memcpy(device_${var}$, c_loc(${var}$), &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        HIP_MEMCPY_HOST_TO_DEVICE, hip_allocation_status)
        if (hip_allocation_status /= 0) error stop 'hipMemcpy (H->D) failed for ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to get a vector from the GPU device
#:def hip_get_vector(*varlist)
    #:for var in varlist
        call hip_memcpy(c_loc(${var}$), device_${var}$, &
                        int(size(${var}$) * storage_size(${var}$)/8, c_size_t), &
                        HIP_MEMCPY_DEVICE_TO_HOST, hip_allocation_status)
        if (hip_allocation_status /= 0) error stop 'hipMemcpy (D->H) failed for ${var}$'
    #:endfor
#:enddef

#! Create boilerplate to free memory in the GPU device
#:def hip_deallocate(*varlist)
    #:for var in varlist
        call hip_free(device_${var}$, hip_allocation_status)
        if (hip_allocation_status /= 0) error stop 'hipFree failed deallocating ${var}$'
    #:endfor
#:enddef

#:def hip_interfaces()
interface
    pure subroutine hip_malloc(devPtr, size, stat) bind(c,name="mfi_hip_malloc")
        import
        type(c_ptr), intent(out) :: devPtr
        integer(c_size_t), value, intent(in) :: size
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine hip_free(devPtr, stat) bind(c,name="mfi_hip_free")
        import
        type(c_ptr), value, intent(in) :: devPtr
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine hip_memcpy(dst, src, count, kind, stat) bind(c,name="mfi_hip_memcpy")
        import
        type(c_ptr), value, intent(in) :: dst
        type(c_ptr), value, intent(in) :: src
        integer(c_size_t), value, intent(in) :: count
        integer(c_int), value, intent(in) :: kind
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine rocblas_create(handle, stat) bind(c,name="mfi_hip_create")
        import
        type(c_ptr), intent(out) :: handle
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine rocblas_destroy(handle, stat) bind(c,name="mfi_hip_destroy")
        import
        type(c_ptr), value, intent(in) :: handle
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine rocblas_set_pointer_mode(handle, mode, stat) bind(c,name="mfi_hip_set_pointer_mode")
        import
        type(c_ptr), value, intent(in) :: handle
        integer(c_int), value, intent(in) :: mode
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine mfi_hip_get_thread_handle(out_handle, stat) bind(c,name="mfi_hip_get_thread_handle")
        import
        type(c_ptr), intent(out) :: out_handle
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine mfi_hip_finalize_all(stat) bind(c,name="mfi_hip_finalize_all")
        import
        integer(c_int), intent(out) :: stat
    end subroutine

    pure function mfi_hip_is_active_c() bind(c,name="mfi_hip_is_active")
        import
        integer(c_int) :: mfi_hip_is_active_c
    end function

    pure subroutine mfi_hip_force_gpu_c() bind(c,name="mfi_hip_force_gpu")
        import
    end subroutine

    pure subroutine mfi_hip_force_cpu_c() bind(c,name="mfi_hip_force_cpu")
        import
    end subroutine

    pure subroutine mfi_hip_lazy_init_c() bind(c,name="mfi_hip_lazy_init")
        import
    end subroutine

    pure subroutine mfi_hip_handle_get_c(handle, stat) bind(c,name="mfi_hip_get_thread_handle")
        import
        type(c_ptr), intent(out) :: handle
        integer(c_int), intent(out) :: stat
    end subroutine

    pure subroutine mfi_hip_set_threads_c(n) bind(c,name="mfi_hip_set_threads")
        import
        integer(c_int), value, intent(in) :: n
    end subroutine

    pure function mfi_hip_handle_count_c() bind(c,name="mfi_hip_handle_count")
        import
        integer(c_int) :: mfi_hip_handle_count_c
    end function
end interface

!> rocBLAS operation constants
integer(c_int), parameter :: ROCBLAS_OP_N = 111
integer(c_int), parameter :: ROCBLAS_OP_T = 112
integer(c_int), parameter :: ROCBLAS_OP_C = 113

!> rocBLAS fill mode constants (rocBLAS uses the standard order)
integer(c_int), parameter :: ROCBLAS_FILL_MODE_UPPER = 121
integer(c_int), parameter :: ROCBLAS_FILL_MODE_LOWER = 122

!> rocBLAS fill mode for TRSM and TRMM
integer(c_int), parameter :: ROCBLAS_TRSM_FILL_UPPER = 121
integer(c_int), parameter :: ROCBLAS_TRSM_FILL_LOWER = 122

!> rocBLAS diagonal constants
integer(c_int), parameter :: ROCBLAS_DIAG_NON_UNIT = 131
integer(c_int), parameter :: ROCBLAS_DIAG_UNIT = 132

!> rocBLAS side constants
integer(c_int), parameter :: ROCBLAS_SIDE_LEFT = 141
integer(c_int), parameter :: ROCBLAS_SIDE_RIGHT = 142

!> HIP memory copy direction constants
integer(c_int), parameter :: HIP_MEMCPY_HOST_TO_DEVICE = 1
integer(c_int), parameter :: HIP_MEMCPY_DEVICE_TO_HOST = 2

!> rocBLAS pointer mode
integer(c_int), parameter :: ROCBLAS_POINTER_MODE_HOST = 0

!> rocBLAS status constants
integer(c_int), parameter :: ROCBLAS_STATUS_SUCCESS = 0
integer(c_int), parameter :: ROCBLAS_STATUS_INVALID_HANDLE = 1
integer(c_int), parameter :: ROCBLAS_STATUS_NOT_IMPLEMENTED = 2
integer(c_int), parameter :: ROCBLAS_STATUS_INVALID_POINTER = 3
integer(c_int), parameter :: ROCBLAS_STATUS_INVALID_SIZE = 4
integer(c_int), parameter :: ROCBLAS_STATUS_MEMORY_ERROR = 5
integer(c_int), parameter :: ROCBLAS_STATUS_INTERNAL_ERROR = 6
integer(c_int), parameter :: ROCBLAS_STATUS_PERF_DEGRADED = 7
integer(c_int), parameter :: ROCBLAS_STATUS_SIZE_QUERY_MISMATCH = 8
integer(c_int), parameter :: ROCBLAS_STATUS_SIZE_INCREASED = 9
integer(c_int), parameter :: ROCBLAS_STATUS_SIZE_UNCHANGED = 10
integer(c_int), parameter :: ROCBLAS_STATUS_INVALID_VALUE = 11
integer(c_int), parameter :: ROCBLAS_STATUS_CONTINUE = 12
integer(c_int), parameter :: ROCBLAS_STATUS_CHECK_NUMERICS_FAIL = 13
integer(c_int), parameter :: ROCBLAS_STATUS_EXCLUDED_FROM_BUILD = 14
integer(c_int), parameter :: ROCBLAS_STATUS_ARCH_MISMATCH = 15
#:enddef
#:endmute
