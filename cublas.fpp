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
    ! CUDA Runtime API - Memory management
    function cuda_malloc(devPtr, size) bind(c,name="cudaMalloc") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), intent(out) :: devPtr
        integer(c_size_t), value, intent(in) :: size
        integer(c_int) :: stat
    end function

    function cuda_free(devPtr) bind(c,name="cudaFree") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: devPtr
        integer(c_int) :: stat
    end function

    subroutine cudaMemcpy(dst, src, count, kind) bind(c,name="cudaMemcpy")
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: dst
        type(c_ptr), value, intent(in) :: src
        integer(c_size_t), value, intent(in) :: count
        integer(c_int), value, intent(in) :: kind
    end subroutine

    ! cuBLAS v2 API - Handle management
    function cublasCreate(handle) bind(c,name="cublasCreate_v2") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), intent(out) :: handle
        integer(c_int) :: stat
    end function

    function cublasDestroy(handle) bind(c,name="cublasDestroy_v2") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: handle
        integer(c_int) :: stat
    end function

    function cublasSetPointerMode(handle, mode) bind(c,name="cublasSetPointerMode_v2") result(stat)
        use, intrinsic :: iso_c_binding
        type(c_ptr), value, intent(in) :: handle
        integer(c_int), value, intent(in) :: mode
        integer(c_int) :: stat
    end function
end interface
#:enddef
#:endmute
