#:mute
#! Allocate memory using modern cudaMalloc
#:def allocate(*varlist)
    #:for var in varlist
        type(c_ptr) :: device_${var}$
        integer(c_int) :: stat_${var}$
    #:endfor
    #:for var in varlist
        ! total_bytes = number_of_elements * storage_size
        call cuda_malloc(device_${var}$, &
            int(size(${var}$), c_size_t) * int(storage_size(${var}$)/8, c_size_t), stat_${var}$)
    #:endfor
#:enddef

#! Copy matrix to GPU: storage_size(var)/8 converts bits to bytes
#:def set_matrix(*varlist)
    #:for var in varlist
        call cublas_set_matrix(int(size(${var}$,1), c_int), int(size(${var}$,2), c_int), &
                            int(storage_size(${var}$)/8, c_int), &
                            ${var}$, int(max(1,size(${var}$,1)), c_int), &
                            device_${var}$, int(max(1,size(${var}$,1)), c_int))
    #:endfor
#:enddef

#:def get_matrix(*varlist)
    #:for var in varlist
        call cublas_get_matrix(int(size(${var}$,1), c_int), int(size(${var}$,2), c_int), &
                            int(storage_size(${var}$)/8, c_int), &
                            device_${var}$, int(max(1,size(${var}$,1)), c_int), &
                            ${var}$, int(max(1,size(${var}$,1)), c_int))
    #:endfor
#:enddef

#:def set_vector(*varlist)
    #:for var in varlist
        call cublas_set_vector(int(size(${var}$), c_int), int(storage_size(${var}$)/8, c_int), &
                            ${var}$, 1_c_int, device_${var}$, 1_c_int)
    #:endfor
#:enddef

#:def get_vector(*varlist)
    #:for var in varlist
        call cublas_get_vector(int(size(${var}$), c_int), int(storage_size(${var}$)/8, c_int), &
                            device_${var}$, 1_c_int, ${var}$, 1_c_int)
    #:endfor
#:enddef

#:def deallocate(*varlist)
    #:for var in varlist
        call cuda_free(device_${var}$, stat_${var}$)
    #:endfor
#:enddef

#:def cublas_interfaces()
#:if defined('MFI_EXTENSIONS') and defined('MFI_USE_CUBLAS')
interface
    ! CUDA Runtime Memory Management (Replaces cublasAlloc/Free)

    ! Define as a SUBROUTINE to allow INTENT(OUT) in a PURE context
    pure subroutine cuda_malloc(ptr, bytes, stat) bind(c, name="cudaMalloc")
        import :: c_ptr, c_size_t, c_int
        type(c_ptr), intent(out) :: ptr
        integer(c_size_t), value :: bytes
        ! This captures the C function return value
        integer(c_int), intent(out) :: stat 
    end subroutine
    
    pure subroutine cuda_free(ptr, stat) bind(c, name="cudaFree")
        import :: c_ptr, c_int
        type(c_ptr), value, intent(in) :: ptr
        integer(c_int), intent(out)    :: stat
    end subroutine

    ! Utility functions use type(*) for generic data support
    pure subroutine cublas_set_matrix(rows, cols, elemSize, A, lda, B, ldb) bind(c, name="cublasSetMatrix")
        import :: c_int, c_ptr
        integer(c_int), value :: rows, cols, elemSize, lda, ldb
        type(*), intent(in) :: A(*)
        type(c_ptr), value :: B
    end subroutine

    pure subroutine cublas_get_matrix(rows, cols, elemSize, A, lda, B, ldb) bind(c, name="cublasGetMatrix")
        import :: c_int, c_ptr
        integer(c_int), value :: rows, cols, elemSize, lda, ldb
        type(c_ptr), value :: A
        type(*), intent(inout) :: B(*)
    end subroutine

    pure subroutine cublas_set_vector(n, elemSize, x, incx, y, incy) bind(c, name="cublasSetVector")
        import :: c_int, c_ptr
        integer(c_int), value :: n, elemSize, incx, incy
        type(*), intent(in) :: x
        type(c_ptr), value :: y
    end subroutine

    pure subroutine cublas_get_vector(n, elemSize, x, incx, y, incy) bind(c, name="cublasGetVector")
        import :: c_int, c_ptr
        integer(c_int), value :: n, elemSize, incx, incy
        type(c_ptr), value :: x
        type(*), intent(inout) :: y
    end subroutine
end interface
#:endif
#:enddef
#:endmute

