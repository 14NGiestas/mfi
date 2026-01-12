#:mute
#! Create boilerplate to allocate memory in GPU device
#:def allocate(*varlist)
    #:for var in varlist
        type(c_ptr) :: device_${var}$
    #:endfor
    #:for var in varlist
        call cublas_alloc(size(${var}$),wp,device_${var}$)
    #:endfor
#:enddef

#! Create boilerplate to copy a matrix to GPU device
#:def set_matrix(*varlist)
    #:for var in varlist
        call cublas_set_matrix(size(${var}$,1),size(${var}$,2),wp,    &
                            ${var}$,       max(1,size(${var}$,1)), &
                            device_${var}$,max(1,size(${var}$,1)))
    #:endfor
#:enddef

#! Create boilerplate to get a matrix from GPU device
#:def get_matrix(*varlist)
    #:for var in varlist
        call cublas_get_matrix(size(${var}$,1),size(${var}$,2),wp,    &
                            device_${var}$,max(1,size(${var}$,1)), &
                            ${var}$,       max(1,size(${var}$,1)))
    #:endfor
#:enddef

#! Create boilerplate to copy a vector content to GPU device
#:def set_vector(*varlist)
    #:for var in varlist
        call cublas_set_vector(size(${var}$),wp,${var}$,1,device_${var}$,1)
    #:endfor
#:enddef

#! Create boilerplate to get a vector from GPU device
#:def get_vector(*varlist)
    #:for var in varlist
        call cublas_get_vector(size(${var}$),wp,device_${var}$,1,${var}$,1)
    #:endfor
#:enddef

#! Create boilerplate to free memory in GPU device
#:def deallocate(*varlist)
    #:for var in varlist
        call cublas_free(device_${var}$)
    #:endfor
#:enddef

#:def cublas_interfaces()

#:if defined('MFI_EXTENSIONS') and defined('MFI_USE_CUBLAS')
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
        type(c_ptr),    intent(in) :: x
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
#:endif
#:enddef
#:endmute
