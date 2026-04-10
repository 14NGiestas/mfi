module f77_blas_cublas
use iso_c_binding
implicit none
interface
!> sgemm GPU (cublas v2) version
pure function cublasSgemm(handle, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) result(stat) &
    bind(C,name="cublasSgemm_v2")
    import
    integer, parameter :: wp = c_float
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: transb
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: k
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    integer(c_int),value, intent(in) :: ldc
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    type(c_ptr), value, intent(in) :: beta
    type(c_ptr), value, intent(in) :: c
    integer(c_int) :: stat
end function
!> sgemv GPU (cublas v2) version
pure function cublasSgemv(handle, trans, m, n, alpha, a, lda, x, incx, beta, y, incy) result(stat) &
    bind(C,name="cublasSgemv_v2")
    import
    integer, parameter :: wp = c_float
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: trans
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: incx
    integer(c_int),value, intent(in) :: incy
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: x
    type(c_ptr), value, intent(in) :: beta
    type(c_ptr), value, intent(in) :: y
    integer(c_int) :: stat
end function
!> strmm GPU (cublas v2) version
pure function cublasStrmm(handle, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) result(stat) &
    bind(C,name="cublasStrmm_v2")
    import
    integer, parameter :: wp = c_float
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: side
    integer(c_int),value, intent(in) :: uplo
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: diag
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int) :: stat
end function
!> strsm GPU (cublas v2) version
pure function cublasStrsm(handle, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) result(stat) &
    bind(C,name="cublasStrsm_v2")
    import
    integer, parameter :: wp = c_float
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: side
    integer(c_int),value, intent(in) :: uplo
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: diag
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int) :: stat
end function
!> dgemm GPU (cublas v2) version
pure function cublasDgemm(handle, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) result(stat) &
    bind(C,name="cublasDgemm_v2")
    import
    integer, parameter :: wp = c_double
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: transb
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: k
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    integer(c_int),value, intent(in) :: ldc
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    type(c_ptr), value, intent(in) :: beta
    type(c_ptr), value, intent(in) :: c
    integer(c_int) :: stat
end function
!> dgemv GPU (cublas v2) version
pure function cublasDgemv(handle, trans, m, n, alpha, a, lda, x, incx, beta, y, incy) result(stat) &
    bind(C,name="cublasDgemv_v2")
    import
    integer, parameter :: wp = c_double
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: trans
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: incx
    integer(c_int),value, intent(in) :: incy
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: x
    type(c_ptr), value, intent(in) :: beta
    type(c_ptr), value, intent(in) :: y
    integer(c_int) :: stat
end function
!> dtrmm GPU (cublas v2) version
pure function cublasDtrmm(handle, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) result(stat) &
    bind(C,name="cublasDtrmm_v2")
    import
    integer, parameter :: wp = c_double
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: side
    integer(c_int),value, intent(in) :: uplo
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: diag
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int) :: stat
end function
!> dtrsm GPU (cublas v2) version
pure function cublasDtrsm(handle, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) result(stat) &
    bind(C,name="cublasDtrsm_v2")
    import
    integer, parameter :: wp = c_double
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: side
    integer(c_int),value, intent(in) :: uplo
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: diag
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int) :: stat
end function
!> cgemm GPU (cublas v2) version
pure function cublasCgemm(handle, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) result(stat) &
    bind(C,name="cublasCgemm_v2")
    import
    integer, parameter :: wp = c_float
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: transb
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: k
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    integer(c_int),value, intent(in) :: ldc
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    type(c_ptr), value, intent(in) :: beta
    type(c_ptr), value, intent(in) :: c
    integer(c_int) :: stat
end function
!> cgemv GPU (cublas v2) version
pure function cublasCgemv(handle, trans, m, n, alpha, a, lda, x, incx, beta, y, incy) result(stat) &
    bind(C,name="cublasCgemv_v2")
    import
    integer, parameter :: wp = c_float
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: trans
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: incx
    integer(c_int),value, intent(in) :: incy
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: x
    type(c_ptr), value, intent(in) :: beta
    type(c_ptr), value, intent(in) :: y
    integer(c_int) :: stat
end function
!> ctrmm GPU (cublas v2) version
pure function cublasCtrmm(handle, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) result(stat) &
    bind(C,name="cublasCtrmm_v2")
    import
    integer, parameter :: wp = c_float
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: side
    integer(c_int),value, intent(in) :: uplo
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: diag
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int) :: stat
end function
!> ctrsm GPU (cublas v2) version
pure function cublasCtrsm(handle, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) result(stat) &
    bind(C,name="cublasCtrsm_v2")
    import
    integer, parameter :: wp = c_float
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: side
    integer(c_int),value, intent(in) :: uplo
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: diag
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int) :: stat
end function
!> zgemm GPU (cublas v2) version
pure function cublasZgemm(handle, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) result(stat) &
    bind(C,name="cublasZgemm_v2")
    import
    integer, parameter :: wp = c_double
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: transb
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: k
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    integer(c_int),value, intent(in) :: ldc
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    type(c_ptr), value, intent(in) :: beta
    type(c_ptr), value, intent(in) :: c
    integer(c_int) :: stat
end function
!> zgemv GPU (cublas v2) version
pure function cublasZgemv(handle, trans, m, n, alpha, a, lda, x, incx, beta, y, incy) result(stat) &
    bind(C,name="cublasZgemv_v2")
    import
    integer, parameter :: wp = c_double
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: trans
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: incx
    integer(c_int),value, intent(in) :: incy
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: x
    type(c_ptr), value, intent(in) :: beta
    type(c_ptr), value, intent(in) :: y
    integer(c_int) :: stat
end function
!> ztrmm GPU (cublas v2) version
pure function cublasZtrmm(handle, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) result(stat) &
    bind(C,name="cublasZtrmm_v2")
    import
    integer, parameter :: wp = c_double
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: side
    integer(c_int),value, intent(in) :: uplo
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: diag
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int) :: stat
end function
!> ztrsm GPU (cublas v2) version
pure function cublasZtrsm(handle, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) result(stat) &
    bind(C,name="cublasZtrsm_v2")
    import
    integer, parameter :: wp = c_double
    type(c_ptr), value, intent(in) :: handle
    integer(c_int),value, intent(in) :: side
    integer(c_int),value, intent(in) :: uplo
    integer(c_int),value, intent(in) :: transa
    integer(c_int),value, intent(in) :: diag
    integer(c_int),value, intent(in) :: m
    integer(c_int),value, intent(in) :: n
    integer(c_int),value, intent(in) :: lda
    integer(c_int),value, intent(in) :: ldb
    type(c_ptr), value, intent(in) :: alpha
    type(c_ptr), value, intent(in) :: a
    type(c_ptr), value, intent(in) :: b
    integer(c_int) :: stat
end function
end interface
end module

