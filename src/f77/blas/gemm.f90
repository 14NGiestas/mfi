module f77_blas_gemm
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for GEMM.
!> Supports s, d, c, z.
!> See also: [[mfi_gemm]], [[sgemm]], [[dgemm]], [[cgemm]], [[zgemm]].
interface f77_gemm
!> Original interface for SGEMM
!> See also: [[mfi_gemm]], [[gemm]].
! sgemm CPU version
pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: b(ldb,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
!> Original interface for DGEMM
!> See also: [[mfi_gemm]], [[gemm]].
! dgemm CPU version
pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: b(ldb,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
!> Original interface for CGEMM
!> See also: [[mfi_gemm]], [[gemm]].
! cgemm CPU version
pure subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: b(ldb,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
!> Original interface for ZGEMM
!> See also: [[mfi_gemm]], [[gemm]].
! zgemm CPU version
pure subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: b(ldb,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
end module

