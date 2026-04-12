module f77_blas_trmm
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for TRMM.
!> Supports s, d, c, z.
!> See also: [[mfi_trmm]], [[strmm]], [[dtrmm]], [[ctrmm]], [[ztrmm]].
interface f77_trmm
!> Original interface for STRMM
!> See also: [[mfi_trmm]], [[trmm]].
pure subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
!> Original interface for DTRMM
!> See also: [[mfi_trmm]], [[trmm]].
pure subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
!> Original interface for CTRMM
!> See also: [[mfi_trmm]], [[trmm]].
pure subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
!> Original interface for ZTRMM
!> See also: [[mfi_trmm]], [[trmm]].
pure subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
end interface
end module

