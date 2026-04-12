module f77_blas_tbsv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for TBSV.
!> Supports s, d, c, z.
!> See also: [[mfi_tbsv]], [[stbsv]], [[dtbsv]], [[ctbsv]], [[ztbsv]].
interface f77_tbsv
!> Original interface for STBSV
!> See also: [[mfi_tbsv]], [[tbsv]].
pure subroutine stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
!> Original interface for DTBSV
!> See also: [[mfi_tbsv]], [[tbsv]].
pure subroutine dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
!> Original interface for CTBSV
!> See also: [[mfi_tbsv]], [[tbsv]].
pure subroutine ctbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
!> Original interface for ZTBSV
!> See also: [[mfi_tbsv]], [[tbsv]].
pure subroutine ztbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
end module

