module f77_blas_tbmv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for TBMV.
!> Supports s, d, c, z.
!> See also: [[mfi_tbmv]], [[stbmv]], [[dtbmv]], [[ctbmv]], [[ztbmv]].
interface f77_tbmv
!> Original interface for STBMV
!> See also: [[mfi_tbmv]], [[tbmv]].
pure subroutine stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
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
!> Original interface for DTBMV
!> See also: [[mfi_tbmv]], [[tbmv]].
pure subroutine dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
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
!> Original interface for CTBMV
!> See also: [[mfi_tbmv]], [[tbmv]].
pure subroutine ctbmv(uplo, trans, diag, n, k, a, lda, x, incx)
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
!> Original interface for ZTBMV
!> See also: [[mfi_tbmv]], [[tbmv]].
pure subroutine ztbmv(uplo, trans, diag, n, k, a, lda, x, incx)
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

