module f77_blas_trmv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for TRMV.
!> Supports s, d, c, z.
!> See also: [[mfi_trmv]], [[strmv]], [[dtrmv]], [[ctrmv]], [[ztrmv]].
interface f77_trmv
!> Original interface for STRMV
!> See also: [[mfi_trmv]], [[trmv]].
pure subroutine strmv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
!> Original interface for DTRMV
!> See also: [[mfi_trmv]], [[trmv]].
pure subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
!> Original interface for CTRMV
!> See also: [[mfi_trmv]], [[trmv]].
pure subroutine ctrmv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
!> Original interface for ZTRMV
!> See also: [[mfi_trmv]], [[trmv]].
pure subroutine ztrmv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
end module

