module f77_blas_trsv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for TRSV.
!> Supports s, d, c, z.
!> See also: [[mfi_trsv]], [[strsv]], [[dtrsv]], [[ctrsv]], [[ztrsv]].
interface f77_trsv
!> Original interface for STRSV
!> See also: [[mfi_trsv]], [[trsv]].
pure subroutine strsv(uplo, trans, diag, n, a, lda, x, incx)
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
!> Original interface for DTRSV
!> See also: [[mfi_trsv]], [[trsv]].
pure subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
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
!> Original interface for CTRSV
!> See also: [[mfi_trsv]], [[trsv]].
pure subroutine ctrsv(uplo, trans, diag, n, a, lda, x, incx)
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
!> Original interface for ZTRSV
!> See also: [[mfi_trsv]], [[trsv]].
pure subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
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

