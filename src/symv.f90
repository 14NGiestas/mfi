module f77_blas_symv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SYMV.
!> Supports s, d.
!> See also: [[mfi_symv]], [[ssymv]], [[dsymv]].
interface f77_symv
!> Original interface for SSYMV
!> See also: [[mfi_symv]], [[symv]].
pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DSYMV
!> See also: [[mfi_symv]], [[symv]].
pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

