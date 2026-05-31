module f77_blas_syr2
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SYR2.
!> Supports s, d.
!> See also: [[mfi_syr2]], [[ssyr2]], [[dsyr2]].
interface f77_syr2
!> Original interface for SSYR2
!> See also: [[mfi_syr2]], [[syr2]].
pure subroutine ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    real(REAL32), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DSYR2
!> See also: [[mfi_syr2]], [[syr2]].
pure subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    real(REAL64), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

