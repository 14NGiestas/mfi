module f77_blas_ger
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for GER.
!> Supports s, d.
!> See also: [[mfi_ger]], [[sger]], [[dger]].
interface f77_ger
!> Original interface for SGER
!> See also: [[mfi_ger]], [[ger]].
pure subroutine sger(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DGER
!> See also: [[mfi_ger]], [[ger]].
pure subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

