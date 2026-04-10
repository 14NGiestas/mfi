module f77_blas_spr2
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SPR2.
!> Supports s, d.
!> See also: [[mfi_spr2]], [[sspr2]], [[dspr2]].
interface f77_spr2
!> Original interface for SSPR2
!> See also: [[mfi_spr2]], [[spr2]].
pure subroutine sspr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    real(REAL32), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DSPR2
!> See also: [[mfi_spr2]], [[spr2]].
pure subroutine dspr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    real(REAL64), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

