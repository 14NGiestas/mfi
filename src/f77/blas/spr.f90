module f77_blas_spr
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SPR.
!> Supports s, d.
!> See also: [[mfi_spr]], [[sspr]], [[dspr]].
interface f77_spr
!> Original interface for SSPR
!> See also: [[mfi_spr]], [[spr]].
pure subroutine sspr(uplo, n, alpha, x, incx, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for DSPR
!> See also: [[mfi_spr]], [[spr]].
pure subroutine dspr(uplo, n, alpha, x, incx, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
end module

