module f77_blas_spmv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SPMV.
!> Supports s, d.
!> See also: [[mfi_spmv]], [[sspmv]], [[dspmv]].
interface f77_spmv
!> Original interface for SSPMV
!> See also: [[mfi_spmv]], [[spmv]].
pure subroutine sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: ap(*)
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DSPMV
!> See also: [[mfi_spmv]], [[spmv]].
pure subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: ap(*)
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

