module f77_blas_rotm
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for ROTM.
!> Supports s, d.
!> See also: [[mfi_rotm]], [[srotm]], [[drotm]].
interface f77_rotm
!> Original interface for SROTM
!> See also: [[mfi_rotm]], [[rotm]].
pure subroutine srotm(n, x, incx, y, incy, param)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    real(REAL32), intent(in) :: param(5)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DROTM
!> See also: [[mfi_rotm]], [[rotm]].
pure subroutine drotm(n, x, incx, y, incy, param)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    real(REAL64), intent(in) :: param(5)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

