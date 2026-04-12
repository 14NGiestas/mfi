module f77_blas_dot
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for DOT.
!> Supports s, d.
!> See also: [[mfi_dot]], [[sdot]], [[ddot]].
interface f77_dot
!> Original interface for SDOT
!> See also: [[mfi_dot]], [[dot]].
pure function sdot(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32) :: sdot
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
!> Original interface for DDOT
!> See also: [[mfi_dot]], [[dot]].
pure function ddot(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64) :: ddot
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
end interface
end module

