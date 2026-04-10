module mfi_blas_dot
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for DOT.
!> Supports s, d.
!> See also:
!> [[f77_dot:sdot]], [[f77_dot:ddot]].
interface mfi_dot
    module procedure :: mfi_sdot
    module procedure :: mfi_ddot
end interface

contains

!> Modern interface for [[f77_dot:f77_dot]].
!> See also: [[mfi_dot]], [[f77_dot]].
pure function mfi_sdot(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32) :: mfi_sdot
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: y(:)
    integer :: n
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    if (present(incy)) then
        local_incy = incy
    else
        local_incy = 1
    end if
    N = size(X)
    mfi_sdot = f77_dot(n,x,local_incx,y,local_incy)
end function
!> Modern interface for [[f77_dot:f77_dot]].
!> See also: [[mfi_dot]], [[f77_dot]].
pure function mfi_ddot(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64) :: mfi_ddot
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: y(:)
    integer :: n
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    if (present(incy)) then
        local_incy = incy
    else
        local_incy = 1
    end if
    N = size(X)
    mfi_ddot = f77_dot(n,x,local_incx,y,local_incy)
end function
end module

