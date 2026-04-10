module mfi_blas_rotm
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for ROTM.
!> Supports s, d.
!> See also:
!> [[f77_rotm:srotm]], [[f77_rotm:drotm]].
interface mfi_rotm
    module procedure :: mfi_srotm
    module procedure :: mfi_drotm
end interface

contains

!> Modern interface for [[f77_rotm:f77_rotm]].
!> See also: [[mfi_rotm]], [[f77_rotm]].
pure subroutine mfi_srotm(x, y, param, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: x(:)
    real(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in) :: param(5)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    call f77_rotm(n,x,local_incx,y,local_incy,param)
end subroutine
!> Modern interface for [[f77_rotm:f77_rotm]].
!> See also: [[mfi_rotm]], [[f77_rotm]].
pure subroutine mfi_drotm(x, y, param, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: x(:)
    real(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in) :: param(5)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    call f77_rotm(n,x,local_incx,y,local_incy,param)
end subroutine
end module

