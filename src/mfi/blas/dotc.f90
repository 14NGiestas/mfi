module mfi_blas_dotc
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for DOTC.
!> Supports c, z.
!> See also:
!> [[f77_dotc:cdotc]], [[f77_dotc:zdotc]].
interface mfi_dotc
    module procedure :: mfi_cdotc
    module procedure :: mfi_zdotc
end interface

contains

!> Modern interface for [[f77_dotc:f77_dotc]].
!> See also: [[mfi_dotc]], [[f77_dotc]].
pure function mfi_cdotc(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32) :: mfi_cdotc
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: y(:)
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
    mfi_cdotc = f77_dotc(n,x,local_incx,y,local_incy)
end function
!> Modern interface for [[f77_dotc:f77_dotc]].
!> See also: [[mfi_dotc]], [[f77_dotc]].
pure function mfi_zdotc(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64) :: mfi_zdotc
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: y(:)
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
    mfi_zdotc = f77_dotc(n,x,local_incx,y,local_incy)
end function
end module

