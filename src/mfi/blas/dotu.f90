module mfi_blas_dotu
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for DOTU.
!> Supports c, z.
!> See also:
!> [[f77_dotu:cdotu]], [[f77_dotu:zdotu]].
interface mfi_dotu
    module procedure :: mfi_cdotu
    module procedure :: mfi_zdotu
end interface

contains

!> Modern interface for [[f77_dotu:f77_dotu]].
!> See also: [[mfi_dotu]], [[f77_dotu]].
pure function mfi_cdotu(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32) :: mfi_cdotu
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
    mfi_cdotu = f77_dotu(n,x,local_incx,y,local_incy)
end function
!> Modern interface for [[f77_dotu:f77_dotu]].
!> See also: [[mfi_dotu]], [[f77_dotu]].
pure function mfi_zdotu(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64) :: mfi_zdotu
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
    mfi_zdotu = f77_dotu(n,x,local_incx,y,local_incy)
end function
end module

