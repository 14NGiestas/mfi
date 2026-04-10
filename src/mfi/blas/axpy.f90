module mfi_blas_axpy
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for AXPY.
!> Supports s, d, c, z.
!> See also:
!> [[f77_axpy:saxpy]], [[f77_axpy:daxpy]], [[f77_axpy:caxpy]], [[f77_axpy:zaxpy]].
interface mfi_axpy
    module procedure :: mfi_saxpy
    module procedure :: mfi_daxpy
    module procedure :: mfi_caxpy
    module procedure :: mfi_zaxpy
end interface

contains

!> Modern interface for [[f77_axpy:f77_axpy]].
!> See also: [[mfi_axpy]], [[f77_axpy]].
pure subroutine mfi_saxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in), optional :: a
    real(REAL32) :: local_a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
    if (present(a)) then
        local_a = a
    else
        local_a = 1.0_wp
    end if
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
    call f77_axpy(n,local_a,x,local_incx,y,local_incy)
end subroutine
!> Modern interface for [[f77_axpy:f77_axpy]].
!> See also: [[mfi_axpy]], [[f77_axpy]].
pure subroutine mfi_daxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in), optional :: a
    real(REAL64) :: local_a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
    if (present(a)) then
        local_a = a
    else
        local_a = 1.0_wp
    end if
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
    call f77_axpy(n,local_a,x,local_incx,y,local_incy)
end subroutine
!> Modern interface for [[f77_axpy:f77_axpy]].
!> See also: [[mfi_axpy]], [[f77_axpy]].
pure subroutine mfi_caxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
    complex(REAL32), intent(in), optional :: a
    complex(REAL32) :: local_a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
    if (present(a)) then
        local_a = a
    else
        local_a = 1.0_wp
    end if
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
    call f77_axpy(n,local_a,x,local_incx,y,local_incy)
end subroutine
!> Modern interface for [[f77_axpy:f77_axpy]].
!> See also: [[mfi_axpy]], [[f77_axpy]].
pure subroutine mfi_zaxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
    complex(REAL64), intent(in), optional :: a
    complex(REAL64) :: local_a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
    if (present(a)) then
        local_a = a
    else
        local_a = 1.0_wp
    end if
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
    call f77_axpy(n,local_a,x,local_incx,y,local_incy)
end subroutine
end module

