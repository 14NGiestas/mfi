module mfi_blas_swap
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for SWAP.
!> Supports s, d, c, z.
!> See also:
!> [[f77_swap:sswap]], [[f77_swap:dswap]], [[f77_swap:cswap]], [[f77_swap:zswap]].
interface mfi_swap
    module procedure :: mfi_sswap
    module procedure :: mfi_dswap
    module procedure :: mfi_cswap
    module procedure :: mfi_zswap
end interface

contains

!> Modern interface for [[f77_swap:f77_swap]].
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine mfi_sswap(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: y(:)
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
    call f77_swap(n,x,local_incx,y,local_incy)
end subroutine
!> Modern interface for [[f77_swap:f77_swap]].
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine mfi_dswap(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: y(:)
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
    call f77_swap(n,x,local_incx,y,local_incy)
end subroutine
!> Modern interface for [[f77_swap:f77_swap]].
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine mfi_cswap(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
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
    call f77_swap(n,x,local_incx,y,local_incy)
end subroutine
!> Modern interface for [[f77_swap:f77_swap]].
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine mfi_zswap(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
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
    call f77_swap(n,x,local_incx,y,local_incy)
end subroutine
end module

