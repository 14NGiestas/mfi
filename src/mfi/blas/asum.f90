module mfi_blas_asum
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for ASUM.
!> Supports s, d, sc, dz.
!> See also:
!> [[f77_asum:sasum]], [[f77_asum:dasum]], [[f77_asum:scasum]], [[f77_asum:dzasum]].
interface mfi_asum
    module procedure :: mfi_sasum
    module procedure :: mfi_dasum
    module procedure :: mfi_scasum
    module procedure :: mfi_dzasum
end interface

contains

!> Modern interface for [[f77_asum:f77_asum]].
!> See also: [[mfi_asum]], [[f77_asum]].
pure function mfi_sasum(x, incx)
    real(REAL32) :: mfi_sasum
    real(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_sasum = f77_asum(n, x, local_incx)
end function
!> Modern interface for [[f77_asum:f77_asum]].
!> See also: [[mfi_asum]], [[f77_asum]].
pure function mfi_dasum(x, incx)
    real(REAL64) :: mfi_dasum
    real(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dasum = f77_asum(n, x, local_incx)
end function
!> Modern interface for [[f77_asum:f77_asum]].
!> See also: [[mfi_asum]], [[f77_asum]].
pure function mfi_scasum(x, incx)
    real(REAL32) :: mfi_scasum
    complex(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_scasum = f77_asum(n, x, local_incx)
end function
!> Modern interface for [[f77_asum:f77_asum]].
!> See also: [[mfi_asum]], [[f77_asum]].
pure function mfi_dzasum(x, incx)
    real(REAL64) :: mfi_dzasum
    complex(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dzasum = f77_asum(n, x, local_incx)
end function
end module

