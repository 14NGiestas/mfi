module mfi_blas_copy
    use iso_fortran_env
    use f77_blas
#if defined(MFI_CUBLAS)
    use iso_c_binding
    use mfi_blas_cublas
#endif
#if defined(MFI_EXTENSIONS)
    use mfi_blas_extensions
#endif
    implicit none

!> Generic modern interface for COPY.
!> Supports s, d, c, z.
!> See also:
!> [[f77_copy:scopy]], [[f77_copy:dcopy]], [[f77_copy:ccopy]], [[f77_copy:zcopy]].
interface mfi_copy
    module procedure :: mfi_scopy
    module procedure :: mfi_dcopy
    module procedure :: mfi_ccopy
    module procedure :: mfi_zcopy
end interface

contains

!> Modern interface for [[f77_copy:f77_copy]].
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine mfi_scopy(x, y, incx, incy)
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
    call f77_copy(n,x,local_incx,y,local_incy)
end subroutine
!> Modern interface for [[f77_copy:f77_copy]].
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine mfi_dcopy(x, y, incx, incy)
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
    call f77_copy(n,x,local_incx,y,local_incy)
end subroutine
!> Modern interface for [[f77_copy:f77_copy]].
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine mfi_ccopy(x, y, incx, incy)
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
    call f77_copy(n,x,local_incx,y,local_incy)
end subroutine
!> Modern interface for [[f77_copy:f77_copy]].
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine mfi_zcopy(x, y, incx, incy)
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
    call f77_copy(n,x,local_incx,y,local_incy)
end subroutine
end module

