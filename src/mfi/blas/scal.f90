module mfi_blas_scal
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for SCAL.
!> Supports s, d, c, z, cs, zd.
!> See also:
!> [[f77_scal:sscal]], [[f77_scal:dscal]], [[f77_scal:cscal]], [[f77_scal:zscal]], [[f77_scal:csscal]], [[f77_scal:zdscal]].
interface mfi_scal
    module procedure :: mfi_sscal
    module procedure :: mfi_dscal
    module procedure :: mfi_cscal
    module procedure :: mfi_zscal
    module procedure :: mfi_csscal
    module procedure :: mfi_zdscal
end interface

contains

!> Modern interface for [[f77_scal:f77_scal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_SSCAL scales a vector by a constant.
pure subroutine mfi_sscal(a, x, incx)
    real(REAL32), intent(inout) :: x(:)
    real(REAL32), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_scal(n,a,x,local_incx)
end subroutine
!> Modern interface for [[f77_scal:f77_scal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_DSCAL scales a vector by a constant.
pure subroutine mfi_dscal(a, x, incx)
    real(REAL64), intent(inout) :: x(:)
    real(REAL64), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_scal(n,a,x,local_incx)
end subroutine
!> Modern interface for [[f77_scal:f77_scal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_CSCAL scales a vector by a constant.
pure subroutine mfi_cscal(a, x, incx)
    complex(REAL32), intent(inout) :: x(:)
    complex(REAL32), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_scal(n,a,x,local_incx)
end subroutine
!> Modern interface for [[f77_scal:f77_scal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_ZSCAL scales a vector by a constant.
pure subroutine mfi_zscal(a, x, incx)
    complex(REAL64), intent(inout) :: x(:)
    complex(REAL64), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_scal(n,a,x,local_incx)
end subroutine
!> Modern interface for [[f77_scal:f77_scal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_CSSCAL scales a vector by a constant.
pure subroutine mfi_csscal(a, x, incx)
    complex(REAL32), intent(inout) :: x(:)
    real(REAL32), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_scal(n,a,x,local_incx)
end subroutine
!> Modern interface for [[f77_scal:f77_scal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_ZDSCAL scales a vector by a constant.
pure subroutine mfi_zdscal(a, x, incx)
    complex(REAL64), intent(inout) :: x(:)
    real(REAL64), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_scal(n,a,x,local_incx)
end subroutine
end module

