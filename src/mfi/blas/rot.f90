module mfi_blas_rot
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

!> Generic modern interface for ROT.
!> Supports s, d, c, z, cs, zd.
!> See also:
!> [[f77_rot:srot]], [[f77_rot:drot]], [[f77_rot:crot]], [[f77_rot:zrot]], [[f77_rot:csrot]], [[f77_rot:zdrot]].
interface mfi_rot
    module procedure :: mfi_srot
    module procedure :: mfi_drot
    module procedure :: mfi_crot
    module procedure :: mfi_zrot
    module procedure :: mfi_csrot
    module procedure :: mfi_zdrot
end interface

contains

!> Modern interface for [[f77_rot:f77_rot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - s*xi
!>```
pure subroutine mfi_srot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: x(:)
    real(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in) :: c
    real(REAL32), intent(in) :: s
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
    n = size(x)
    call f77_rot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Modern interface for [[f77_rot:f77_rot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - s*xi
!>```
pure subroutine mfi_drot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: x(:)
    real(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in) :: c
    real(REAL64), intent(in) :: s
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
    n = size(x)
    call f77_rot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Modern interface for [[f77_rot:f77_rot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_crot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in) :: c
    complex(REAL32), intent(in) :: s
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
    n = size(x)
    call f77_rot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Modern interface for [[f77_rot:f77_rot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_zrot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in) :: c
    complex(REAL64), intent(in) :: s
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
    n = size(x)
    call f77_rot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Modern interface for [[f77_rot:f77_rot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_csrot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in) :: c
    real(REAL32), intent(in) :: s
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
    n = size(x)
    call f77_rot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Modern interface for [[f77_rot:f77_rot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_zdrot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in) :: c
    real(REAL64), intent(in) :: s
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
    n = size(x)
    call f77_rot(n,x,local_incx,y,local_incy,c,s)
end subroutine
end module

