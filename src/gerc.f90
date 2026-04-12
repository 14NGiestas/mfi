module mfi_blas_gerc
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

!> Generic modern interface for GERC.
!> Supports c, z.
!> See also:
!> [[f77_gerc:cgerc]], [[f77_gerc:zgerc]].
interface mfi_gerc
    module procedure :: mfi_cgerc
    module procedure :: mfi_zgerc
end interface

contains

!> Modern interface for [[f77_gerc:f77_gerc]].
!> See also: [[mfi_gerc]], [[f77_gerc]].
pure subroutine mfi_cgerc(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: y(:)
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
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
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call f77_gerc(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
!> Modern interface for [[f77_gerc:f77_gerc]].
!> See also: [[mfi_gerc]], [[f77_gerc]].
pure subroutine mfi_zgerc(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: y(:)
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
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
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call f77_gerc(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
end module

