module mfi_blas_hpr
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

!> Generic modern interface for HPR.
!> Supports c, z.
!> See also:
!> [[f77_hpr:chpr]], [[f77_hpr:zhpr]].
interface mfi_hpr
    module procedure :: mfi_chpr
    module procedure :: mfi_zhpr
end interface

contains

!> Modern interface for [[f77_hpr:f77_hpr]].
!> See also: [[mfi_hpr]], [[f77_hpr]].
pure subroutine mfi_chpr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
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
    n = size(x)
    call f77_hpr(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
!> Modern interface for [[f77_hpr:f77_hpr]].
!> See also: [[mfi_hpr]], [[f77_hpr]].
pure subroutine mfi_zhpr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
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
    n = size(x)
    call f77_hpr(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
end module

