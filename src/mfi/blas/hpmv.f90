module mfi_blas_hpmv
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

!> Generic modern interface for HPMV.
!> Supports c, z.
!> See also:
!> [[f77_hpmv:chpmv]], [[f77_hpmv:zhpmv]].
interface mfi_hpmv
    module procedure :: mfi_chpmv
    module procedure :: mfi_zhpmv
end interface

contains

!> Modern interface for [[f77_hpmv:f77_hpmv]].
!> See also: [[mfi_hpmv]], [[f77_hpmv]].
pure subroutine mfi_chpmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: ap(:)
    complex(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
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
    n = size(x)
    call f77_hpmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_hpmv:f77_hpmv]].
!> See also: [[mfi_hpmv]], [[f77_hpmv]].
pure subroutine mfi_zhpmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: ap(:)
    complex(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
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
    n = size(x)
    call f77_hpmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
end module

