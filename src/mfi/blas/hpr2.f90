module mfi_blas_hpr2
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for HPR2.
!> Supports c, z.
!> See also:
!> [[f77_hpr2:chpr2]], [[f77_hpr2:zhpr2]].
interface mfi_hpr2
    module procedure :: mfi_chpr2
    module procedure :: mfi_zhpr2
end interface

contains

!> Modern interface for [[f77_hpr2:f77_hpr2]].
!> See also: [[mfi_hpr2]], [[f77_hpr2]].
pure subroutine mfi_chpr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: y(:)
    complex(REAL32), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
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
    call f77_hpr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
!> Modern interface for [[f77_hpr2:f77_hpr2]].
!> See also: [[mfi_hpr2]], [[f77_hpr2]].
pure subroutine mfi_zhpr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: y(:)
    complex(REAL64), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
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
    call f77_hpr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
end module

