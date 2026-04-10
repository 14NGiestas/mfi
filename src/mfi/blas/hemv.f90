module mfi_blas_hemv
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for HEMV.
!> Supports c, z.
!> See also:
!> [[f77_hemv:chemv]], [[f77_hemv:zhemv]].
interface mfi_hemv
    module procedure :: mfi_chemv
    module procedure :: mfi_zhemv
end interface

contains

!> Modern interface for [[f77_hemv:f77_hemv]].
!> See also: [[mfi_hemv]], [[f77_hemv]].
pure subroutine mfi_chemv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: a(:,:)
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
    integer :: n, lda
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_hemv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_hemv:f77_hemv]].
!> See also: [[mfi_hemv]], [[f77_hemv]].
pure subroutine mfi_zhemv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: a(:,:)
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
    integer :: n, lda
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_hemv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
end module

