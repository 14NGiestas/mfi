module mfi_blas_spmv
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for SPMV.
!> Supports s, d.
!> See also:
!> [[f77_spmv:sspmv]], [[f77_spmv:dspmv]].
interface mfi_spmv
    module procedure :: mfi_sspmv
    module procedure :: mfi_dspmv
end interface

contains

!> Modern interface for [[f77_spmv:f77_spmv]].
!> See also: [[mfi_spmv]], [[f77_spmv]].
pure subroutine mfi_sspmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: ap(:)
    real(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
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
    call f77_spmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_spmv:f77_spmv]].
!> See also: [[mfi_spmv]], [[f77_spmv]].
pure subroutine mfi_dspmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: ap(:)
    real(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
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
    call f77_spmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
end module

