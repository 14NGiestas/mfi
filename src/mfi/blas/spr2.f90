module mfi_blas_spr2
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for SPR2.
!> Supports s, d.
!> See also:
!> [[f77_spr2:sspr2]], [[f77_spr2:dspr2]].
interface mfi_spr2
    module procedure :: mfi_sspr2
    module procedure :: mfi_dspr2
end interface

contains

!> Modern interface for [[f77_spr2:f77_spr2]].
!> See also: [[mfi_spr2]], [[f77_spr2]].
pure subroutine mfi_sspr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: y(:)
    real(REAL32), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
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
    call f77_spr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
!> Modern interface for [[f77_spr2:f77_spr2]].
!> See also: [[mfi_spr2]], [[f77_spr2]].
pure subroutine mfi_dspr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: y(:)
    real(REAL64), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
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
    call f77_spr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
end module

