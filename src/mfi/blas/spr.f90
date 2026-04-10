module mfi_blas_spr
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for SPR.
!> Supports s, d.
!> See also:
!> [[f77_spr:sspr]], [[f77_spr:dspr]].
interface mfi_spr
    module procedure :: mfi_sspr
    module procedure :: mfi_dspr
end interface

contains

!> Modern interface for [[f77_spr:f77_spr]].
!> See also: [[mfi_spr]], [[f77_spr]].
pure subroutine mfi_sspr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
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
    call f77_spr(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
!> Modern interface for [[f77_spr:f77_spr]].
!> See also: [[mfi_spr]], [[f77_spr]].
pure subroutine mfi_dspr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
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
    call f77_spr(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
end module

