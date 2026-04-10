module mfi_blas_tpsv
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for TPSV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_tpsv:stpsv]], [[f77_tpsv:dtpsv]], [[f77_tpsv:ctpsv]], [[f77_tpsv:ztpsv]].
interface mfi_tpsv
    module procedure :: mfi_stpsv
    module procedure :: mfi_dtpsv
    module procedure :: mfi_ctpsv
    module procedure :: mfi_ztpsv
end interface

contains

!> Modern interface for [[f77_tpsv:f77_tpsv]].
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
pure subroutine mfi_stpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: ap(:)
    real(REAL32), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
!> Modern interface for [[f77_tpsv:f77_tpsv]].
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
pure subroutine mfi_dtpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: ap(:)
    real(REAL64), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
!> Modern interface for [[f77_tpsv:f77_tpsv]].
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
pure subroutine mfi_ctpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: ap(:)
    complex(REAL32), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
!> Modern interface for [[f77_tpsv:f77_tpsv]].
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
pure subroutine mfi_ztpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: ap(:)
    complex(REAL64), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
end module

