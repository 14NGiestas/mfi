module mfi_blas_tpmv
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for TPMV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_tpmv:stpmv]], [[f77_tpmv:dtpmv]], [[f77_tpmv:ctpmv]], [[f77_tpmv:ztpmv]].
interface mfi_tpmv
    module procedure :: mfi_stpmv
    module procedure :: mfi_dtpmv
    module procedure :: mfi_ctpmv
    module procedure :: mfi_ztpmv
end interface

contains

!> Modern interface for [[f77_tpmv:f77_tpmv]].
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
pure subroutine mfi_stpmv(ap, x, uplo, trans, diag, incx)
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
    call f77_tpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
!> Modern interface for [[f77_tpmv:f77_tpmv]].
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
pure subroutine mfi_dtpmv(ap, x, uplo, trans, diag, incx)
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
    call f77_tpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
!> Modern interface for [[f77_tpmv:f77_tpmv]].
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
pure subroutine mfi_ctpmv(ap, x, uplo, trans, diag, incx)
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
    call f77_tpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
!> Modern interface for [[f77_tpmv:f77_tpmv]].
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
pure subroutine mfi_ztpmv(ap, x, uplo, trans, diag, incx)
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
    call f77_tpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
end module

