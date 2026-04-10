module mfi_blas_tbmv
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for TBMV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_tbmv:stbmv]], [[f77_tbmv:dtbmv]], [[f77_tbmv:ctbmv]], [[f77_tbmv:ztbmv]].
interface mfi_tbmv
    module procedure :: mfi_stbmv
    module procedure :: mfi_dtbmv
    module procedure :: mfi_ctbmv
    module procedure :: mfi_ztbmv
end interface

contains

!> Modern interface for [[f77_tbmv:f77_tbmv]].
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
pure subroutine mfi_stbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
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
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_tbmv:f77_tbmv]].
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
pure subroutine mfi_dtbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
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
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_tbmv:f77_tbmv]].
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
pure subroutine mfi_ctbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
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
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_tbmv:f77_tbmv]].
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
pure subroutine mfi_ztbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
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
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
end module

