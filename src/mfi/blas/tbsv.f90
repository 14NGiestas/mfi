module mfi_blas_tbsv
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for TBSV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_tbsv:stbsv]], [[f77_tbsv:dtbsv]], [[f77_tbsv:ctbsv]], [[f77_tbsv:ztbsv]].
interface mfi_tbsv
    module procedure :: mfi_stbsv
    module procedure :: mfi_dtbsv
    module procedure :: mfi_ctbsv
    module procedure :: mfi_ztbsv
end interface

contains

!> Modern interface for [[f77_tbsv:f77_tbsv]].
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
pure subroutine mfi_stbsv(a, x, uplo, trans, diag, incx)
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
    call f77_tbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_tbsv:f77_tbsv]].
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
pure subroutine mfi_dtbsv(a, x, uplo, trans, diag, incx)
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
    call f77_tbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_tbsv:f77_tbsv]].
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
pure subroutine mfi_ctbsv(a, x, uplo, trans, diag, incx)
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
    call f77_tbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_tbsv:f77_tbsv]].
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
pure subroutine mfi_ztbsv(a, x, uplo, trans, diag, incx)
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
    call f77_tbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
end module

