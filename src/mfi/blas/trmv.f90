module mfi_blas_trmv
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

!> Generic modern interface for TRMV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_trmv:strmv]], [[f77_trmv:dtrmv]], [[f77_trmv:ctrmv]], [[f77_trmv:ztrmv]].
interface mfi_trmv
    module procedure :: mfi_strmv
    module procedure :: mfi_dtrmv
    module procedure :: mfi_ctrmv
    module procedure :: mfi_ztrmv
end interface

contains

!> Modern interface for [[f77_trmv:f77_trmv]].
!> See also: [[mfi_trmv]], [[f77_trmv]].
pure subroutine mfi_strmv(a, x, uplo, trans, diag, incx)
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
    integer :: n, lda
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_trmv:f77_trmv]].
!> See also: [[mfi_trmv]], [[f77_trmv]].
pure subroutine mfi_dtrmv(a, x, uplo, trans, diag, incx)
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
    integer :: n, lda
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_trmv:f77_trmv]].
!> See also: [[mfi_trmv]], [[f77_trmv]].
pure subroutine mfi_ctrmv(a, x, uplo, trans, diag, incx)
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
    integer :: n, lda
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_trmv:f77_trmv]].
!> See also: [[mfi_trmv]], [[f77_trmv]].
pure subroutine mfi_ztrmv(a, x, uplo, trans, diag, incx)
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
    integer :: n, lda
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
end module

