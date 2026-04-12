module mfi_blas_trsv
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

!> Generic modern interface for TRSV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_trsv:strsv]], [[f77_trsv:dtrsv]], [[f77_trsv:ctrsv]], [[f77_trsv:ztrsv]].
interface mfi_trsv
    module procedure :: mfi_strsv
    module procedure :: mfi_dtrsv
    module procedure :: mfi_ctrsv
    module procedure :: mfi_ztrsv
end interface

contains

!> Modern interface for [[f77_trsv:f77_trsv]].
!> See also: [[mfi_trsv]], [[f77_trsv]].
pure subroutine mfi_strsv(a, x, uplo, trans, diag, incx)
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
    call f77_trsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_trsv:f77_trsv]].
!> See also: [[mfi_trsv]], [[f77_trsv]].
pure subroutine mfi_dtrsv(a, x, uplo, trans, diag, incx)
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
    call f77_trsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_trsv:f77_trsv]].
!> See also: [[mfi_trsv]], [[f77_trsv]].
pure subroutine mfi_ctrsv(a, x, uplo, trans, diag, incx)
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
    call f77_trsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
!> Modern interface for [[f77_trsv:f77_trsv]].
!> See also: [[mfi_trsv]], [[f77_trsv]].
pure subroutine mfi_ztrsv(a, x, uplo, trans, diag, incx)
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
    call f77_trsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
end module

