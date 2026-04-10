module mfi_blas_syr
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

!> Generic modern interface for SYR.
!> Supports s, d.
!> See also:
!> [[f77_syr:ssyr]], [[f77_syr:dsyr]].
interface mfi_syr
    module procedure :: mfi_ssyr
    module procedure :: mfi_dsyr
end interface

contains

!> Modern interface for [[f77_syr:f77_syr]].
!> See also: [[mfi_syr]], [[f77_syr]].
pure subroutine mfi_ssyr(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
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
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_syr(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
!> Modern interface for [[f77_syr:f77_syr]].
!> See also: [[mfi_syr]], [[f77_syr]].
pure subroutine mfi_dsyr(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
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
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_syr(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
end module

