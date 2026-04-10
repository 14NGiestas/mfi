module mfi_blas_syr2k
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

!> Generic modern interface for SYR2K.
!> Supports s, d.
!> See also:
!> [[f77_syr2k:ssyr2k]], [[f77_syr2k:dsyr2k]].
interface mfi_syr2k
    module procedure :: mfi_ssyr2k
    module procedure :: mfi_dsyr2k
end interface

contains

!> Modern interface for [[f77_syr2k:f77_syr2k]].
!> See also: [[mfi_syr2k]], [[f77_syr2k]].
pure subroutine mfi_ssyr2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(in) :: b(:,:)
    real(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
    integer :: n, k, lda, ldb, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call f77_syr2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
!> Modern interface for [[f77_syr2k:f77_syr2k]].
!> See also: [[mfi_syr2k]], [[f77_syr2k]].
pure subroutine mfi_dsyr2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(in) :: b(:,:)
    real(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
    integer :: n, k, lda, ldb, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call f77_syr2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
end module

