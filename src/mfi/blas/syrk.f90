module mfi_blas_syrk
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for SYRK.
!> Supports s, d.
!> See also:
!> [[f77_syrk:ssyrk]], [[f77_syrk:dsyrk]].
interface mfi_syrk
    module procedure :: mfi_ssyrk
    module procedure :: mfi_dsyrk
end interface

contains

!> Modern interface for [[f77_syrk:f77_syrk]].
!> See also: [[mfi_syrk]], [[f77_syrk]].
pure subroutine mfi_ssyrk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
    integer :: n, k, lda, ldc
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
    ldc = max(1,size(c,1))
    call f77_syrk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
!> Modern interface for [[f77_syrk:f77_syrk]].
!> See also: [[mfi_syrk]], [[f77_syrk]].
pure subroutine mfi_dsyrk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
    integer :: n, k, lda, ldc
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
    ldc = max(1,size(c,1))
    call f77_syrk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
end module

