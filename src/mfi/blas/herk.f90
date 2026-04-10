module mfi_blas_herk
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for HERK.
!> Supports c, z.
!> See also:
!> [[f77_herk:cherk]], [[f77_herk:zherk]].
interface mfi_herk
    module procedure :: mfi_cherk
    module procedure :: mfi_zherk
end interface

contains

!> Modern interface for [[f77_herk:f77_herk]].
!> See also: [[mfi_herk]], [[f77_herk]].
pure subroutine mfi_cherk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
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
    call f77_herk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
!> Modern interface for [[f77_herk:f77_herk]].
!> See also: [[mfi_herk]], [[f77_herk]].
pure subroutine mfi_zherk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
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
    call f77_herk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
end module

