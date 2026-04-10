module mfi_blas_symm
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for SYMM.
!> Supports s, d.
!> See also:
!> [[f77_symm:ssymm]], [[f77_symm:dsymm]].
interface mfi_symm
    module procedure :: mfi_ssymm
    module procedure :: mfi_dsymm
end interface

contains

!> Modern interface for [[f77_symm:f77_symm]].
!> See also: [[mfi_symm]], [[f77_symm]].
pure subroutine mfi_ssymm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(in) :: b(:,:)
    real(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
    integer :: m, n, lda, ldb, ldc
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call f77_symm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
!> Modern interface for [[f77_symm:f77_symm]].
!> See also: [[mfi_symm]], [[f77_symm]].
pure subroutine mfi_dsymm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(in) :: b(:,:)
    real(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
    integer :: m, n, lda, ldb, ldc
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call f77_symm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
end module

