module mfi_blas_hemm
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

!> Generic modern interface for HEMM.
!> Supports c, z.
!> See also:
!> [[f77_hemm:chemm]], [[f77_hemm:zhemm]].
interface mfi_hemm
    module procedure :: mfi_chemm
    module procedure :: mfi_zhemm
end interface

contains

!> Modern interface for [[f77_hemm:f77_hemm]].
!> See also: [[mfi_hemm]], [[f77_hemm]].
pure subroutine mfi_chemm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(in) :: b(:,:)
    complex(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
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
    call f77_hemm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
!> Modern interface for [[f77_hemm:f77_hemm]].
!> See also: [[mfi_hemm]], [[f77_hemm]].
pure subroutine mfi_zhemm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(in) :: b(:,:)
    complex(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
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
    call f77_hemm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
end module

