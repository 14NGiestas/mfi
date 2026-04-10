module f77_blas_her2k
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for HER2K.
!> Supports c, z.
!> See also: [[mfi_her2k]], [[cher2k]], [[zher2k]].
interface f77_her2k
!> Original interface for CHER2K
!> See also: [[mfi_her2k]], [[her2k]].
pure subroutine cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: b(ldb,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
!> Original interface for ZHER2K
!> See also: [[mfi_her2k]], [[her2k]].
pure subroutine zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: b(ldb,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
end module

