module f77_blas_herk
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for HERK.
!> Supports c, z.
!> See also: [[mfi_herk]], [[cherk]], [[zherk]].
interface f77_herk
!> Original interface for CHERK
!> See also: [[mfi_herk]], [[herk]].
pure subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
!> Original interface for ZHERK
!> See also: [[mfi_herk]], [[herk]].
pure subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
end interface
end module

