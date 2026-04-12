module f77_blas_symm
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SYMM.
!> Supports s, d.
!> See also: [[mfi_symm]], [[ssymm]], [[dsymm]].
interface f77_symm
!> Original interface for SSYMM
!> See also: [[mfi_symm]], [[symm]].
pure subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: b(ldb,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
!> Original interface for DSYMM
!> See also: [[mfi_symm]], [[symm]].
pure subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: b(ldb,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
end module

