module f77_blas_syrk
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SYRK.
!> Supports s, d.
!> See also: [[mfi_syrk]], [[ssyrk]], [[dsyrk]].
interface f77_syrk
!> Original interface for SSYRK
!> See also: [[mfi_syrk]], [[syrk]].
pure subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
!> Original interface for DSYRK
!> See also: [[mfi_syrk]], [[syrk]].
pure subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
end interface
end module

