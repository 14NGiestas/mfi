module f77_blas_syr
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SYR.
!> Supports s, d.
!> See also: [[mfi_syr]], [[ssyr]], [[dsyr]].
interface f77_syr
!> Original interface for SSYR
!> See also: [[mfi_syr]], [[syr]].
pure subroutine ssyr(uplo, n, alpha, x, incx, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
!> Original interface for DSYR
!> See also: [[mfi_syr]], [[syr]].
pure subroutine dsyr(uplo, n, alpha, x, incx, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
end module

