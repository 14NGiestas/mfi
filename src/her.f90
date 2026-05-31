module f77_blas_her
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for HER.
!> Supports c, z.
!> See also: [[mfi_her]], [[cher]], [[zher]].
interface f77_her
!> Original interface for CHER
!> See also: [[mfi_her]], [[her]].
pure subroutine cher(uplo, n, alpha, x, incx, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
!> Original interface for ZHER
!> See also: [[mfi_her]], [[her]].
pure subroutine zher(uplo, n, alpha, x, incx, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
end module

