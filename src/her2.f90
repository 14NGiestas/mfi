module f77_blas_her2
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for HER2.
!> Supports c, z.
!> See also: [[mfi_her2]], [[cher2]], [[zher2]].
interface f77_her2
!> Original interface for CHER2
!> See also: [[mfi_her2]], [[her2]].
pure subroutine cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    complex(REAL32), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZHER2
!> See also: [[mfi_her2]], [[her2]].
pure subroutine zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    complex(REAL64), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

