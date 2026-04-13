module f77_blas_hbmv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for HBMV.
!> Supports c, z.
!> See also: [[mfi_hbmv]], [[chbmv]], [[zhbmv]].
interface f77_hbmv
!> Original interface for CHBMV
!> See also: [[mfi_hbmv]], [[hbmv]].
pure subroutine chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZHBMV
!> See also: [[mfi_hbmv]], [[hbmv]].
pure subroutine zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

