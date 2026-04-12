module f77_blas_gbmv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for GBMV.
!> Supports s, d, c, z.
!> See also: [[mfi_gbmv]], [[sgbmv]], [[dgbmv]], [[cgbmv]], [[zgbmv]].
interface f77_gbmv
!> Original interface for SGBMV
!> See also: [[mfi_gbmv]], [[gbmv]].
pure subroutine sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DGBMV
!> See also: [[mfi_gbmv]], [[gbmv]].
pure subroutine dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for CGBMV
!> See also: [[mfi_gbmv]], [[gbmv]].
pure subroutine cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZGBMV
!> See also: [[mfi_gbmv]], [[gbmv]].
pure subroutine zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

