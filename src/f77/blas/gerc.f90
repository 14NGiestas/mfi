module f77_blas_gerc
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for GERC.
!> Supports c, z.
!> See also: [[mfi_gerc]], [[cgerc]], [[zgerc]].
interface f77_gerc
!> Original interface for CGERC
!> See also: [[mfi_gerc]], [[gerc]].
pure subroutine cgerc(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZGERC
!> See also: [[mfi_gerc]], [[gerc]].
pure subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

