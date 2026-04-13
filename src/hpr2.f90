module f77_blas_hpr2
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for HPR2.
!> Supports c, z.
!> See also: [[mfi_hpr2]], [[chpr2]], [[zhpr2]].
interface f77_hpr2
!> Original interface for CHPR2
!> See also: [[mfi_hpr2]], [[hpr2]].
pure subroutine chpr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    complex(REAL32), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZHPR2
!> See also: [[mfi_hpr2]], [[hpr2]].
pure subroutine zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    complex(REAL64), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

