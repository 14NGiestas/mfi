module f77_blas_hpr
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for HPR.
!> Supports c, z.
!> See also: [[mfi_hpr]], [[chpr]], [[zhpr]].
interface f77_hpr
!> Original interface for CHPR
!> See also: [[mfi_hpr]], [[hpr]].
pure subroutine chpr(uplo, n, alpha, x, incx, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for ZHPR
!> See also: [[mfi_hpr]], [[hpr]].
pure subroutine zhpr(uplo, n, alpha, x, incx, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
end module

