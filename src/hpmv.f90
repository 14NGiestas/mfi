module f77_blas_hpmv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for HPMV.
!> Supports c, z.
!> See also: [[mfi_hpmv]], [[chpmv]], [[zhpmv]].
interface f77_hpmv
!> Original interface for CHPMV
!> See also: [[mfi_hpmv]], [[hpmv]].
pure subroutine chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: ap(*)
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZHPMV
!> See also: [[mfi_hpmv]], [[hpmv]].
pure subroutine zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: ap(*)
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

