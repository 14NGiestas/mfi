module f77_blas_dotu
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for DOTU.
!> Supports c, z.
!> See also: [[mfi_dotu]], [[cdotu]], [[zdotu]].
interface f77_dotu
!> Original interface for CDOTU
!> See also: [[mfi_dotu]], [[dotu]].
pure function cdotu(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32) :: cdotu
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
!> Original interface for ZDOTU
!> See also: [[mfi_dotu]], [[dotu]].
pure function zdotu(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64) :: zdotu
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
end interface
end module

