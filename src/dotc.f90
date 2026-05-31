module f77_blas_dotc
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for DOTC.
!> Supports c, z.
!> See also: [[mfi_dotc]], [[cdotc]], [[zdotc]].
interface f77_dotc
!> Original interface for CDOTC
!> See also: [[mfi_dotc]], [[dotc]].
pure function cdotc(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32) :: cdotc
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
!> Original interface for ZDOTC
!> See also: [[mfi_dotc]], [[dotc]].
pure function zdotc(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64) :: zdotc
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
end interface
end module

