module f77_blas_copy
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for COPY.
!> Supports s, d, c, z.
!> See also: [[mfi_copy]], [[scopy]], [[dcopy]], [[ccopy]], [[zcopy]].
interface f77_copy
!> Original interface for SCOPY
!> See also: [[mfi_copy]], [[copy]].
pure subroutine scopy(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DCOPY
!> See also: [[mfi_copy]], [[copy]].
pure subroutine dcopy(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for CCOPY
!> See also: [[mfi_copy]], [[copy]].
pure subroutine ccopy(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZCOPY
!> See also: [[mfi_copy]], [[copy]].
pure subroutine zcopy(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

