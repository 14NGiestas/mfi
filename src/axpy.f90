module f77_blas_axpy
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for AXPY.
!> Supports s, d, c, z.
!> See also: [[mfi_axpy]], [[saxpy]], [[daxpy]], [[caxpy]], [[zaxpy]].
interface f77_axpy
!> Original interface for SAXPY
!> See also: [[mfi_axpy]], [[axpy]].
pure subroutine saxpy(n, a, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: a
    real(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DAXPY
!> See also: [[mfi_axpy]], [[axpy]].
pure subroutine daxpy(n, a, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: a
    real(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for CAXPY
!> See also: [[mfi_axpy]], [[axpy]].
pure subroutine caxpy(n, a, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: a
    complex(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZAXPY
!> See also: [[mfi_axpy]], [[axpy]].
pure subroutine zaxpy(n, a, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: a
    complex(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
end module

