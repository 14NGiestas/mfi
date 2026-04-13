module f77_blas_swap
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SWAP.
!> Supports s, d, c, z.
!> See also: [[mfi_swap]], [[sswap]], [[dswap]], [[cswap]], [[zswap]].
interface f77_swap
!> Original interface for SSWAP
!> See also: [[mfi_swap]], [[swap]].
pure subroutine sswap(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DSWAP
!> See also: [[mfi_swap]], [[swap]].
pure subroutine dswap(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for CSWAP
!> See also: [[mfi_swap]], [[swap]].
pure subroutine cswap(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZSWAP
!> See also: [[mfi_swap]], [[swap]].
pure subroutine zswap(n, x, incx, y, incy)
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

