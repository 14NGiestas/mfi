module f77_blas_scal
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for SCAL.
!> Supports s, d, c, z, cs, zd.
!> See also: [[mfi_scal]], [[sscal]], [[dscal]], [[cscal]], [[zscal]], [[csscal]], [[zdscal]].
interface f77_scal
!> Original interface for SSCAL
!> See also: [[mfi_scal]], [[scal]].
!> SSCAL scales a vector by a constant.
pure subroutine sscal(n, a, x, incx)
    import :: REAL32
    real(REAL32), intent(inout) :: x(*)
    real(REAL32), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for DSCAL
!> See also: [[mfi_scal]], [[scal]].
!> DSCAL scales a vector by a constant.
pure subroutine dscal(n, a, x, incx)
    import :: REAL64
    real(REAL64), intent(inout) :: x(*)
    real(REAL64), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for CSCAL
!> See also: [[mfi_scal]], [[scal]].
!> CSCAL scales a vector by a constant.
pure subroutine cscal(n, a, x, incx)
    import :: REAL32
    complex(REAL32), intent(inout) :: x(*)
    complex(REAL32), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for ZSCAL
!> See also: [[mfi_scal]], [[scal]].
!> ZSCAL scales a vector by a constant.
pure subroutine zscal(n, a, x, incx)
    import :: REAL64
    complex(REAL64), intent(inout) :: x(*)
    complex(REAL64), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for CSSCAL
!> See also: [[mfi_scal]], [[scal]].
!> CSSCAL scales a vector by a constant.
pure subroutine csscal(n, a, x, incx)
    import :: REAL32
    complex(REAL32), intent(inout) :: x(*)
    real(REAL32), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for ZDSCAL
!> See also: [[mfi_scal]], [[scal]].
!> ZDSCAL scales a vector by a constant.
pure subroutine zdscal(n, a, x, incx)
    import :: REAL64
    complex(REAL64), intent(inout) :: x(*)
    real(REAL64), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
end module

