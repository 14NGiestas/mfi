module f77_blas_asum
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for ASUM.
!> Supports s, d, sc, dz.
!> See also: [[mfi_asum]], [[sasum]], [[dasum]], [[scasum]], [[dzasum]].
interface f77_asum
!> Original interface for SASUM
!> See also: [[mfi_asum]], [[asum]].
pure function sasum(n, x, incx)
    import :: REAL32
    real(REAL32) :: sasum
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for DASUM
!> See also: [[mfi_asum]], [[asum]].
pure function dasum(n, x, incx)
    import :: REAL64
    real(REAL64) :: dasum
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for SCASUM
!> See also: [[mfi_asum]], [[asum]].
pure function scasum(n, x, incx)
    import :: REAL32
    real(REAL32) :: scasum
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for DZASUM
!> See also: [[mfi_asum]], [[asum]].
pure function dzasum(n, x, incx)
    import :: REAL64
    real(REAL64) :: dzasum
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
end interface
end module

