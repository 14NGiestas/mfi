module f77_blas_rotmg
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for ROTMG.
!> Supports s, d.
!> See also: [[mfi_rotmg]], [[srotmg]], [[drotmg]].
interface f77_rotmg
!> Original interface for SROTMG
!> See also: [[mfi_rotmg]], [[rotmg]].
pure subroutine srotmg(d1, d2, x1, y1, param)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: y1
    real(REAL32), intent(out) :: param(5)
    real(REAL32), intent(inout) :: d1
    real(REAL32), intent(inout) :: d2
    real(REAL32), intent(inout) :: x1
end subroutine
!> Original interface for DROTMG
!> See also: [[mfi_rotmg]], [[rotmg]].
pure subroutine drotmg(d1, d2, x1, y1, param)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: y1
    real(REAL64), intent(out) :: param(5)
    real(REAL64), intent(inout) :: d1
    real(REAL64), intent(inout) :: d2
    real(REAL64), intent(inout) :: x1
end subroutine
end interface
end module

