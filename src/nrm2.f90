module f77_blas_nrm2
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for NRM2.
!> Supports s, d, sc, dz.
!> See also: [[mfi_nrm2]], [[snrm2]], [[dnrm2]], [[scnrm2]], [[dznrm2]].
interface f77_nrm2
!> Original interface for SNRM2
!> See also: [[mfi_nrm2]], [[nrm2]].
pure function snrm2(n, x, incx)
    import :: REAL32
    real(REAL32) :: snrm2
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for DNRM2
!> See also: [[mfi_nrm2]], [[nrm2]].
pure function dnrm2(n, x, incx)
    import :: REAL64
    real(REAL64) :: dnrm2
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for SCNRM2
!> See also: [[mfi_nrm2]], [[nrm2]].
pure function scnrm2(n, x, incx)
    import :: REAL32
    real(REAL32) :: scnrm2
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for DZNRM2
!> See also: [[mfi_nrm2]], [[nrm2]].
pure function dznrm2(n, x, incx)
    import :: REAL64
    real(REAL64) :: dznrm2
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
end interface
end module

