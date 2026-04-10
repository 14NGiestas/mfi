module mfi_blas_nrm2
    use iso_fortran_env
    use f77_blas
#if defined(MFI_CUBLAS)
    use iso_c_binding
    use mfi_blas_cublas
#endif
#if defined(MFI_EXTENSIONS)
    use mfi_blas_extensions
#endif
    implicit none

!> Generic modern interface for NRM2.
!> Supports s, d, sc, dz.
!> See also:
!> [[f77_nrm2:snrm2]], [[f77_nrm2:dnrm2]], [[f77_nrm2:scnrm2]], [[f77_nrm2:dznrm2]].
interface mfi_nrm2
    module procedure :: mfi_snrm2
    module procedure :: mfi_dnrm2
    module procedure :: mfi_scnrm2
    module procedure :: mfi_dznrm2
end interface

contains

!> Modern interface for [[f77_nrm2:f77_nrm2]].
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function mfi_snrm2(x, incx)
    real(REAL32) :: mfi_snrm2
    real(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_snrm2 = f77_nrm2(n, x, local_incx)
end function
!> Modern interface for [[f77_nrm2:f77_nrm2]].
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function mfi_dnrm2(x, incx)
    real(REAL64) :: mfi_dnrm2
    real(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dnrm2 = f77_nrm2(n, x, local_incx)
end function
!> Modern interface for [[f77_nrm2:f77_nrm2]].
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function mfi_scnrm2(x, incx)
    real(REAL32) :: mfi_scnrm2
    complex(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_scnrm2 = f77_nrm2(n, x, local_incx)
end function
!> Modern interface for [[f77_nrm2:f77_nrm2]].
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function mfi_dznrm2(x, incx)
    real(REAL64) :: mfi_dznrm2
    complex(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dznrm2 = f77_nrm2(n, x, local_incx)
end function
end module

