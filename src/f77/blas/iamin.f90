module f77_blas_iamin
#if defined(MFI_EXTENSIONS)
    use iso_fortran_env
    use iso_c_binding
    implicit none
#if defined(MFI_LINK_EXTERNAL)
!> Generic old style interface for IAMIN.
!> Supports s, d, c, z.
!> See also: [[mfi_iamin]], [[isamin]], [[idamin]], [[icamin]], [[izamin]].
interface f77_iamin
!> Original interface for ISAMIN
!> See also: [[mfi_iamin]], [[iamin]].
pure function isamin(n, x, incx)
    import :: REAL32
    integer :: isamin
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for IDAMIN
!> See also: [[mfi_iamin]], [[iamin]].
pure function idamin(n, x, incx)
    import :: REAL64
    integer :: idamin
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for ICAMIN
!> See also: [[mfi_iamin]], [[iamin]].
pure function icamin(n, x, incx)
    import :: REAL32
    integer :: icamin
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for IZAMIN
!> See also: [[mfi_iamin]], [[iamin]].
pure function izamin(n, x, incx)
    import :: REAL64
    integer :: izamin
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
end interface
#else
interface f77_iamin
    procedure :: isamin
    procedure :: idamin
    procedure :: icamin
    procedure :: izamin
end interface
contains
pure function isamin(n, x, incx)
    integer :: isamin
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        isamin = 0
        return
    end if
    isamin = minloc(x(1:n:incx),dim=1)
end function
pure function idamin(n, x, incx)
    integer :: idamin
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        idamin = 0
        return
    end if
    idamin = minloc(x(1:n:incx),dim=1)
end function
pure function icamin(n, x, incx)
    integer :: icamin
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        icamin = 0
        return
    end if
    icamin = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
end function
pure function izamin(n, x, incx)
    integer :: izamin
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        izamin = 0
        return
    end if
    izamin = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
end function
#endif
#endif
end module

!> cuBLAS v2 interfaces
