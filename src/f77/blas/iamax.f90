module f77_blas_iamax
#if defined(MFI_EXTENSIONS)
    use iso_fortran_env
    use iso_c_binding
    implicit none
#if defined(MFI_LINK_EXTERNAL)
!> Generic old style interface for IAMAX.
!> Supports s, d, c, z.
!> See also: [[mfi_iamax]], [[isamax]], [[idamax]], [[icamax]], [[izamax]].
interface f77_iamax
!> Original interface for ISAMAX
!> See also: [[mfi_iamax]], [[iamax]].
pure function isamax(n, x, incx)
    import :: REAL32
    integer :: isamax
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for IDAMAX
!> See also: [[mfi_iamax]], [[iamax]].
pure function idamax(n, x, incx)
    import :: REAL64
    integer :: idamax
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for ICAMAX
!> See also: [[mfi_iamax]], [[iamax]].
pure function icamax(n, x, incx)
    import :: REAL32
    integer :: icamax
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for IZAMAX
!> See also: [[mfi_iamax]], [[iamax]].
pure function izamax(n, x, incx)
    import :: REAL64
    integer :: izamax
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
end interface
#else
interface f77_iamax
    procedure :: isamax
    procedure :: idamax
    procedure :: icamax
    procedure :: izamax
end interface
contains
pure function isamax(n, x, incx)
    integer :: isamax
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        isamax = 0
        return
    end if
    isamax = minloc(x(1:n:incx),dim=1)
end function
pure function idamax(n, x, incx)
    integer :: idamax
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        idamax = 0
        return
    end if
    idamax = minloc(x(1:n:incx),dim=1)
end function
pure function icamax(n, x, incx)
    integer :: icamax
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        icamax = 0
        return
    end if
    icamax = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
end function
pure function izamax(n, x, incx)
    integer :: izamax
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        izamax = 0
        return
    end if
    izamax = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
end function
#endif
#endif
end module

