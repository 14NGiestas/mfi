module f77_blas_iamin
    use iso_fortran_env
    use iso_c_binding
    implicit none
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
end module

!> cuBLAS v2 interfaces
