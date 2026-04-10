module f77_blas_rot
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for ROT.
!> Supports s, d, c, z, cs, zd.
!> See also: [[mfi_rot]], [[srot]], [[drot]], [[crot]], [[zrot]], [[csrot]], [[zdrot]].
interface f77_rot
!> Original interface for SROT
!> See also: [[mfi_rot]], [[rot]].
!> SROT applies a plane rotation.
pure subroutine srot(n, x, incx, y, incy, c, s)
    import :: REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL32), intent(in) :: c
    real(REAL32), intent(in) :: s
end subroutine
!> Original interface for DROT
!> See also: [[mfi_rot]], [[rot]].
!> DROT applies a plane rotation.
pure subroutine drot(n, x, incx, y, incy, c, s)
    import :: REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL64), intent(in) :: c
    real(REAL64), intent(in) :: s
end subroutine
!> Original interface for CROT
!> See also: [[mfi_rot]], [[rot]].
!> CROT applies a plane rotation.
pure subroutine crot(n, x, incx, y, incy, c, s)
    import :: REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL32), intent(in) :: c
    complex(REAL32), intent(in) :: s
end subroutine
!> Original interface for ZROT
!> See also: [[mfi_rot]], [[rot]].
!> ZROT applies a plane rotation.
pure subroutine zrot(n, x, incx, y, incy, c, s)
    import :: REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL64), intent(in) :: c
    complex(REAL64), intent(in) :: s
end subroutine
!> Original interface for CSROT
!> See also: [[mfi_rot]], [[rot]].
!> CSROT applies a plane rotation.
pure subroutine csrot(n, x, incx, y, incy, c, s)
    import :: REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL32), intent(in) :: c
    real(REAL32), intent(in) :: s
end subroutine
!> Original interface for ZDROT
!> See also: [[mfi_rot]], [[rot]].
!> ZDROT applies a plane rotation.
pure subroutine zdrot(n, x, incx, y, incy, c, s)
    import :: REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL64), intent(in) :: c
    real(REAL64), intent(in) :: s
end subroutine
end interface
end module

