module mfi_blas_lamch
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for LAMCH.
!> Supports s, d.
!> See also:
!> [[f77_lamch:slamch]], [[f77_lamch:dlamch]].
interface mfi_lamch
    module procedure :: mfi_slamch
    module procedure :: mfi_dlamch
end interface

contains

!> Modern interface for [[f77_lamch:f77_lamch]].
!> See also: [[mfi_lamch]], [[f77_lamch]].
pure function mfi_slamch(cmach, kind) result(res)
    integer, parameter :: wp = REAL32
    character, intent(in) :: cmach
    real(REAL32), intent(in) :: kind
    !! Just a kind placeholder
    real(REAL32) :: res
    res = slamch(cmach)
end function
!> Modern interface for [[f77_lamch:f77_lamch]].
!> See also: [[mfi_lamch]], [[f77_lamch]].
pure function mfi_dlamch(cmach, kind) result(res)
    integer, parameter :: wp = REAL64
    character, intent(in) :: cmach
    real(REAL64), intent(in) :: kind
    !! Just a kind placeholder
    real(REAL64) :: res
    res = dlamch(cmach)
end function
end module


!> cuBLAS interfaces and constants
