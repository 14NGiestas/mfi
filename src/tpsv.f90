module f77_blas_tpsv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for TPSV.
!> Supports s, d, c, z.
!> See also: [[mfi_tpsv]], [[stpsv]], [[dtpsv]], [[ctpsv]], [[ztpsv]].
interface f77_tpsv
!> Original interface for STPSV
!> See also: [[mfi_tpsv]], [[tpsv]].
pure subroutine stpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: ap(*)
    real(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for DTPSV
!> See also: [[mfi_tpsv]], [[tpsv]].
pure subroutine dtpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: ap(*)
    real(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for CTPSV
!> See also: [[mfi_tpsv]], [[tpsv]].
pure subroutine ctpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: ap(*)
    complex(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for ZTPSV
!> See also: [[mfi_tpsv]], [[tpsv]].
pure subroutine ztpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: ap(*)
    complex(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
end module

