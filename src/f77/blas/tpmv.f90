module f77_blas_tpmv
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for TPMV.
!> Supports s, d, c, z.
!> See also: [[mfi_tpmv]], [[stpmv]], [[dtpmv]], [[ctpmv]], [[ztpmv]].
interface f77_tpmv
!> Original interface for STPMV
!> See also: [[mfi_tpmv]], [[tpmv]].
pure subroutine stpmv(uplo, trans, diag, n, ap, x, incx)
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
!> Original interface for DTPMV
!> See also: [[mfi_tpmv]], [[tpmv]].
pure subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
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
!> Original interface for CTPMV
!> See also: [[mfi_tpmv]], [[tpmv]].
pure subroutine ctpmv(uplo, trans, diag, n, ap, x, incx)
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
!> Original interface for ZTPMV
!> See also: [[mfi_tpmv]], [[tpmv]].
pure subroutine ztpmv(uplo, trans, diag, n, ap, x, incx)
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

