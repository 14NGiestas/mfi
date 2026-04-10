module mfi_blas_ger
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for GER.
!> Supports s, d.
!> See also:
!> [[f77_ger:sger]], [[f77_ger:dger]].
interface mfi_ger
    module procedure :: mfi_sger
    module procedure :: mfi_dger
end interface

contains

!> Modern interface for [[f77_ger:f77_ger]].
!> See also: [[mfi_ger]], [[f77_ger]].
pure subroutine mfi_sger(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: y(:)
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    if (present(incy)) then
        local_incy = incy
    else
        local_incy = 1
    end if
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call f77_ger(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
!> Modern interface for [[f77_ger:f77_ger]].
!> See also: [[mfi_ger]], [[f77_ger]].
pure subroutine mfi_dger(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: y(:)
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    if (present(incy)) then
        local_incy = incy
    else
        local_incy = 1
    end if
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call f77_ger(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
end module

