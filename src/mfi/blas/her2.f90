module mfi_blas_her2
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for HER2.
!> Supports c, z.
!> See also:
!> [[f77_her2:cher2]], [[f77_her2:zher2]].
interface mfi_her2
    module procedure :: mfi_cher2
    module procedure :: mfi_zher2
end interface

contains

!> Modern interface for [[f77_her2:f77_her2]].
!> See also: [[mfi_her2]], [[f77_her2]].
pure subroutine mfi_cher2(a, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: y(:)
    complex(REAL32), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_her2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
!> Modern interface for [[f77_her2:f77_her2]].
!> See also: [[mfi_her2]], [[f77_her2]].
pure subroutine mfi_zher2(a, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: y(:)
    complex(REAL64), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_her2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
end module

