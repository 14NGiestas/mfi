module mfi_blas_her
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for HER.
!> Supports c, z.
!> See also:
!> [[f77_her:cher]], [[f77_her:zher]].
interface mfi_her
    module procedure :: mfi_cher
    module procedure :: mfi_zher
end interface

contains

!> Modern interface for [[f77_her:f77_her]].
!> See also: [[mfi_her]], [[f77_her]].
pure subroutine mfi_cher(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_her(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
!> Modern interface for [[f77_her:f77_her]].
!> See also: [[mfi_her]], [[f77_her]].
pure subroutine mfi_zher(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_her(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
end module

