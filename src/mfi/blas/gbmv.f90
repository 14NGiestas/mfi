module mfi_blas_gbmv
    use iso_fortran_env
    use f77_blas
    implicit none

!> Generic modern interface for GBMV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_gbmv:sgbmv]], [[f77_gbmv:dgbmv]], [[f77_gbmv:cgbmv]], [[f77_gbmv:zgbmv]].
interface mfi_gbmv
    module procedure :: mfi_sgbmv
    module procedure :: mfi_dgbmv
    module procedure :: mfi_cgbmv
    module procedure :: mfi_zgbmv
end interface

contains

!> Modern interface for [[f77_gbmv:f77_gbmv]].
!> See also: [[mfi_gbmv]], [[f77_gbmv]].
pure subroutine mfi_sgbmv(a, x, y, kl, m, alpha, beta, trans, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
    integer, intent(in), optional :: kl
    integer :: local_kl
    integer, intent(in), optional :: m
    integer :: local_m
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, ku, lda
    n = size(a,2)
    lda = max(1,size(a,1))
    if (present(kl)) then
        local_kl = kl
    else
        local_kl = (lda-1)/2
    end if
    if (present(m)) then
        local_m = m
    else
        local_m = n
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
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
    ku = lda-local_kl-1
    call f77_gbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gbmv:f77_gbmv]].
!> See also: [[mfi_gbmv]], [[f77_gbmv]].
pure subroutine mfi_dgbmv(a, x, y, kl, m, alpha, beta, trans, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
    integer, intent(in), optional :: kl
    integer :: local_kl
    integer, intent(in), optional :: m
    integer :: local_m
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, ku, lda
    n = size(a,2)
    lda = max(1,size(a,1))
    if (present(kl)) then
        local_kl = kl
    else
        local_kl = (lda-1)/2
    end if
    if (present(m)) then
        local_m = m
    else
        local_m = n
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
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
    ku = lda-local_kl-1
    call f77_gbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gbmv:f77_gbmv]].
!> See also: [[mfi_gbmv]], [[f77_gbmv]].
pure subroutine mfi_cgbmv(a, x, y, kl, m, alpha, beta, trans, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
    integer, intent(in), optional :: kl
    integer :: local_kl
    integer, intent(in), optional :: m
    integer :: local_m
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, ku, lda
    n = size(a,2)
    lda = max(1,size(a,1))
    if (present(kl)) then
        local_kl = kl
    else
        local_kl = (lda-1)/2
    end if
    if (present(m)) then
        local_m = m
    else
        local_m = n
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
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
    ku = lda-local_kl-1
    call f77_gbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gbmv:f77_gbmv]].
!> See also: [[mfi_gbmv]], [[f77_gbmv]].
pure subroutine mfi_zgbmv(a, x, y, kl, m, alpha, beta, trans, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
    integer, intent(in), optional :: kl
    integer :: local_kl
    integer, intent(in), optional :: m
    integer :: local_m
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, ku, lda
    n = size(a,2)
    lda = max(1,size(a,1))
    if (present(kl)) then
        local_kl = kl
    else
        local_kl = (lda-1)/2
    end if
    if (present(m)) then
        local_m = m
    else
        local_m = n
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
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
    ku = lda-local_kl-1
    call f77_gbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
end module

