module mfi_blas
use iso_fortran_env
use f77_blas
use f77_blas, only: mfi_rotmg => f77_rotmg
implicit none

! BLAS level 1
!$:mfi_interface('?asum',  DEFAULT_TYPES)
interface mfi_axpy
    module procedure mfi_saxpy
    module procedure mfi_daxpy
    module procedure mfi_caxpy
    module procedure mfi_zaxpy
end interface
interface mfi_copy
    module procedure mfi_scopy
    module procedure mfi_dcopy
    module procedure mfi_ccopy
    module procedure mfi_zcopy
end interface
!$:mfi_interface('?dot',   REAL_TYPES)
!$:mfi_interface('sdsdot', REAL_TYPES)
interface mfi_dotu
    module procedure mfi_cdotu
    module procedure mfi_zdotu
end interface
interface mfi_dotc
    module procedure mfi_cdotc
    module procedure mfi_zdotc
end interface
!$:mfi_interface('?nrm2',  DEFAULT_TYPES)
!$:mfi_interface('?rot',   DEFAULT_TYPES)
!$:mfi_interface('?rotg',  DEFAULT_TYPES)
interface mfi_rotm
    module procedure mfi_srotm
    module procedure mfi_drotm
end interface
!$:mfi_interface('?rotmg', REAL_TYPES)
!$:f77_interface('?scal')
interface mfi_swap
    module procedure mfi_sswap
    module procedure mfi_dswap
    module procedure mfi_cswap
    module procedure mfi_zswap
end interface

! BLAS level 2
interface mfi_gbmv
    module procedure mfi_sgbmv
    module procedure mfi_dgbmv
    module procedure mfi_cgbmv
    module procedure mfi_zgbmv
end interface
interface mfi_gemv
    module procedure mfi_sgemv
    module procedure mfi_dgemv
    module procedure mfi_cgemv
    module procedure mfi_zgemv
end interface
interface mfi_ger
    module procedure mfi_sger
    module procedure mfi_dger
end interface
interface mfi_gerc
    module procedure mfi_cgerc
    module procedure mfi_zgerc
end interface
interface mfi_geru
    module procedure mfi_cgeru
    module procedure mfi_zgeru
end interface
interface mfi_hbmv
    module procedure mfi_chbmv
    module procedure mfi_zhbmv
end interface
interface mfi_hemv
    module procedure mfi_chemv
    module procedure mfi_zhemv
end interface
interface mfi_her
    module procedure mfi_cher
    module procedure mfi_zher
end interface
interface mfi_her2
    module procedure mfi_cher2
    module procedure mfi_zher2
end interface
interface mfi_hpmv
    module procedure mfi_chpmv
    module procedure mfi_zhpmv
end interface
interface mfi_hpr
    module procedure mfi_chpr
    module procedure mfi_zhpr
end interface
interface mfi_hpr2
    module procedure mfi_chpr2
    module procedure mfi_zhpr2
end interface
interface mfi_sbmv
    module procedure mfi_ssbmv
    module procedure mfi_dsbmv
end interface
interface mfi_spmv
    module procedure mfi_sspmv
    module procedure mfi_dspmv
end interface
interface mfi_spr
    module procedure mfi_sspr
    module procedure mfi_dspr
end interface
interface mfi_spr2
    module procedure mfi_sspr2
    module procedure mfi_dspr2
end interface
interface mfi_symv
    module procedure mfi_ssymv
    module procedure mfi_dsymv
end interface
interface mfi_syr
    module procedure mfi_ssyr
    module procedure mfi_dsyr
end interface
interface mfi_syr2
    module procedure mfi_ssyr2
    module procedure mfi_dsyr2
end interface
interface mfi_tbmv
    module procedure mfi_stbmv
    module procedure mfi_dtbmv
    module procedure mfi_ctbmv
    module procedure mfi_ztbmv
end interface
interface mfi_tbsv
    module procedure mfi_stbsv
    module procedure mfi_dtbsv
    module procedure mfi_ctbsv
    module procedure mfi_ztbsv
end interface
interface mfi_tpmv
    module procedure mfi_stpmv
    module procedure mfi_dtpmv
    module procedure mfi_ctpmv
    module procedure mfi_ztpmv
end interface
interface mfi_tpsv
    module procedure mfi_stpsv
    module procedure mfi_dtpsv
    module procedure mfi_ctpsv
    module procedure mfi_ztpsv
end interface
interface mfi_trmv
    module procedure mfi_strmv
    module procedure mfi_dtrmv
    module procedure mfi_ctrmv
    module procedure mfi_ztrmv
end interface
interface mfi_trsv
    module procedure mfi_strsv
    module procedure mfi_dtrsv
    module procedure mfi_ctrsv
    module procedure mfi_ztrsv
end interface

! BLAS level 3
interface mfi_gemm
    module procedure mfi_sgemm
    module procedure mfi_dgemm
    module procedure mfi_cgemm
    module procedure mfi_zgemm
end interface
interface mfi_hemm
    module procedure mfi_chemm
    module procedure mfi_zhemm
end interface
interface mfi_herk
    module procedure mfi_cherk
    module procedure mfi_zherk
end interface
interface mfi_her2k
    module procedure mfi_cher2k
    module procedure mfi_zher2k
end interface
interface mfi_symm
    module procedure mfi_ssymm
    module procedure mfi_dsymm
end interface
interface mfi_syrk
    module procedure mfi_ssyrk
    module procedure mfi_dsyrk
end interface
interface mfi_syr2k
    module procedure mfi_ssyr2k
    module procedure mfi_dsyr2k
end interface
interface mfi_trmm
    module procedure mfi_strmm
    module procedure mfi_dtrmm
    module procedure mfi_ctrmm
    module procedure mfi_ztrmm
end interface
interface mfi_trsm
    module procedure mfi_strsm
    module procedure mfi_dtrsm
    module procedure mfi_ctrsm
    module procedure mfi_ztrsm
end interface

! Extensions
! BLAS level 1 - Utils / Extensions
interface mfi_iamax
    module procedure mfi_isamax
    module procedure mfi_idamax
    module procedure mfi_icamax
    module procedure mfi_izamax
end interface
interface mfi_iamin
    module procedure mfi_isamin
    module procedure mfi_idamin
    module procedure mfi_icamin
    module procedure mfi_izamin
end interface

contains

! BLAS level 1
!$:mfi_interface('?asum',  DEFAULT_TYPES)
pure subroutine mfi_saxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in), optional :: a
    real(wp) :: local_a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
    if (present(a)) then
        local_a = a
    else
        local_a = 1.0_wp
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
    N = size(X)
    call f77_axpy(n,local_a,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_daxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in), optional :: a
    real(wp) :: local_a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
    if (present(a)) then
        local_a = a
    else
        local_a = 1.0_wp
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
    N = size(X)
    call f77_axpy(n,local_a,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_caxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    complex(wp), intent(in), optional :: a
    complex(wp) :: local_a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
    if (present(a)) then
        local_a = a
    else
        local_a = 1.0_wp
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
    N = size(X)
    call f77_axpy(n,local_a,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_zaxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    complex(wp), intent(in), optional :: a
    complex(wp) :: local_a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
    if (present(a)) then
        local_a = a
    else
        local_a = 1.0_wp
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
    N = size(X)
    call f77_axpy(n,local_a,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_scopy(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_copy(n,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_dcopy(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_copy(n,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_ccopy(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_copy(n,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_zcopy(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_copy(n,x,local_incx,y,local_incy)
end subroutine
!$:mfi_interface('?dot',   REAL_TYPES)
!$:mfi_interface('sdsdot', REAL_TYPES)
pure function mfi_cdotu(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp) :: mfi_cdotu
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    integer :: n
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
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
    N = size(X)
    mfi_cdotu = f77_dotu(n,x,local_incx,y,local_incy)
end function
pure function mfi_zdotu(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp) :: mfi_zdotu
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    integer :: n
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
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
    N = size(X)
    mfi_zdotu = f77_dotu(n,x,local_incx,y,local_incy)
end function
pure function mfi_cdotc(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp) :: mfi_cdotc
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    integer :: n
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
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
    N = size(X)
    mfi_cdotc = f77_dotc(n,x,local_incx,y,local_incy)
end function
pure function mfi_zdotc(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp) :: mfi_zdotc
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    integer :: n
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
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
    N = size(X)
    mfi_zdotc = f77_dotc(n,x,local_incx,y,local_incy)
end function
!$:mfi_interface('?nrm2',  DEFAULT_TYPES)
!$:mfi_interface('?rot',   DEFAULT_TYPES)
!$:mfi_interface('?rotg',  DEFAULT_TYPES)
pure subroutine mfi_srotm(x, y, param, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: param(5)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_rotm(n,x,local_incx,y,local_incy,param)
end subroutine
pure subroutine mfi_drotm(x, y, param, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: param(5)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_rotm(n,x,local_incx,y,local_incy,param)
end subroutine
!$:mfi_implement('?rotmg', REAL_TYPES, rotmg)
!$:f77_interface('?scal')
pure subroutine mfi_sswap(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_swap(n,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_dswap(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_swap(n,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_cswap(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_swap(n,x,local_incx,y,local_incy)
end subroutine
pure subroutine mfi_zswap(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    N = size(X)
    call f77_swap(n,x,local_incx,y,local_incy)
end subroutine

! BLAS level 2
pure subroutine mfi_sgbmv(a, x, y, kl, m, alpha, beta, trans, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
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
pure subroutine mfi_dgbmv(a, x, y, kl, m, alpha, beta, trans, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
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
pure subroutine mfi_cgbmv(a, x, y, kl, m, alpha, beta, trans, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
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
pure subroutine mfi_zgbmv(a, x, y, kl, m, alpha, beta, trans, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
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
pure subroutine mfi_sgemv(a, x, y, trans, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
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
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_dgemv(a, x, y, trans, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
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
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_cgemv(a, x, y, trans, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
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
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_zgemv(a, x, y, trans, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
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
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_sger(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: y(:)
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
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
pure subroutine mfi_dger(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: y(:)
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
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
pure subroutine mfi_cgerc(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
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
    call f77_gerc(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
pure subroutine mfi_zgerc(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
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
    call f77_gerc(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
pure subroutine mfi_cgeru(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
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
    call f77_geru(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
pure subroutine mfi_zgeru(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
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
    call f77_geru(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
pure subroutine mfi_chbmv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, k, lda
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
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_hbmv(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_zhbmv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, k, lda
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
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_hbmv(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_chemv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_hemv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_zhemv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_hemv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_cher(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: a(:,:)
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
pure subroutine mfi_zher(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: a(:,:)
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
pure subroutine mfi_cher2(a, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    complex(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
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
pure subroutine mfi_zher2(a, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    complex(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
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
pure subroutine mfi_chpmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: ap(:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    n = size(x)
    call f77_hpmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_zhpmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: ap(:)
    complex(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    n = size(x)
    call f77_hpmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_chpr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
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
    n = size(x)
    call f77_hpr(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
pure subroutine mfi_zhpr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
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
    n = size(x)
    call f77_hpr(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
pure subroutine mfi_chpr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    complex(wp), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    n = size(x)
    call f77_hpr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
pure subroutine mfi_zhpr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    complex(wp), intent(in) :: y(:)
    complex(wp), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    n = size(x)
    call f77_hpr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
pure subroutine mfi_ssbmv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, k, lda
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
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_sbmv(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_dsbmv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n, k, lda
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
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_sbmv(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_sspmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: ap(:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    n = size(x)
    call f77_spmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_dspmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: ap(:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    n = size(x)
    call f77_spmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_sspr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
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
    n = size(x)
    call f77_spr(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
pure subroutine mfi_dspr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
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
    n = size(x)
    call f77_spr(local_uplo,n,local_alpha,x,local_incx,ap)
end subroutine
pure subroutine mfi_sspr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: y(:)
    real(wp), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    n = size(x)
    call f77_spr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
pure subroutine mfi_dspr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: y(:)
    real(wp), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: n
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
    n = size(x)
    call f77_spr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
end subroutine
pure subroutine mfi_ssymv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_symv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_dsymv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
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
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_symv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
pure subroutine mfi_ssyr(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: a(:,:)
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
    call f77_syr(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
pure subroutine mfi_dsyr(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(inout) :: a(:,:)
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
    call f77_syr(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
pure subroutine mfi_ssyr2(a, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: y(:)
    real(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
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
    call f77_syr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
pure subroutine mfi_dsyr2(a, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: y(:)
    real(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
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
    call f77_syr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
pure subroutine mfi_stbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_dtbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_ctbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_ztbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_stbsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_dtbsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_ctbsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_ztbsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, k, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    k = size(a,1)-1
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_tbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_stpmv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: ap(:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
pure subroutine mfi_dtpmv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: ap(:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
pure subroutine mfi_ctpmv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: ap(:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
pure subroutine mfi_ztpmv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: ap(:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
pure subroutine mfi_stpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: ap(:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
pure subroutine mfi_dtpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: ap(:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
pure subroutine mfi_ctpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: ap(:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
pure subroutine mfi_ztpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: ap(:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call f77_tpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
end subroutine
pure subroutine mfi_strmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_dtrmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_ctrmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_ztrmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_strsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_dtrsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_ctrsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine
pure subroutine mfi_ztrsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: x(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call f77_trsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
end subroutine

! BLAS level 3
pure subroutine mfi_sgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: b(:,:)
    real(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: m, n, k, lda, ldb, ldc
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(transb)) then
        local_transb = transb
    else
        local_transb = 'N'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_dgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: b(:,:)
    real(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: m, n, k, lda, ldb, ldc
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(transb)) then
        local_transb = transb
    else
        local_transb = 'N'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_cgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: b(:,:)
    complex(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer :: m, n, k, lda, ldb, ldc
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(transb)) then
        local_transb = transb
    else
        local_transb = 'N'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_zgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: b(:,:)
    complex(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer :: m, n, k, lda, ldb, ldc
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(transb)) then
        local_transb = transb
    else
        local_transb = 'N'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_chemm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: b(:,:)
    complex(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer :: m, n, lda, ldb, ldc
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call f77_hemm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_zhemm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: b(:,:)
    complex(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    complex(wp), intent(in), optional :: beta
    complex(wp) :: local_beta
    integer :: m, n, lda, ldb, ldc
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call f77_hemm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_cherk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: n, k, lda, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call f77_herk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
pure subroutine mfi_zherk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: n, k, lda, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call f77_herk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
pure subroutine mfi_cher2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: b(:,:)
    complex(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: n, k, lda, ldb, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call f77_her2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_zher2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: b(:,:)
    complex(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: n, k, lda, ldb, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call f77_her2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_ssymm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: b(:,:)
    real(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: m, n, lda, ldb, ldc
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call f77_symm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_dsymm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: b(:,:)
    real(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: m, n, lda, ldb, ldc
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call f77_symm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_ssyrk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: n, k, lda, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call f77_syrk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
pure subroutine mfi_dsyrk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: n, k, lda, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldc = max(1,size(c,1))
    call f77_syrk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
end subroutine
pure subroutine mfi_ssyr2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: b(:,:)
    real(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: n, k, lda, ldb, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call f77_syr2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_dsyr2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(in) :: b(:,:)
    real(wp), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    real(wp), intent(in), optional :: beta
    real(wp) :: local_beta
    integer :: n, k, lda, ldb, ldc
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
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
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    n = size(c,2)
    if (local_trans == 'N' .or. local_trans == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    call f77_syr2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
pure subroutine mfi_strmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call f77_trmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
pure subroutine mfi_dtrmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call f77_trmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
pure subroutine mfi_ctrmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call f77_trmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
pure subroutine mfi_ztrmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call f77_trmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
pure subroutine mfi_strsm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call f77_trsm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
pure subroutine mfi_dtrsm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(wp), intent(in), optional :: alpha
    real(wp) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call f77_trsm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
pure subroutine mfi_ctrsm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call f77_trsm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
pure subroutine mfi_ztrsm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(wp), intent(in), optional :: alpha
    complex(wp) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    call f77_trsm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine

! Extensions
! BLAS level 1 - Utils / Extensions
pure function mfi_isamax(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_isamax
    real(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_isamax = f77_iamax(n,x,local_incx)
end function
pure function mfi_idamax(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_idamax
    real(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_idamax = f77_iamax(n,x,local_incx)
end function
pure function mfi_icamax(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_icamax
    complex(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_icamax = f77_iamax(n,x,local_incx)
end function
pure function mfi_izamax(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_izamax
    complex(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_izamax = f77_iamax(n,x,local_incx)
end function
pure function mfi_isamin(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_isamin
    real(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_isamin = f77_iamin(n,x,local_incx)
end function
pure function mfi_idamin(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_idamin
    real(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_idamin = f77_iamin(n,x,local_incx)
end function
pure function mfi_icamin(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_icamin
    complex(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_icamin = f77_iamin(n,x,local_incx)
end function
pure function mfi_izamin(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_izamin
    complex(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_izamin = f77_iamin(n,x,local_incx)
end function

end module
