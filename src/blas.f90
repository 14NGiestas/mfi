!> Modern fortran interfaces for BLAS
module mfi_blas
use iso_fortran_env
use f77_blas
use f77_blas, only: mfi_rotmg => f77_rotmg
implicit none

! BLAS level 1
interface mfi_asum
    module procedure mfi_sasum
    module procedure mfi_dasum
    module procedure mfi_scasum
    module procedure mfi_dzasum
end interface
interface mfi_nrm2
    module procedure mfi_snrm2
    module procedure mfi_dnrm2
    module procedure mfi_scnrm2
    module procedure mfi_dznrm2
end interface
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
interface mfi_rot
    module procedure mfi_srot
    module procedure mfi_drot
    module procedure mfi_crot
    module procedure mfi_zrot
    module procedure mfi_csrot
    module procedure mfi_zdrot
end interface
!$:mfi_interface('?rotg',  DEFAULT_TYPES)
interface mfi_rotm
    module procedure mfi_srotm
    module procedure mfi_drotm
end interface
!$:mfi_interface('?rotmg', REAL_TYPES)
interface mfi_scal
    module procedure mfi_sscal
    module procedure mfi_dscal
    module procedure mfi_cscal
    module procedure mfi_zscal
    module procedure mfi_csscal
    module procedure mfi_zdscal
end interface
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
interface mfi_lamch
    module procedure mfi_slamch
    module procedure mfi_dlamch
end interface

contains

! BLAS level 1
pure function mfi_snrm2(x, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    real(wp) :: mfi_snrm2
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_snrm2 = snrm2(n, x, local_incx)
end function
pure function mfi_dnrm2(x, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    real(wp) :: mfi_dnrm2
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dnrm2 = dnrm2(n, x, local_incx)
end function
pure function mfi_scnrm2(x, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    real(wp) :: mfi_scnrm2
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_scnrm2 = scnrm2(n, x, local_incx)
end function
pure function mfi_dznrm2(x, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    real(wp) :: mfi_dznrm2
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dznrm2 = dznrm2(n, x, local_incx)
end function
pure function mfi_sasum(x, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    real(wp) :: mfi_sasum
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_sasum = sasum(n, x, local_incx)
end function
pure function mfi_dasum(x, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    real(wp) :: mfi_dasum
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dasum = dasum(n, x, local_incx)
end function
pure function mfi_scasum(x, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    real(wp) :: mfi_scasum
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_scasum = scasum(n, x, local_incx)
end function
pure function mfi_dzasum(x, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    real(wp) :: mfi_dzasum
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dzasum = dzasum(n, x, local_incx)
end function
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
    call saxpy(n,local_a,x,local_incx,y,local_incy)
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
    call daxpy(n,local_a,x,local_incx,y,local_incy)
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
    call caxpy(n,local_a,x,local_incx,y,local_incy)
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
    call zaxpy(n,local_a,x,local_incx,y,local_incy)
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
    call scopy(n,x,local_incx,y,local_incy)
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
    call dcopy(n,x,local_incx,y,local_incy)
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
    call ccopy(n,x,local_incx,y,local_incy)
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
    call zcopy(n,x,local_incx,y,local_incy)
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
    mfi_cdotu = cdotu(n,x,local_incx,y,local_incy)
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
    mfi_zdotu = zdotu(n,x,local_incx,y,local_incy)
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
    mfi_cdotc = cdotc(n,x,local_incx,y,local_incy)
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
    mfi_zdotc = zdotc(n,x,local_incx,y,local_incy)
end function
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - s*xi
!>```
pure subroutine mfi_srot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: c
    real(wp), intent(in) :: s
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
    n = size(x)
    call srot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - s*xi
!>```
pure subroutine mfi_drot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: x(:)
    real(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: c
    real(wp), intent(in) :: s
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
    n = size(x)
    call drot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_crot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: x(:)
    complex(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: c
    complex(wp), intent(in) :: s
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
    n = size(x)
    call crot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_zrot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: x(:)
    complex(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: c
    complex(wp), intent(in) :: s
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
    n = size(x)
    call zrot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_csrot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: x(:)
    complex(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: c
    real(wp), intent(in) :: s
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
    n = size(x)
    call csrot(n,x,local_incx,y,local_incy,c,s)
end subroutine
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_zdrot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: x(:)
    complex(wp), intent(inout) :: y(:)
    real(wp), intent(in) :: c
    real(wp), intent(in) :: s
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
    n = size(x)
    call zdrot(n,x,local_incx,y,local_incy,c,s)
end subroutine
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
    call srotm(n,x,local_incx,y,local_incy,param)
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
    call drotm(n,x,local_incx,y,local_incy,param)
end subroutine
!$:mfi_implement('?rotmg', REAL_TYPES, rotmg)
!> MFI_SSCAL scales a vector by a constant.
pure subroutine mfi_sscal(x, a, incx)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: x(:)
    real(wp), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call sscal(n,a,x,local_incx)
end subroutine
!> MFI_DSCAL scales a vector by a constant.
pure subroutine mfi_dscal(x, a, incx)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: x(:)
    real(wp), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call dscal(n,a,x,local_incx)
end subroutine
!> MFI_CSCAL scales a vector by a constant.
pure subroutine mfi_cscal(x, a, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: x(:)
    complex(wp), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call cscal(n,a,x,local_incx)
end subroutine
!> MFI_ZSCAL scales a vector by a constant.
pure subroutine mfi_zscal(x, a, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: x(:)
    complex(wp), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call zscal(n,a,x,local_incx)
end subroutine
!> MFI_CSSCAL scales a vector by a constant.
pure subroutine mfi_csscal(x, a, incx)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: x(:)
    real(wp), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call csscal(n,a,x,local_incx)
end subroutine
!> MFI_ZDSCAL scales a vector by a constant.
pure subroutine mfi_zdscal(x, a, incx)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: x(:)
    real(wp), intent(in) :: a
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    call zdscal(n,a,x,local_incx)
end subroutine
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
    call sswap(n,x,local_incx,y,local_incy)
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
    call dswap(n,x,local_incx,y,local_incy)
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
    call cswap(n,x,local_incx,y,local_incy)
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
    call zswap(n,x,local_incx,y,local_incy)
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
    call sgbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call dgbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call cgbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call zgbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call sgemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call dgemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call cgemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call zgemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call sger(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call dger(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call cgerc(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call zgerc(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call cgeru(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call zgeru(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call chbmv(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call zhbmv(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call chemv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call zhemv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call cher(local_uplo,n,local_alpha,x,local_incx,a,lda)
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
    call zher(local_uplo,n,local_alpha,x,local_incx,a,lda)
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
    call cher2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call zher2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call chpmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
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
    call zhpmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
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
    call chpr(local_uplo,n,local_alpha,x,local_incx,ap)
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
    call zhpr(local_uplo,n,local_alpha,x,local_incx,ap)
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
    call chpr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
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
    call zhpr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
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
    call ssbmv(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call dsbmv(local_uplo,n,k,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call sspmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
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
    call dspmv(local_uplo,n,local_alpha,ap,x,local_incx,local_beta,y,local_incy)
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
    call sspr(local_uplo,n,local_alpha,x,local_incx,ap)
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
    call dspr(local_uplo,n,local_alpha,x,local_incx,ap)
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
    call sspr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
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
    call dspr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,ap)
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
    call ssymv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call dsymv(local_uplo,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
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
    call ssyr(local_uplo,n,local_alpha,x,local_incx,a,lda)
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
    call dsyr(local_uplo,n,local_alpha,x,local_incx,a,lda)
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
    call ssyr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call dsyr2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
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
    call stbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
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
    call dtbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
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
    call ctbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
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
    call ztbmv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
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
    call stbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
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
    call dtbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
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
    call ctbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
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
    call ztbsv(local_uplo,local_trans,local_diag,n,k,a,lda,x,local_incx)
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
    call stpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
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
    call dtpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
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
    call ctpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
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
    call ztpmv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
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
    call stpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
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
    call dtpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
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
    call ctpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
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
    call ztpsv(local_uplo,local_trans,local_diag,n,ap,x,local_incx)
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
    call strmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
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
    call dtrmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
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
    call ctrmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
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
    call ztrmv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
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
    call strsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
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
    call dtrsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
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
    call ctrsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
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
    call ztrsv(local_uplo,local_trans,local_diag,n,a,lda,x,local_incx)
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
    call sgemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call dgemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call cgemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call zgemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call chemm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call zhemm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call cherk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
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
    call zherk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
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
    call cher2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call zher2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call ssymm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call dsymm(local_side,local_uplo,m,n,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call ssyrk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
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
    call dsyrk(local_uplo,local_trans,n,k,local_alpha,a,lda,local_beta,c,ldc)
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
    call ssyr2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call dsyr2k(local_uplo,local_trans,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
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
    call strmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
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
    call dtrmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
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
    call ctrmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
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
    call ztrmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
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
    call strsm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
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
    call dtrsm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
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
    call ctrsm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
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
    call ztrsm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
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
    mfi_isamax = isamax(n,x,local_incx)
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
    mfi_idamax = idamax(n,x,local_incx)
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
    mfi_icamax = icamax(n,x,local_incx)
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
    mfi_izamax = izamax(n,x,local_incx)
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
    mfi_isamin = isamin(n,x,local_incx)
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
    mfi_idamin = idamin(n,x,local_incx)
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
    mfi_icamin = icamin(n,x,local_incx)
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
    mfi_izamin = izamin(n,x,local_incx)
end function
pure function mfi_slamch(cmach, kind) result(res)
    integer, parameter :: wp = REAL32
    character, intent(in) :: cmach
    real(wp), intent(in) :: kind
    !! Just a kind placeholder
    real(wp) :: res
    res = slamch(cmach)
end function
pure function mfi_dlamch(cmach, kind) result(res)
    integer, parameter :: wp = REAL64
    character, intent(in) :: cmach
    real(wp), intent(in) :: kind
    !! Just a kind placeholder
    real(wp) :: res
    res = dlamch(cmach)
end function

end module
