!> Modern fortran interfaces for BLAS
module mfi_blas
use iso_fortran_env
use f77_blas
use f77_blas, only: mfi_rotg  => f77_rotg
use f77_blas, only: mfi_rotmg => f77_rotmg
implicit none

!> Generic modern interface for COPY.
!> Supports s, d, c, z.
!> See also:
!> [[f77_copy:scopy]], [[f77_copy:dcopy]], [[f77_copy:ccopy]], [[f77_copy:zcopy]].
interface mfi_copy
    module procedure :: mfi_scopy
    module procedure :: mfi_dcopy
    module procedure :: mfi_ccopy
    module procedure :: mfi_zcopy
end interface
!> Generic modern interface for SWAP.
!> Supports s, d, c, z.
!> See also:
!> [[f77_swap:sswap]], [[f77_swap:dswap]], [[f77_swap:cswap]], [[f77_swap:zswap]].
interface mfi_swap
    module procedure :: mfi_sswap
    module procedure :: mfi_dswap
    module procedure :: mfi_cswap
    module procedure :: mfi_zswap
end interface
!> Generic modern interface for AXPY.
!> Supports s, d, c, z.
!> See also:
!> [[f77_axpy:saxpy]], [[f77_axpy:daxpy]], [[f77_axpy:caxpy]], [[f77_axpy:zaxpy]].
interface mfi_axpy
    module procedure :: mfi_saxpy
    module procedure :: mfi_daxpy
    module procedure :: mfi_caxpy
    module procedure :: mfi_zaxpy
end interface
!> Generic modern interface for DOT.
!> Supports s, d.
!> See also:
!> [[f77_dot:sdot]], [[f77_dot:ddot]].
interface mfi_dot
    module procedure :: mfi_sdot
    module procedure :: mfi_ddot
end interface
!> Generic modern interface for DOTC.
!> Supports c, z.
!> See also:
!> [[f77_dotc:cdotc]], [[f77_dotc:zdotc]].
interface mfi_dotc
    module procedure :: mfi_cdotc
    module procedure :: mfi_zdotc
end interface
!> Generic modern interface for DOTU.
!> Supports c, z.
!> See also:
!> [[f77_dotu:cdotu]], [[f77_dotu:zdotu]].
interface mfi_dotu
    module procedure :: mfi_cdotu
    module procedure :: mfi_zdotu
end interface
!> Generic modern interface for ASUM.
!> Supports s, d, sc, dz.
!> See also:
!> [[f77_asum:sasum]], [[f77_asum:dasum]], [[f77_asum:scasum]], [[f77_asum:dzasum]].
interface mfi_asum
    module procedure :: mfi_sasum
    module procedure :: mfi_dasum
    module procedure :: mfi_scasum
    module procedure :: mfi_dzasum
end interface
!> Generic modern interface for NRM2.
!> Supports s, d, sc, dz.
!> See also:
!> [[f77_nrm2:snrm2]], [[f77_nrm2:dnrm2]], [[f77_nrm2:scnrm2]], [[f77_nrm2:dznrm2]].
interface mfi_nrm2
    module procedure :: mfi_snrm2
    module procedure :: mfi_dnrm2
    module procedure :: mfi_scnrm2
    module procedure :: mfi_dznrm2
end interface
!> Generic modern interface for ROT.
!> Supports s, d, c, z, cs, zd.
!> See also:
!> [[f77_rot:srot]], [[f77_rot:drot]], [[f77_rot:crot]], [[f77_rot:zrot]], [[f77_rot:csrot]], [[f77_rot:zdrot]].
interface mfi_rot
    module procedure :: mfi_srot
    module procedure :: mfi_drot
    module procedure :: mfi_crot
    module procedure :: mfi_zrot
    module procedure :: mfi_csrot
    module procedure :: mfi_zdrot
end interface
!> Generic modern interface for ROTM.
!> Supports s, d.
!> See also:
!> [[f77_rotm:srotm]], [[f77_rotm:drotm]].
interface mfi_rotm
    module procedure :: mfi_srotm
    module procedure :: mfi_drotm
end interface
!> Generic modern interface for SCAL.
!> Supports s, d, c, z, cs, zd.
!> See also:
!> [[f77_scal:sscal]], [[f77_scal:dscal]], [[f77_scal:cscal]], [[f77_scal:zscal]], [[f77_scal:csscal]], [[f77_scal:zdscal]].
interface mfi_scal
    module procedure :: mfi_sscal
    module procedure :: mfi_dscal
    module procedure :: mfi_cscal
    module procedure :: mfi_zscal
    module procedure :: mfi_csscal
    module procedure :: mfi_zdscal
end interface
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
!> Generic modern interface for GEMV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_gemv:sgemv]], [[f77_gemv:dgemv]], [[f77_gemv:cgemv]], [[f77_gemv:zgemv]].
interface mfi_gemv
    module procedure :: mfi_sgemv
    module procedure :: mfi_dgemv
    module procedure :: mfi_cgemv
    module procedure :: mfi_zgemv
end interface
!> Generic modern interface for GER.
!> Supports s, d.
!> See also:
!> [[f77_ger:sger]], [[f77_ger:dger]].
interface mfi_ger
    module procedure :: mfi_sger
    module procedure :: mfi_dger
end interface
!> Generic modern interface for GERC.
!> Supports c, z.
!> See also:
!> [[f77_gerc:cgerc]], [[f77_gerc:zgerc]].
interface mfi_gerc
    module procedure :: mfi_cgerc
    module procedure :: mfi_zgerc
end interface
!> Generic modern interface for GERU.
!> Supports c, z.
!> See also:
!> [[f77_geru:cgeru]], [[f77_geru:zgeru]].
interface mfi_geru
    module procedure :: mfi_cgeru
    module procedure :: mfi_zgeru
end interface
!> Generic modern interface for HBMV.
!> Supports c, z.
!> See also:
!> [[f77_hbmv:chbmv]], [[f77_hbmv:zhbmv]].
interface mfi_hbmv
    module procedure :: mfi_chbmv
    module procedure :: mfi_zhbmv
end interface
!> Generic modern interface for HEMV.
!> Supports c, z.
!> See also:
!> [[f77_hemv:chemv]], [[f77_hemv:zhemv]].
interface mfi_hemv
    module procedure :: mfi_chemv
    module procedure :: mfi_zhemv
end interface
!> Generic modern interface for HER.
!> Supports c, z.
!> See also:
!> [[f77_her:cher]], [[f77_her:zher]].
interface mfi_her
    module procedure :: mfi_cher
    module procedure :: mfi_zher
end interface
!> Generic modern interface for HER2.
!> Supports c, z.
!> See also:
!> [[f77_her2:cher2]], [[f77_her2:zher2]].
interface mfi_her2
    module procedure :: mfi_cher2
    module procedure :: mfi_zher2
end interface
!> Generic modern interface for HPMV.
!> Supports c, z.
!> See also:
!> [[f77_hpmv:chpmv]], [[f77_hpmv:zhpmv]].
interface mfi_hpmv
    module procedure :: mfi_chpmv
    module procedure :: mfi_zhpmv
end interface
!> Generic modern interface for HPR.
!> Supports c, z.
!> See also:
!> [[f77_hpr:chpr]], [[f77_hpr:zhpr]].
interface mfi_hpr
    module procedure :: mfi_chpr
    module procedure :: mfi_zhpr
end interface
!> Generic modern interface for HPR2.
!> Supports c, z.
!> See also:
!> [[f77_hpr2:chpr2]], [[f77_hpr2:zhpr2]].
interface mfi_hpr2
    module procedure :: mfi_chpr2
    module procedure :: mfi_zhpr2
end interface
!> Generic modern interface for SBMV.
!> Supports s, d.
!> See also:
!> [[f77_sbmv:ssbmv]], [[f77_sbmv:dsbmv]].
interface mfi_sbmv
    module procedure :: mfi_ssbmv
    module procedure :: mfi_dsbmv
end interface
!> Generic modern interface for SPMV.
!> Supports s, d.
!> See also:
!> [[f77_spmv:sspmv]], [[f77_spmv:dspmv]].
interface mfi_spmv
    module procedure :: mfi_sspmv
    module procedure :: mfi_dspmv
end interface
!> Generic modern interface for SPR.
!> Supports s, d.
!> See also:
!> [[f77_spr:sspr]], [[f77_spr:dspr]].
interface mfi_spr
    module procedure :: mfi_sspr
    module procedure :: mfi_dspr
end interface
!> Generic modern interface for SPR2.
!> Supports s, d.
!> See also:
!> [[f77_spr2:sspr2]], [[f77_spr2:dspr2]].
interface mfi_spr2
    module procedure :: mfi_sspr2
    module procedure :: mfi_dspr2
end interface
!> Generic modern interface for SYMV.
!> Supports s, d.
!> See also:
!> [[f77_symv:ssymv]], [[f77_symv:dsymv]].
interface mfi_symv
    module procedure :: mfi_ssymv
    module procedure :: mfi_dsymv
end interface
!> Generic modern interface for SYR.
!> Supports s, d.
!> See also:
!> [[f77_syr:ssyr]], [[f77_syr:dsyr]].
interface mfi_syr
    module procedure :: mfi_ssyr
    module procedure :: mfi_dsyr
end interface
!> Generic modern interface for SYR2.
!> Supports s, d.
!> See also:
!> [[f77_syr2:ssyr2]], [[f77_syr2:dsyr2]].
interface mfi_syr2
    module procedure :: mfi_ssyr2
    module procedure :: mfi_dsyr2
end interface
!> Generic modern interface for TBMV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_tbmv:stbmv]], [[f77_tbmv:dtbmv]], [[f77_tbmv:ctbmv]], [[f77_tbmv:ztbmv]].
interface mfi_tbmv
    module procedure :: mfi_stbmv
    module procedure :: mfi_dtbmv
    module procedure :: mfi_ctbmv
    module procedure :: mfi_ztbmv
end interface
!> Generic modern interface for TBSV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_tbsv:stbsv]], [[f77_tbsv:dtbsv]], [[f77_tbsv:ctbsv]], [[f77_tbsv:ztbsv]].
interface mfi_tbsv
    module procedure :: mfi_stbsv
    module procedure :: mfi_dtbsv
    module procedure :: mfi_ctbsv
    module procedure :: mfi_ztbsv
end interface
!> Generic modern interface for TPMV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_tpmv:stpmv]], [[f77_tpmv:dtpmv]], [[f77_tpmv:ctpmv]], [[f77_tpmv:ztpmv]].
interface mfi_tpmv
    module procedure :: mfi_stpmv
    module procedure :: mfi_dtpmv
    module procedure :: mfi_ctpmv
    module procedure :: mfi_ztpmv
end interface
!> Generic modern interface for TPSV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_tpsv:stpsv]], [[f77_tpsv:dtpsv]], [[f77_tpsv:ctpsv]], [[f77_tpsv:ztpsv]].
interface mfi_tpsv
    module procedure :: mfi_stpsv
    module procedure :: mfi_dtpsv
    module procedure :: mfi_ctpsv
    module procedure :: mfi_ztpsv
end interface
!> Generic modern interface for TRMV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_trmv:strmv]], [[f77_trmv:dtrmv]], [[f77_trmv:ctrmv]], [[f77_trmv:ztrmv]].
interface mfi_trmv
    module procedure :: mfi_strmv
    module procedure :: mfi_dtrmv
    module procedure :: mfi_ctrmv
    module procedure :: mfi_ztrmv
end interface
!> Generic modern interface for TRSV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_trsv:strsv]], [[f77_trsv:dtrsv]], [[f77_trsv:ctrsv]], [[f77_trsv:ztrsv]].
interface mfi_trsv
    module procedure :: mfi_strsv
    module procedure :: mfi_dtrsv
    module procedure :: mfi_ctrsv
    module procedure :: mfi_ztrsv
end interface
!> Generic modern interface for GEMM.
!> Supports s, d, c, z.
!> See also:
!> [[f77_gemm:sgemm]], [[f77_gemm:dgemm]], [[f77_gemm:cgemm]], [[f77_gemm:zgemm]].
interface mfi_gemm
    module procedure :: mfi_sgemm
    module procedure :: mfi_dgemm
    module procedure :: mfi_cgemm
    module procedure :: mfi_zgemm
end interface
!> Generic modern interface for HEMM.
!> Supports c, z.
!> See also:
!> [[f77_hemm:chemm]], [[f77_hemm:zhemm]].
interface mfi_hemm
    module procedure :: mfi_chemm
    module procedure :: mfi_zhemm
end interface
!> Generic modern interface for HERK.
!> Supports c, z.
!> See also:
!> [[f77_herk:cherk]], [[f77_herk:zherk]].
interface mfi_herk
    module procedure :: mfi_cherk
    module procedure :: mfi_zherk
end interface
!> Generic modern interface for HER2K.
!> Supports c, z.
!> See also:
!> [[f77_her2k:cher2k]], [[f77_her2k:zher2k]].
interface mfi_her2k
    module procedure :: mfi_cher2k
    module procedure :: mfi_zher2k
end interface
!> Generic modern interface for SYMM.
!> Supports s, d.
!> See also:
!> [[f77_symm:ssymm]], [[f77_symm:dsymm]].
interface mfi_symm
    module procedure :: mfi_ssymm
    module procedure :: mfi_dsymm
end interface
!> Generic modern interface for SYRK.
!> Supports s, d.
!> See also:
!> [[f77_syrk:ssyrk]], [[f77_syrk:dsyrk]].
interface mfi_syrk
    module procedure :: mfi_ssyrk
    module procedure :: mfi_dsyrk
end interface
!> Generic modern interface for SYR2K.
!> Supports s, d.
!> See also:
!> [[f77_syr2k:ssyr2k]], [[f77_syr2k:dsyr2k]].
interface mfi_syr2k
    module procedure :: mfi_ssyr2k
    module procedure :: mfi_dsyr2k
end interface
!> Generic modern interface for TRMM.
!> Supports s, d, c, z.
!> See also:
!> [[f77_trmm:strmm]], [[f77_trmm:dtrmm]], [[f77_trmm:ctrmm]], [[f77_trmm:ztrmm]].
interface mfi_trmm
    module procedure :: mfi_strmm
    module procedure :: mfi_dtrmm
    module procedure :: mfi_ctrmm
    module procedure :: mfi_ztrmm
end interface
!> Generic modern interface for TRSM.
!> Supports s, d, c, z.
!> See also:
!> [[f77_trsm:strsm]], [[f77_trsm:dtrsm]], [[f77_trsm:ctrsm]], [[f77_trsm:ztrsm]].
interface mfi_trsm
    module procedure :: mfi_strsm
    module procedure :: mfi_dtrsm
    module procedure :: mfi_ctrsm
    module procedure :: mfi_ztrsm
end interface
!> Generic modern interface for LAMCH.
!> Supports s, d.
!> See also:
!> [[f77_lamch:slamch]], [[f77_lamch:dlamch]].
interface mfi_lamch
    module procedure :: mfi_slamch
    module procedure :: mfi_dlamch
end interface

! Extensions
! BLAS level 1 - Utils / Extensions
!> Generic modern interface for IAMAX.
!> Supports s, d, c, z.
!> See also:
!> [[f77_iamax:isamax]], [[f77_iamax:idamax]], [[f77_iamax:icamax]], [[f77_iamax:izamax]].
interface mfi_iamax
    module procedure :: mfi_isamax
    module procedure :: mfi_idamax
    module procedure :: mfi_icamax
    module procedure :: mfi_izamax
end interface
!> Generic modern interface for IAMIN.
!> Supports s, d, c, z.
!> See also:
!> [[f77_iamin:isamin]], [[f77_iamin:idamin]], [[f77_iamin:icamin]], [[f77_iamin:izamin]].
interface mfi_iamin
    module procedure :: mfi_isamin
    module procedure :: mfi_idamin
    module procedure :: mfi_icamin
    module procedure :: mfi_izamin
end interface

contains


!> Modern interface for [[f77_copy:scopy]].
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine mfi_scopy(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: y(:)
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
!> Modern interface for [[f77_copy:dcopy]].
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine mfi_dcopy(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: y(:)
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
!> Modern interface for [[f77_copy:ccopy]].
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine mfi_ccopy(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
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
!> Modern interface for [[f77_copy:zcopy]].
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine mfi_zcopy(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
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
!> Modern interface for [[f77_swap:sswap]].
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine mfi_sswap(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: y(:)
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
!> Modern interface for [[f77_swap:dswap]].
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine mfi_dswap(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: y(:)
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
!> Modern interface for [[f77_swap:cswap]].
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine mfi_cswap(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
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
!> Modern interface for [[f77_swap:zswap]].
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine mfi_zswap(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
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
!> Modern interface for [[f77_axpy:saxpy]].
!> See also: [[mfi_axpy]], [[f77_axpy]].
pure subroutine mfi_saxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in), optional :: a
    real(REAL32) :: local_a
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
!> Modern interface for [[f77_axpy:daxpy]].
!> See also: [[mfi_axpy]], [[f77_axpy]].
pure subroutine mfi_daxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in), optional :: a
    real(REAL64) :: local_a
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
!> Modern interface for [[f77_axpy:caxpy]].
!> See also: [[mfi_axpy]], [[f77_axpy]].
pure subroutine mfi_caxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
    complex(REAL32), intent(in), optional :: a
    complex(REAL32) :: local_a
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
!> Modern interface for [[f77_axpy:zaxpy]].
!> See also: [[mfi_axpy]], [[f77_axpy]].
pure subroutine mfi_zaxpy(x, y, a, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
    complex(REAL64), intent(in), optional :: a
    complex(REAL64) :: local_a
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
!> Modern interface for [[f77_dot:sdot]].
!> See also: [[mfi_dot]], [[f77_dot]].
pure function mfi_sdot(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32) :: mfi_sdot
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: y(:)
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
    mfi_sdot = sdot(n,x,local_incx,y,local_incy)
end function
!> Modern interface for [[f77_dot:ddot]].
!> See also: [[mfi_dot]], [[f77_dot]].
pure function mfi_ddot(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64) :: mfi_ddot
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: y(:)
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
    mfi_ddot = ddot(n,x,local_incx,y,local_incy)
end function
!> Modern interface for [[f77_dotc:cdotc]].
!> See also: [[mfi_dotc]], [[f77_dotc]].
pure function mfi_cdotc(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32) :: mfi_cdotc
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: y(:)
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
!> Modern interface for [[f77_dotc:zdotc]].
!> See also: [[mfi_dotc]], [[f77_dotc]].
pure function mfi_zdotc(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64) :: mfi_zdotc
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: y(:)
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
!> Modern interface for [[f77_dotu:cdotu]].
!> See also: [[mfi_dotu]], [[f77_dotu]].
pure function mfi_cdotu(x, y, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32) :: mfi_cdotu
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: y(:)
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
!> Modern interface for [[f77_dotu:zdotu]].
!> See also: [[mfi_dotu]], [[f77_dotu]].
pure function mfi_zdotu(x, y, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64) :: mfi_zdotu
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: y(:)
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
!> Modern interface for [[f77_asum:sasum]].
!> See also: [[mfi_asum]], [[f77_asum]].
pure function mfi_sasum(x, incx)
    real(REAL32) :: mfi_sasum
    real(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_sasum = sasum(n, x, local_incx)
end function
!> Modern interface for [[f77_asum:dasum]].
!> See also: [[mfi_asum]], [[f77_asum]].
pure function mfi_dasum(x, incx)
    real(REAL64) :: mfi_dasum
    real(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dasum = dasum(n, x, local_incx)
end function
!> Modern interface for [[f77_asum:scasum]].
!> See also: [[mfi_asum]], [[f77_asum]].
pure function mfi_scasum(x, incx)
    real(REAL32) :: mfi_scasum
    complex(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_scasum = scasum(n, x, local_incx)
end function
!> Modern interface for [[f77_asum:dzasum]].
!> See also: [[mfi_asum]], [[f77_asum]].
pure function mfi_dzasum(x, incx)
    real(REAL64) :: mfi_dzasum
    complex(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dzasum = dzasum(n, x, local_incx)
end function
!> Modern interface for [[f77_nrm2:snrm2]].
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function mfi_snrm2(x, incx)
    real(REAL32) :: mfi_snrm2
    real(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_snrm2 = snrm2(n, x, local_incx)
end function
!> Modern interface for [[f77_nrm2:dnrm2]].
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function mfi_dnrm2(x, incx)
    real(REAL64) :: mfi_dnrm2
    real(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dnrm2 = dnrm2(n, x, local_incx)
end function
!> Modern interface for [[f77_nrm2:scnrm2]].
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function mfi_scnrm2(x, incx)
    real(REAL32) :: mfi_scnrm2
    complex(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_scnrm2 = scnrm2(n, x, local_incx)
end function
!> Modern interface for [[f77_nrm2:dznrm2]].
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function mfi_dznrm2(x, incx)
    real(REAL64) :: mfi_dznrm2
    complex(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_dznrm2 = dznrm2(n, x, local_incx)
end function
!> Modern interface for [[f77_rot:srot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - s*xi
!>```
pure subroutine mfi_srot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: x(:)
    real(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in) :: c
    real(REAL32), intent(in) :: s
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
!> Modern interface for [[f77_rot:drot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - s*xi
!>```
pure subroutine mfi_drot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: x(:)
    real(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in) :: c
    real(REAL64), intent(in) :: s
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
!> Modern interface for [[f77_rot:crot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_crot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in) :: c
    complex(REAL32), intent(in) :: s
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
!> Modern interface for [[f77_rot:zrot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_zrot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in) :: c
    complex(REAL64), intent(in) :: s
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
!> Modern interface for [[f77_rot:csrot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_csrot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: x(:)
    complex(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in) :: c
    real(REAL32), intent(in) :: s
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
!> Modern interface for [[f77_rot:zdrot]].
!> See also: [[mfi_rot]], [[f77_rot]].
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
!>```
pure subroutine mfi_zdrot(x, y, c, s, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: x(:)
    complex(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in) :: c
    real(REAL64), intent(in) :: s
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
!> Modern interface for [[f77_rotm:srotm]].
!> See also: [[mfi_rotm]], [[f77_rotm]].
pure subroutine mfi_srotm(x, y, param, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: x(:)
    real(REAL32), intent(inout) :: y(:)
    real(REAL32), intent(in) :: param(5)
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
!> Modern interface for [[f77_rotm:drotm]].
!> See also: [[mfi_rotm]], [[f77_rotm]].
pure subroutine mfi_drotm(x, y, param, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: x(:)
    real(REAL64), intent(inout) :: y(:)
    real(REAL64), intent(in) :: param(5)
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
!> Modern interface for [[f77_scal:sscal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_SSCAL scales a vector by a constant.
pure subroutine mfi_sscal(a, x, incx)
    real(REAL32), intent(inout) :: x(:)
    real(REAL32), intent(in) :: a
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
!> Modern interface for [[f77_scal:dscal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_DSCAL scales a vector by a constant.
pure subroutine mfi_dscal(a, x, incx)
    real(REAL64), intent(inout) :: x(:)
    real(REAL64), intent(in) :: a
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
!> Modern interface for [[f77_scal:cscal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_CSCAL scales a vector by a constant.
pure subroutine mfi_cscal(a, x, incx)
    complex(REAL32), intent(inout) :: x(:)
    complex(REAL32), intent(in) :: a
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
!> Modern interface for [[f77_scal:zscal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_ZSCAL scales a vector by a constant.
pure subroutine mfi_zscal(a, x, incx)
    complex(REAL64), intent(inout) :: x(:)
    complex(REAL64), intent(in) :: a
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
!> Modern interface for [[f77_scal:csscal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_CSSCAL scales a vector by a constant.
pure subroutine mfi_csscal(a, x, incx)
    complex(REAL32), intent(inout) :: x(:)
    real(REAL32), intent(in) :: a
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
!> Modern interface for [[f77_scal:zdscal]].
!> See also: [[mfi_scal]], [[f77_scal]].
!> MFI_ZDSCAL scales a vector by a constant.
pure subroutine mfi_zdscal(a, x, incx)
    complex(REAL64), intent(inout) :: x(:)
    real(REAL64), intent(in) :: a
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
!> Modern interface for [[f77_gbmv:sgbmv]].
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
    call sgbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gbmv:dgbmv]].
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
    call dgbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gbmv:cgbmv]].
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
    call cgbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gbmv:zgbmv]].
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
    call zgbmv(local_trans,local_m,n,local_kl,ku,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gemv:sgemv]].
!> See also: [[mfi_gemv]], [[f77_gemv]].
pure subroutine mfi_sgemv(a, x, y, trans, alpha, beta, incx, incy)
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
!> Modern interface for [[f77_gemv:dgemv]].
!> See also: [[mfi_gemv]], [[f77_gemv]].
pure subroutine mfi_dgemv(a, x, y, trans, alpha, beta, incx, incy)
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
!> Modern interface for [[f77_gemv:cgemv]].
!> See also: [[mfi_gemv]], [[f77_gemv]].
pure subroutine mfi_cgemv(a, x, y, trans, alpha, beta, incx, incy)
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
!> Modern interface for [[f77_gemv:zgemv]].
!> See also: [[mfi_gemv]], [[f77_gemv]].
pure subroutine mfi_zgemv(a, x, y, trans, alpha, beta, incx, incy)
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
!> Modern interface for [[f77_ger:sger]].
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
    call sger(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
!> Modern interface for [[f77_ger:dger]].
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
    call dger(m,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
!> Modern interface for [[f77_gerc:cgerc]].
!> See also: [[mfi_gerc]], [[f77_gerc]].
pure subroutine mfi_cgerc(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: y(:)
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
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
!> Modern interface for [[f77_gerc:zgerc]].
!> See also: [[mfi_gerc]], [[f77_gerc]].
pure subroutine mfi_zgerc(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: y(:)
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
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
!> Modern interface for [[f77_geru:cgeru]].
!> See also: [[mfi_geru]], [[f77_geru]].
pure subroutine mfi_cgeru(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: y(:)
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
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
!> Modern interface for [[f77_geru:zgeru]].
!> See also: [[mfi_geru]], [[f77_geru]].
pure subroutine mfi_zgeru(a, x, y, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: y(:)
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
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
!> Modern interface for [[f77_hbmv:chbmv]].
!> See also: [[mfi_hbmv]], [[f77_hbmv]].
pure subroutine mfi_chbmv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
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
!> Modern interface for [[f77_hbmv:zhbmv]].
!> See also: [[mfi_hbmv]], [[f77_hbmv]].
pure subroutine mfi_zhbmv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
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
!> Modern interface for [[f77_hemv:chemv]].
!> See also: [[mfi_hemv]], [[f77_hemv]].
pure subroutine mfi_chemv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
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
!> Modern interface for [[f77_hemv:zhemv]].
!> See also: [[mfi_hemv]], [[f77_hemv]].
pure subroutine mfi_zhemv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
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
!> Modern interface for [[f77_her:cher]].
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
    call cher(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
!> Modern interface for [[f77_her:zher]].
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
    call zher(local_uplo,n,local_alpha,x,local_incx,a,lda)
end subroutine
!> Modern interface for [[f77_her2:cher2]].
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
    call cher2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
!> Modern interface for [[f77_her2:zher2]].
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
    call zher2(local_uplo,n,local_alpha,x,local_incx,y,local_incy,a,lda)
end subroutine
!> Modern interface for [[f77_hpmv:chpmv]].
!> See also: [[mfi_hpmv]], [[f77_hpmv]].
pure subroutine mfi_chpmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: ap(:)
    complex(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
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
!> Modern interface for [[f77_hpmv:zhpmv]].
!> See also: [[mfi_hpmv]], [[f77_hpmv]].
pure subroutine mfi_zhpmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: ap(:)
    complex(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
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
!> Modern interface for [[f77_hpr:chpr]].
!> See also: [[mfi_hpr]], [[f77_hpr]].
pure subroutine mfi_chpr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(inout) :: ap(:)
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
!> Modern interface for [[f77_hpr:zhpr]].
!> See also: [[mfi_hpr]], [[f77_hpr]].
pure subroutine mfi_zhpr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(inout) :: ap(:)
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
!> Modern interface for [[f77_hpr2:chpr2]].
!> See also: [[mfi_hpr2]], [[f77_hpr2]].
pure subroutine mfi_chpr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(:)
    complex(REAL32), intent(in) :: y(:)
    complex(REAL32), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
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
!> Modern interface for [[f77_hpr2:zhpr2]].
!> See also: [[mfi_hpr2]], [[f77_hpr2]].
pure subroutine mfi_zhpr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(:)
    complex(REAL64), intent(in) :: y(:)
    complex(REAL64), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
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
!> Modern interface for [[f77_sbmv:ssbmv]].
!> See also: [[mfi_sbmv]], [[f77_sbmv]].
pure subroutine mfi_ssbmv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
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
!> Modern interface for [[f77_sbmv:dsbmv]].
!> See also: [[mfi_sbmv]], [[f77_sbmv]].
pure subroutine mfi_dsbmv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
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
!> Modern interface for [[f77_spmv:sspmv]].
!> See also: [[mfi_spmv]], [[f77_spmv]].
pure subroutine mfi_sspmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: ap(:)
    real(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
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
!> Modern interface for [[f77_spmv:dspmv]].
!> See also: [[mfi_spmv]], [[f77_spmv]].
pure subroutine mfi_dspmv(ap, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: ap(:)
    real(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
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
!> Modern interface for [[f77_spr:sspr]].
!> See also: [[mfi_spr]], [[f77_spr]].
pure subroutine mfi_sspr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
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
!> Modern interface for [[f77_spr:dspr]].
!> See also: [[mfi_spr]], [[f77_spr]].
pure subroutine mfi_dspr(ap, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
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
!> Modern interface for [[f77_spr2:sspr2]].
!> See also: [[mfi_spr2]], [[f77_spr2]].
pure subroutine mfi_sspr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: y(:)
    real(REAL32), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
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
!> Modern interface for [[f77_spr2:dspr2]].
!> See also: [[mfi_spr2]], [[f77_spr2]].
pure subroutine mfi_dspr2(ap, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: y(:)
    real(REAL64), intent(inout) :: ap(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
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
!> Modern interface for [[f77_symv:ssymv]].
!> See also: [[mfi_symv]], [[f77_symv]].
pure subroutine mfi_ssymv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
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
!> Modern interface for [[f77_symv:dsymv]].
!> See also: [[mfi_symv]], [[f77_symv]].
pure subroutine mfi_dsymv(a, x, y, uplo, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: y(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
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
!> Modern interface for [[f77_syr:ssyr]].
!> See also: [[mfi_syr]], [[f77_syr]].
pure subroutine mfi_ssyr(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
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
!> Modern interface for [[f77_syr:dsyr]].
!> See also: [[mfi_syr]], [[f77_syr]].
pure subroutine mfi_dsyr(a, x, uplo, alpha, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
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
!> Modern interface for [[f77_syr2:ssyr2]].
!> See also: [[mfi_syr2]], [[f77_syr2]].
pure subroutine mfi_ssyr2(a, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(:)
    real(REAL32), intent(in) :: y(:)
    real(REAL32), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
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
!> Modern interface for [[f77_syr2:dsyr2]].
!> See also: [[mfi_syr2]], [[f77_syr2]].
pure subroutine mfi_dsyr2(a, x, y, uplo, alpha, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(:)
    real(REAL64), intent(in) :: y(:)
    real(REAL64), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
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
!> Modern interface for [[f77_tbmv:stbmv]].
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
pure subroutine mfi_stbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tbmv:dtbmv]].
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
pure subroutine mfi_dtbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tbmv:ctbmv]].
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
pure subroutine mfi_ctbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tbmv:ztbmv]].
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
pure subroutine mfi_ztbmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tbsv:stbsv]].
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
pure subroutine mfi_stbsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tbsv:dtbsv]].
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
pure subroutine mfi_dtbsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tbsv:ctbsv]].
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
pure subroutine mfi_ctbsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tbsv:ztbsv]].
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
pure subroutine mfi_ztbsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tpmv:stpmv]].
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
pure subroutine mfi_stpmv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: ap(:)
    real(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tpmv:dtpmv]].
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
pure subroutine mfi_dtpmv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: ap(:)
    real(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tpmv:ctpmv]].
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
pure subroutine mfi_ctpmv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: ap(:)
    complex(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tpmv:ztpmv]].
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
pure subroutine mfi_ztpmv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: ap(:)
    complex(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tpsv:stpsv]].
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
pure subroutine mfi_stpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: ap(:)
    real(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tpsv:dtpsv]].
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
pure subroutine mfi_dtpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: ap(:)
    real(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tpsv:ctpsv]].
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
pure subroutine mfi_ctpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: ap(:)
    complex(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_tpsv:ztpsv]].
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
pure subroutine mfi_ztpsv(ap, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: ap(:)
    complex(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_trmv:strmv]].
!> See also: [[mfi_trmv]], [[f77_trmv]].
pure subroutine mfi_strmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_trmv:dtrmv]].
!> See also: [[mfi_trmv]], [[f77_trmv]].
pure subroutine mfi_dtrmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_trmv:ctrmv]].
!> See also: [[mfi_trmv]], [[f77_trmv]].
pure subroutine mfi_ctrmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_trmv:ztrmv]].
!> See also: [[mfi_trmv]], [[f77_trmv]].
pure subroutine mfi_ztrmv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_trsv:strsv]].
!> See also: [[mfi_trsv]], [[f77_trsv]].
pure subroutine mfi_strsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_trsv:dtrsv]].
!> See also: [[mfi_trsv]], [[f77_trsv]].
pure subroutine mfi_dtrsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_trsv:ctrsv]].
!> See also: [[mfi_trsv]], [[f77_trsv]].
pure subroutine mfi_ctrsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: x(:)
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
!> Modern interface for [[f77_trsv:ztrsv]].
!> See also: [[mfi_trsv]], [[f77_trsv]].
pure subroutine mfi_ztrsv(a, x, uplo, trans, diag, incx)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: x(:)
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
!> Modern interface for [[f77_gemm:sgemm]].
!> See also: [[mfi_gemm]], [[f77_gemm]].
pure subroutine mfi_sgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(in) :: b(:,:)
    real(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
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
!> Modern interface for [[f77_gemm:dgemm]].
!> See also: [[mfi_gemm]], [[f77_gemm]].
pure subroutine mfi_dgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(in) :: b(:,:)
    real(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
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
!> Modern interface for [[f77_gemm:cgemm]].
!> See also: [[mfi_gemm]], [[f77_gemm]].
pure subroutine mfi_cgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(in) :: b(:,:)
    complex(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
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
!> Modern interface for [[f77_gemm:zgemm]].
!> See also: [[mfi_gemm]], [[f77_gemm]].
pure subroutine mfi_zgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(in) :: b(:,:)
    complex(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
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
!> Modern interface for [[f77_hemm:chemm]].
!> See also: [[mfi_hemm]], [[f77_hemm]].
pure subroutine mfi_chemm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(in) :: b(:,:)
    complex(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
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
!> Modern interface for [[f77_hemm:zhemm]].
!> See also: [[mfi_hemm]], [[f77_hemm]].
pure subroutine mfi_zhemm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(in) :: b(:,:)
    complex(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
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
!> Modern interface for [[f77_herk:cherk]].
!> See also: [[mfi_herk]], [[f77_herk]].
pure subroutine mfi_cherk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: c(:,:)
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
!> Modern interface for [[f77_herk:zherk]].
!> See also: [[mfi_herk]], [[f77_herk]].
pure subroutine mfi_zherk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: c(:,:)
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
!> Modern interface for [[f77_her2k:cher2k]].
!> See also: [[mfi_her2k]], [[f77_her2k]].
pure subroutine mfi_cher2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(in) :: b(:,:)
    complex(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
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
!> Modern interface for [[f77_her2k:zher2k]].
!> See also: [[mfi_her2k]], [[f77_her2k]].
pure subroutine mfi_zher2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(in) :: b(:,:)
    complex(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
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
!> Modern interface for [[f77_symm:ssymm]].
!> See also: [[mfi_symm]], [[f77_symm]].
pure subroutine mfi_ssymm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(in) :: b(:,:)
    real(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
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
!> Modern interface for [[f77_symm:dsymm]].
!> See also: [[mfi_symm]], [[f77_symm]].
pure subroutine mfi_dsymm(a, b, c, side, uplo, alpha, beta)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(in) :: b(:,:)
    real(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
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
!> Modern interface for [[f77_syrk:ssyrk]].
!> See also: [[mfi_syrk]], [[f77_syrk]].
pure subroutine mfi_ssyrk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
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
!> Modern interface for [[f77_syrk:dsyrk]].
!> See also: [[mfi_syrk]], [[f77_syrk]].
pure subroutine mfi_dsyrk(a, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
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
!> Modern interface for [[f77_syr2k:ssyr2k]].
!> See also: [[mfi_syr2k]], [[f77_syr2k]].
pure subroutine mfi_ssyr2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(in) :: b(:,:)
    real(REAL32), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
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
!> Modern interface for [[f77_syr2k:dsyr2k]].
!> See also: [[mfi_syr2k]], [[f77_syr2k]].
pure subroutine mfi_dsyr2k(a, b, c, uplo, trans, alpha, beta)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(in) :: b(:,:)
    real(REAL64), intent(inout) :: c(:,:)
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: uplo
    character :: local_uplo
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
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
!> Modern interface for [[f77_trmm:strmm]].
!> See also: [[mfi_trmm]], [[f77_trmm]].
pure subroutine mfi_strmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
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
!> Modern interface for [[f77_trmm:dtrmm]].
!> See also: [[mfi_trmm]], [[f77_trmm]].
pure subroutine mfi_dtrmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
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
!> Modern interface for [[f77_trmm:ctrmm]].
!> See also: [[mfi_trmm]], [[f77_trmm]].
pure subroutine mfi_ctrmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
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
!> Modern interface for [[f77_trmm:ztrmm]].
!> See also: [[mfi_trmm]], [[f77_trmm]].
pure subroutine mfi_ztrmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
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
!> Modern interface for [[f77_trsm:strsm]].
!> See also: [[mfi_trsm]], [[f77_trsm]].
pure subroutine mfi_strsm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
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
!> Modern interface for [[f77_trsm:dtrsm]].
!> See also: [[mfi_trsm]], [[f77_trsm]].
pure subroutine mfi_dtrsm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
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
!> Modern interface for [[f77_trsm:ctrsm]].
!> See also: [[mfi_trsm]], [[f77_trsm]].
pure subroutine mfi_ctrsm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
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
!> Modern interface for [[f77_trsm:ztrsm]].
!> See also: [[mfi_trsm]], [[f77_trsm]].
pure subroutine mfi_ztrsm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
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
!> Modern interface for [[f77_lamch:slamch]].
!> See also: [[mfi_lamch]], [[f77_lamch]].
pure function mfi_slamch(cmach, kind) result(res)
    integer, parameter :: wp = REAL32
    character, intent(in) :: cmach
    real(REAL32), intent(in) :: kind
    !! Just a kind placeholder
    real(REAL32) :: res
    res = slamch(cmach)
end function
!> Modern interface for [[f77_lamch:dlamch]].
!> See also: [[mfi_lamch]], [[f77_lamch]].
pure function mfi_dlamch(cmach, kind) result(res)
    integer, parameter :: wp = REAL64
    character, intent(in) :: cmach
    real(REAL64), intent(in) :: kind
    !! Just a kind placeholder
    real(REAL64) :: res
    res = dlamch(cmach)
end function

! Extensions
! BLAS level 1 - Utils / Extensions
!> Modern interface for [[f77_iamax:isamax]].
!> See also: [[mfi_iamax]], [[f77_iamax]].
pure function mfi_isamax(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_isamax
    real(REAL32), intent(in) :: x(:)
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
!> Modern interface for [[f77_iamax:idamax]].
!> See also: [[mfi_iamax]], [[f77_iamax]].
pure function mfi_idamax(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_idamax
    real(REAL64), intent(in) :: x(:)
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
!> Modern interface for [[f77_iamax:icamax]].
!> See also: [[mfi_iamax]], [[f77_iamax]].
pure function mfi_icamax(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_icamax
    complex(REAL32), intent(in) :: x(:)
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
!> Modern interface for [[f77_iamax:izamax]].
!> See also: [[mfi_iamax]], [[f77_iamax]].
pure function mfi_izamax(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_izamax
    complex(REAL64), intent(in) :: x(:)
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
!> Modern interface for [[f77_iamin:isamin]].
!> See also: [[mfi_iamin]], [[f77_iamin]].
pure function mfi_isamin(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_isamin
    real(REAL32), intent(in) :: x(:)
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
!> Modern interface for [[f77_iamin:idamin]].
!> See also: [[mfi_iamin]], [[f77_iamin]].
pure function mfi_idamin(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_idamin
    real(REAL64), intent(in) :: x(:)
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
!> Modern interface for [[f77_iamin:icamin]].
!> See also: [[mfi_iamin]], [[f77_iamin]].
pure function mfi_icamin(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_icamin
    complex(REAL32), intent(in) :: x(:)
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
!> Modern interface for [[f77_iamin:izamin]].
!> See also: [[mfi_iamin]], [[f77_iamin]].
pure function mfi_izamin(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_izamin
    complex(REAL64), intent(in) :: x(:)
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

end module
