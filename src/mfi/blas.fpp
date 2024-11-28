#:mute
#:include "common.fpp"
#:include "src/mfi/blas/lamch.fypp"
#:include "src/mfi/blas/asum.fypp"
#:include "src/mfi/blas/axpy.fypp"
#:include "src/mfi/blas/copy_swap.fypp"
#:include "src/mfi/blas/dot_product.fypp"
#:include "src/mfi/blas/rotm.fypp"
#:include "src/mfi/blas/iamin_iamax.fypp"
#:include "src/mfi/blas/gbmv.fypp"
#:include "src/mfi/blas/gemv.fypp"
#:include "src/mfi/blas/ger_gerc_geru.fypp"
#:include "src/mfi/blas/hbmv_sbmv.fypp"
#:include "src/mfi/blas/hemv_symv.fypp"
#:include "src/mfi/blas/her.fypp"
#:include "src/mfi/blas/syr.fypp"
#:include "src/mfi/blas/her_syr2.fypp"
#:include "src/mfi/blas/hpmv_spmv.fypp"
#:include "src/mfi/blas/hpr.fypp"
#:include "src/mfi/blas/spr.fypp"
#:include "src/mfi/blas/hpr_spr2.fypp"
#:include "src/mfi/blas/tbmv_tbsv.fypp"
#:include "src/mfi/blas/tpmv_tpsv.fypp"
#:include "src/mfi/blas/trmv_trsv.fypp"
#:include "src/mfi/blas/gemm.fypp"
#:include "src/mfi/blas/hemm_symm.fypp"
#:include "src/mfi/blas/herk.fypp"
#:include "src/mfi/blas/syrk.fypp"
#:include "src/mfi/blas/her2k.fypp"
#:include "src/mfi/blas/syr2k.fypp"
#:include "src/mfi/blas/trmm_trsm.fypp"
#:endmute
!> Modern fortran interfaces for BLAS
module mfi_blas
use iso_fortran_env
use f77_blas
use f77_blas, only: mfi_rotmg => f77_rotmg
implicit none

! BLAS level 1
$:mfi_interface('?asum',  DEFAULT_TYPES, &
    f=lambda pfx: 'sc' if pfx == 'c' else &
                  'dz' if pfx == 'z' else pfx)
$:mfi_interface('?axpy',  DEFAULT_TYPES)
$:mfi_interface('?copy',  DEFAULT_TYPES)
!$:mfi_interface('?dot',   REAL_TYPES)
!$:mfi_interface('sdsdot', REAL_TYPES)
$:mfi_interface('?dotu',  COMPLEX_TYPES)
$:mfi_interface('?dotc',  COMPLEX_TYPES)
!$:mfi_interface('?nrm2',  DEFAULT_TYPES)
!$:mfi_interface('?rot',   DEFAULT_TYPES)
!$:mfi_interface('?rotg',  DEFAULT_TYPES)
$:mfi_interface('?rotm',  REAL_TYPES)
!$:mfi_interface('?rotmg', REAL_TYPES)
!$:f77_interface('?scal')
$:mfi_interface('?swap',  DEFAULT_TYPES)

! BLAS level 2
$:mfi_interface('?gbmv',  DEFAULT_TYPES)
$:mfi_interface('?gemv',  DEFAULT_TYPES)
$:mfi_interface('?ger',   REAL_TYPES)
$:mfi_interface('?gerc',  COMPLEX_TYPES)
$:mfi_interface('?geru',  COMPLEX_TYPES)
$:mfi_interface('?hbmv',  COMPLEX_TYPES)
$:mfi_interface('?hemv',  COMPLEX_TYPES)
$:mfi_interface('?her',   COMPLEX_TYPES)
$:mfi_interface('?her2',  COMPLEX_TYPES)
$:mfi_interface('?hpmv',  COMPLEX_TYPES)
$:mfi_interface('?hpr',   COMPLEX_TYPES)
$:mfi_interface('?hpr2',  COMPLEX_TYPES)
$:mfi_interface('?sbmv',  REAL_TYPES)
$:mfi_interface('?spmv',  REAL_TYPES)
$:mfi_interface('?spr',   REAL_TYPES)
$:mfi_interface('?spr2',  REAL_TYPES)
$:mfi_interface('?symv',  REAL_TYPES)
$:mfi_interface('?syr',   REAL_TYPES)
$:mfi_interface('?syr2',  REAL_TYPES)
$:mfi_interface('?tbmv',  DEFAULT_TYPES)
$:mfi_interface('?tbsv',  DEFAULT_TYPES)
$:mfi_interface('?tpmv',  DEFAULT_TYPES)
$:mfi_interface('?tpsv',  DEFAULT_TYPES)
$:mfi_interface('?trmv',  DEFAULT_TYPES)
$:mfi_interface('?trsv',  DEFAULT_TYPES)

! BLAS level 3
$:mfi_interface('?gemm',  DEFAULT_TYPES)
$:mfi_interface('?hemm',  COMPLEX_TYPES)
$:mfi_interface('?herk',  COMPLEX_TYPES)
$:mfi_interface('?her2k', COMPLEX_TYPES)
$:mfi_interface('?symm',  REAL_TYPES)
$:mfi_interface('?syrk',  REAL_TYPES)
$:mfi_interface('?syr2k', REAL_TYPES)
$:mfi_interface('?trmm',  DEFAULT_TYPES)
$:mfi_interface('?trsm',  DEFAULT_TYPES)

! Extensions
! BLAS level 1 - Utils / Extensions
$:mfi_interface('i?amax', DEFAULT_TYPES)
#:if defined('MFI_EXTENSIONS')
$:mfi_interface('i?amin', DEFAULT_TYPES)
#:endif
$:mfi_interface('?lamch', REAL_TYPES)

contains

! BLAS level 1
$:mfi_implement('?asum',  DEFAULT_TYPES, asum, &
    f=lambda pfx: 'sc' if pfx == 'c' else &
                  'dz' if pfx == 'z' else pfx)
$:mfi_implement('?axpy',  DEFAULT_TYPES, axpy)
$:mfi_implement('?copy',  DEFAULT_TYPES, copy_swap)
!$:mfi_interface('?dot',   REAL_TYPES)
!$:mfi_interface('sdsdot', REAL_TYPES)
$:mfi_implement('?dotu',  COMPLEX_TYPES, dot_product)
$:mfi_implement('?dotc',  COMPLEX_TYPES, dot_product)
!$:mfi_interface('?nrm2',  DEFAULT_TYPES)
!$:mfi_interface('?rot',   DEFAULT_TYPES)
!$:mfi_interface('?rotg',  DEFAULT_TYPES)
$:mfi_implement('?rotm',  REAL_TYPES, rotm)
!$:mfi_implement('?rotmg', REAL_TYPES, rotmg)
!$:f77_interface('?scal')
$:mfi_implement('?swap',  DEFAULT_TYPES, copy_swap)

! BLAS level 2
$:mfi_implement('?gbmv',  DEFAULT_TYPES, gbmv)
$:mfi_implement('?gemv',  DEFAULT_TYPES, gemv)
$:mfi_implement('?ger',   REAL_TYPES,    ger_gerc_geru)
$:mfi_implement('?gerc',  COMPLEX_TYPES, ger_gerc_geru)
$:mfi_implement('?geru',  COMPLEX_TYPES, ger_gerc_geru)
$:mfi_implement('?hbmv',  COMPLEX_TYPES, hbmv_sbmv)
$:mfi_implement('?hemv',  COMPLEX_TYPES, hemv_symv)
$:mfi_implement('?her',   COMPLEX_TYPES, her)
$:mfi_implement('?her2',  COMPLEX_TYPES, her_syr2)
$:mfi_implement('?hpmv',  COMPLEX_TYPES, hpmv_spmv)
$:mfi_implement('?hpr',   COMPLEX_TYPES, hpr)
$:mfi_implement('?hpr2',  COMPLEX_TYPES, hpr_spr2)
$:mfi_implement('?sbmv',  REAL_TYPES,    hbmv_sbmv)
$:mfi_implement('?spmv',  REAL_TYPES,    hpmv_spmv)
$:mfi_implement('?spr',   REAL_TYPES,    spr)
$:mfi_implement('?spr2',  REAL_TYPES,    hpr_spr2)
$:mfi_implement('?symv',  REAL_TYPES,    hemv_symv)
$:mfi_implement('?syr',   REAL_TYPES,    syr)
$:mfi_implement('?syr2',  REAL_TYPES,    her_syr2)
$:mfi_implement('?tbmv',  DEFAULT_TYPES, tbmv_tbsv)
$:mfi_implement('?tbsv',  DEFAULT_TYPES, tbmv_tbsv)
$:mfi_implement('?tpmv',  DEFAULT_TYPES, tpmv_tpsv)
$:mfi_implement('?tpsv',  DEFAULT_TYPES, tpmv_tpsv)
$:mfi_implement('?trmv',  DEFAULT_TYPES, trmv_trsv)
$:mfi_implement('?trsv',  DEFAULT_TYPES, trmv_trsv)

! BLAS level 3
$:mfi_implement('?gemm',  DEFAULT_TYPES, gemm)
$:mfi_implement('?hemm',  COMPLEX_TYPES, hemm_symm)
$:mfi_implement('?herk',  COMPLEX_TYPES, herk)
$:mfi_implement('?her2k', COMPLEX_TYPES, her2k)
$:mfi_implement('?symm',  REAL_TYPES,    hemm_symm)
$:mfi_implement('?syrk',  REAL_TYPES,    syrk)
$:mfi_implement('?syr2k', REAL_TYPES,    syr2k)
$:mfi_implement('?trmm',  DEFAULT_TYPES, trmm_trsm)
$:mfi_implement('?trsm',  DEFAULT_TYPES, trmm_trsm)

! Extensions
! BLAS level 1 - Utils / Extensions
$:mfi_implement('i?amax', DEFAULT_TYPES, iamin_iamax)
#:if defined('MFI_EXTENSIONS')
$:mfi_implement('i?amin', DEFAULT_TYPES, iamin_iamax)
#:endif
$:mfi_implement('?lamch', REAL_TYPES, lamch)

end module
