#:mute
#:include "common.fpp"
#:include "src/mfi/blas/lamch.fypp"
#:include "src/mfi/blas/asum_nrm2.fypp"
#:include "src/mfi/blas/axpy.fypp"
#:include "src/mfi/blas/scal.fypp"
#:include "src/mfi/blas/copy_swap.fypp"
#:include "src/mfi/blas/dot_product.fypp"
#:include "src/mfi/blas/rot.fypp"
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
#:set COLLECT = [                                              &
    ('?copy', DEFAULT_TYPES,                    copy_swap),    &
    ('?swap', DEFAULT_TYPES,                    copy_swap),    &
    ('?axpy', DEFAULT_TYPES,                    axpy),         &
    ('?dot',  REAL_TYPES,                       dot_product),  &
    ('?dotc', COMPLEX_TYPES,                    dot_product),  &
    ('?dotu', COMPLEX_TYPES,                    dot_product),  &
    ('?asum', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),    &
    ('?nrm2', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),    &
    ('?rot',  DEFAULT_TYPES + MIX_COMPLEX_REAL, rot),          &
    ('?rotm', REAL_TYPES,                       rotm),         &
    ('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL, scal),         &
    ('?gbmv', DEFAULT_TYPES,                    gbmv),         &
    ('?gemv', DEFAULT_TYPES,                    gemv),         &
    ('?ger',  REAL_TYPES,                       ger_gerc_geru),&
    ('?gerc', COMPLEX_TYPES,                    ger_gerc_geru),&
    ('?geru', COMPLEX_TYPES,                    ger_gerc_geru),&
    ('?hbmv', COMPLEX_TYPES,                    hbmv_sbmv),    &
    ('?hemv', COMPLEX_TYPES,                    hemv_symv),    &
    ('?her',  COMPLEX_TYPES,                    her),          &
    ('?her2', COMPLEX_TYPES,                    her_syr2),     &
    ('?hpmv', COMPLEX_TYPES,                    hpmv_spmv),    &
    ('?hpr',  COMPLEX_TYPES,                    hpr),          &
    ('?hpr2', COMPLEX_TYPES,                    hpr_spr2),     &
    ('?sbmv', REAL_TYPES,                       hbmv_sbmv),    &
    ('?spmv', REAL_TYPES,                       hpmv_spmv),    &
    ('?spr',  REAL_TYPES,                       spr),          &
    ('?spr2', REAL_TYPES,                       hpr_spr2),     &
    ('?symv', REAL_TYPES,                       hemv_symv),    &
    ('?syr',  REAL_TYPES,                       syr),          &
    ('?syr2', REAL_TYPES,                       her_syr2),     &
    ('?tbmv', DEFAULT_TYPES,                    tbmv_tbsv),    &
    ('?tbsv', DEFAULT_TYPES,                    tbmv_tbsv),    &
    ('?tpmv', DEFAULT_TYPES,                    tpmv_tpsv),    &
    ('?tpsv', DEFAULT_TYPES,                    tpmv_tpsv),    &
    ('?trmv', DEFAULT_TYPES,                    trmv_trsv),    &
    ('?trsv', DEFAULT_TYPES,                    trmv_trsv),    &
    ('?gemm', DEFAULT_TYPES,                    gemm),         &
    ('?hemm', COMPLEX_TYPES,                    hemm_symm),    &
    ('?herk', COMPLEX_TYPES,                    herk),         &
    ('?her2k',COMPLEX_TYPES,                    her2k),        &
    ('?symm', REAL_TYPES,                       hemm_symm),    &
    ('?syrk', REAL_TYPES,                       syrk),         &
    ('?syr2k',REAL_TYPES,                       syr2k),        &
    ('?trmm', DEFAULT_TYPES,                    trmm_trsm),    &
    ('?trsm', DEFAULT_TYPES,                    trmm_trsm),    &
    ('?lamch',REAL_TYPES,                       lamch),        &
]
#:endmute
!> Modern fortran interfaces for BLAS
module mfi_blas
use iso_fortran_env
use f77_blas
use f77_blas, only: mfi_rotg  => f77_rotg
use f77_blas, only: mfi_rotmg => f77_rotmg
implicit none

#:for name, supported_types, code in COLLECT
$:mfi_interface(name, supported_types)
#:endfor

! Extensions
! BLAS level 1 - Utils / Extensions
#:if defined('MFI_EXTENSIONS')
$:mfi_interface('i?amax', DEFAULT_TYPES)
$:mfi_interface('i?amin', DEFAULT_TYPES)
#:endif


#:if defined('USE_CUBLAS')
$:cublas_interfaces()
#:endif

contains


#:for name, supported_types, code in COLLECT
$:mfi_implement(name, supported_types, code)
#:endfor

! Extensions
! BLAS level 1 - Utils / Extensions
#:if defined('MFI_EXTENSIONS')
$:mfi_implement('i?amax', DEFAULT_TYPES, iamin_iamax)
$:mfi_implement('i?amin', DEFAULT_TYPES, iamin_iamax)
#:endif


end module
