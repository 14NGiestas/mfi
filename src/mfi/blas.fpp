#:mute
#:include "common.fpp"
#:include "cublas.fpp"
#:include "extensions.fpp"
#! MFI BLAS macros
#:include "src/mfi/blas/axpy.fypp"
#:include "src/mfi/blas/asum_nrm2.fypp"
#:include "src/mfi/blas/copy_swap.fypp"
#:include "src/mfi/blas/dot_product.fypp"
#:include "src/mfi/blas/gbmv.fypp"
#:include "src/mfi/blas/gemv.fypp"
#:include "src/mfi/blas/gemm.fypp"
#:include "src/mfi/blas/ger_gerc_geru.fypp"
#:include "src/mfi/blas/hbmv_sbmv.fypp"
#:include "src/mfi/blas/hemv_symv.fypp"
#:include "src/mfi/blas/hemm_symm.fypp"
#:include "src/mfi/blas/her.fypp"
#:include "src/mfi/blas/herk.fypp"
#:include "src/mfi/blas/her2k.fypp"
#:include "src/mfi/blas/her_syr2.fypp"
#:include "src/mfi/blas/hpmv_spmv.fypp"
#:include "src/mfi/blas/hpr.fypp"
#:include "src/mfi/blas/hpr_spr2.fypp"
#:include "src/mfi/blas/iamin_iamax.fypp"
#:include "src/mfi/blas/lamch.fypp"
#:include "src/mfi/blas/rot.fypp"
#:include "src/mfi/blas/rotm.fypp"
#:include "src/mfi/blas/scal.fypp"
#:include "src/mfi/blas/spr.fypp"
#:include "src/mfi/blas/syr.fypp"
#:include "src/mfi/blas/syrk.fypp"
#:include "src/mfi/blas/syr2k.fypp"
#:include "src/mfi/blas/tbmv_tbsv.fypp"
#:include "src/mfi/blas/tpmv_tpsv.fypp"
#:include "src/mfi/blas/trmm_trsm.fypp"
#:include "src/mfi/blas/trmv_trsv.fypp"
#:set _COLLECT = [                                               &
    ('?copy', DEFAULT_TYPES,                    copy_swap),        &
    ('?swap', DEFAULT_TYPES,                    copy_swap),        &
    ('?axpy', DEFAULT_TYPES,                    axpy),             &
    ('?dot',  REAL_TYPES,                       dot_product),      &
    ('?dotc', COMPLEX_TYPES,                    dot_product),      &
    ('?dotu', COMPLEX_TYPES,                    dot_product),      &
    ('?asum', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),        &
    ('?nrm2', REAL_TYPES    + MIX_REAL_COMPLEX, asum_nrm2),        &
    ('?rot',  DEFAULT_TYPES + MIX_COMPLEX_REAL, rot),              &
    ('?rotm', REAL_TYPES,                       rotm),             &
    ('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL, scal),             &
    ('?gbmv', DEFAULT_TYPES,                    gbmv),             &
    ('?gemv', DEFAULT_TYPES,                    gemv),             &
    ('?ger',  REAL_TYPES,                       ger_gerc_geru),    &
    ('?gerc', COMPLEX_TYPES,                    ger_gerc_geru),    &
    ('?geru', COMPLEX_TYPES,                    ger_gerc_geru),    &
    ('?hbmv', COMPLEX_TYPES,                    hbmv_sbmv),        &
    ('?hemv', COMPLEX_TYPES,                    hemv_symv),        &
    ('?her',  COMPLEX_TYPES,                    her),              &
    ('?her2', COMPLEX_TYPES,                    her_syr2),         &
    ('?hpmv', COMPLEX_TYPES,                    hpmv_spmv),        &
    ('?hpr',  COMPLEX_TYPES,                    hpr),              &
    ('?hpr2', COMPLEX_TYPES,                    hpr_spr2),         &
    ('?sbmv', REAL_TYPES,                       hbmv_sbmv),        &
    ('?spmv', REAL_TYPES,                       hpmv_spmv),        &
    ('?spr',  REAL_TYPES,                       spr),              &
    ('?spr2', REAL_TYPES,                       hpr_spr2),         &
    ('?symv', REAL_TYPES,                       hemv_symv),        &
    ('?syr',  REAL_TYPES,                       syr),              &
    ('?syr2', REAL_TYPES,                       her_syr2),         &
    ('?tbmv', DEFAULT_TYPES,                    tbmv_tbsv),        &
    ('?tbsv', DEFAULT_TYPES,                    tbmv_tbsv),        &
    ('?tpmv', DEFAULT_TYPES,                    tpmv_tpsv),        &
    ('?tpsv', DEFAULT_TYPES,                    tpmv_tpsv),        &
    ('?trmv', DEFAULT_TYPES,                    trmv_trsv),        &
    ('?trsv', DEFAULT_TYPES,                    trmv_trsv),        &
    ('?gemm', DEFAULT_TYPES,                    gemm),             &
    ('?hemm', COMPLEX_TYPES,                    hemm_symm),        &
    ('?herk', COMPLEX_TYPES,                    herk),             &
    ('?her2k',COMPLEX_TYPES,                    her2k),            &
    ('?symm', REAL_TYPES,                       hemm_symm),        &
    ('?syrk', REAL_TYPES,                       syrk),             &
    ('?syr2k',REAL_TYPES,                       syr2k),            &
    ('?trmm', DEFAULT_TYPES,                    trmm_trsm),        &
    ('?trsm', DEFAULT_TYPES,                    trmm_trsm),        &
    ('?lamch',REAL_TYPES,                       lamch),            &
]
#:set _modname = lambda name: name.replace('?','')
#:endmute

#! Per-module submodules
#:for name, supported_types, code in _COLLECT
module mfi_blas_${_modname(name)}$
    use iso_fortran_env
    use f77_blas
#if defined(MFI_CUBLAS)
    use iso_c_binding
    use mfi_blas_cublas
#endif
#if defined(MFI_EXTENSIONS)
    use mfi_blas_extensions
#endif
    implicit none

    $:mfi_interface(name, supported_types)

contains

$:mfi_implement(name, supported_types, code)
end module

#:endfor

!> cuBLAS interfaces and constants
module mfi_blas_cublas
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    use iso_c_binding
    implicit none

$:cublas_interfaces()
    type(c_ptr), save :: mfi_cublas_handle = c_null_ptr
#endif
end module

!> Extensions module
module mfi_blas_extensions
#if defined(MFI_EXTENSIONS)
    use iso_fortran_env
    use f77_blas
#if defined(MFI_CUBLAS)
    use iso_c_binding
    use mfi_blas_cublas
#endif
    implicit none

$:mfi_extensions_interfaces()

contains

$:mfi_extensions_implement()
#endif
end module

#! Umbrella module — re-exports all submodules
module mfi_blas
    use iso_fortran_env
    use iso_c_binding
    use f77_blas
    use f77_blas, only: mfi_rotg  => f77_rotg
    use f77_blas, only: mfi_rotmg => f77_rotmg
#:for name, supported_types, code in _COLLECT
    use mfi_blas_${_modname(name)}$
#:endfor
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS) 
    use mfi_blas_cublas
    use mfi_blas_extensions
#elif defined(MFI_EXTENSIONS)
    use mfi_blas_extensions
#endif
    implicit none
end module
