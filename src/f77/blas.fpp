#:mute
#:include "common.fpp"
#! F77 BLAS macros
#:include "src/f77/blas/asum_nrm2.fypp"
#:include "src/f77/blas/axpy.fypp"
#:include "src/f77/blas/copy_swap.fypp"
#:include "src/f77/blas/dot_product.fypp"
#:include "src/f77/blas/gbmv.fypp"
#:include "src/f77/blas/gemv.fypp"
#:include "src/f77/blas/gemm.fypp"
#:include "src/f77/blas/ger_gerc_geru.fypp"
#:include "src/f77/blas/hbmv_sbmv.fypp"
#:include "src/f77/blas/hemv_symv.fypp"
#:include "src/f77/blas/hemm_symm.fypp"
#:include "src/f77/blas/her.fypp"
#:include "src/f77/blas/herk.fypp"
#:include "src/f77/blas/her2k.fypp"
#:include "src/f77/blas/her_syr2.fypp"
#:include "src/f77/blas/hpmv_spmv.fypp"
#:include "src/f77/blas/hpr.fypp"
#:include "src/f77/blas/hpr_spr2.fypp"
#:include "src/f77/blas/iamax_iamin.fypp"
#:include "src/f77/blas/iamin_stub.fypp"
#:include "src/f77/blas/rot.fypp"
#:include "src/f77/blas/rotg.fypp"
#:include "src/f77/blas/rotm.fypp"
#:include "src/f77/blas/rotmg.fypp"
#:include "src/f77/blas/scal.fypp"
#:include "src/f77/blas/spr.fypp"
#:include "src/f77/blas/syr.fypp"
#:include "src/f77/blas/syrk.fypp"
#:include "src/f77/blas/syr2k.fypp"
#:include "src/f77/blas/specific_interfaces.fypp"
#:include "src/f77/blas/tbmv_tbsv.fypp"
#:include "src/f77/blas/tpmv_tpsv.fypp"
#:include "src/f77/blas/trmm_trsm.fypp"
#:include "src/f77/blas/trmv_trsv.fypp"
#:set _COLLECT = [                                              &
    ('?copy', DEFAULT_TYPES,                    copy_swap),       &
    ('?swap', DEFAULT_TYPES,                    copy_swap),       &
    ('?axpy', DEFAULT_TYPES,                    axpy),            &
    ('?dot',  REAL_TYPES,                       dot_product),     &
    ('?dotc', COMPLEX_TYPES,                    dot_product),     &
    ('?dotu', COMPLEX_TYPES,                    dot_product),     &
    ('?asum', REAL_TYPES + MIX_REAL_COMPLEX,    asum_nrm2),       &
    ('?nrm2', REAL_TYPES + MIX_REAL_COMPLEX,    asum_nrm2),       &
    ('?rot',  DEFAULT_TYPES + MIX_COMPLEX_REAL, rot),             &
    ('?rotg', DEFAULT_TYPES,                    rotg),            &
    ('?rotm', REAL_TYPES,                       rotm),            &
    ('?rotmg',REAL_TYPES,                       rotmg),           &
    ('?scal', DEFAULT_TYPES + MIX_COMPLEX_REAL, scal),            &
    ('?gbmv', DEFAULT_TYPES,                    gbmv),            &
    ('?gemv', DEFAULT_TYPES,                    gemv),            &
    ('?ger',  REAL_TYPES,                       ger_gerc_geru),   &
    ('?gerc', COMPLEX_TYPES,                    ger_gerc_geru),   &
    ('?geru', COMPLEX_TYPES,                    ger_gerc_geru),   &
    ('?hbmv', COMPLEX_TYPES,                    hbmv_sbmv),       &
    ('?hemv', COMPLEX_TYPES,                    hemv_symv),       &
    ('?her',  COMPLEX_TYPES,                    her),             &
    ('?her2', COMPLEX_TYPES,                    her_syr2),        &
    ('?hpmv', COMPLEX_TYPES,                    hpmv_spmv),       &
    ('?hpr',  COMPLEX_TYPES,                    hpr),             &
    ('?hpr2', COMPLEX_TYPES,                    hpr_spr2),        &
    ('?sbmv', REAL_TYPES,                       hbmv_sbmv),       &
    ('?spmv', REAL_TYPES,                       hpmv_spmv),       &
    ('?spr',  REAL_TYPES,                       spr),             &
    ('?spr2', REAL_TYPES,                       hpr_spr2),        &
    ('?symv', REAL_TYPES,                       hemv_symv),       &
    ('?syr',  REAL_TYPES,                       syr),             &
    ('?syr2', REAL_TYPES,                       her_syr2),        &
    ('?tbmv', DEFAULT_TYPES,                    tbmv_tbsv),       &
    ('?tbsv', DEFAULT_TYPES,                    tbmv_tbsv),       &
    ('?tpmv', DEFAULT_TYPES,                    tpmv_tpsv),       &
    ('?tpsv', DEFAULT_TYPES,                    tpmv_tpsv),       &
    ('?trmv', DEFAULT_TYPES,                    trmv_trsv),       &
    ('?trsv', DEFAULT_TYPES,                    trmv_trsv),       &
    ('?gemm', DEFAULT_TYPES,                    gemm),            &
    ('?hemm', COMPLEX_TYPES,                    hemm_symm),       &
    ('?herk', COMPLEX_TYPES,                    herk),            &
    ('?her2k',COMPLEX_TYPES,                    her2k),           &
    ('?symm', REAL_TYPES,                       hemm_symm),       &
    ('?syrk', REAL_TYPES,                       syrk),            &
    ('?syr2k',REAL_TYPES,                       syr2k),           &
    ('?trmm', DEFAULT_TYPES,                    trmm_trsm),       &
    ('?trsm', DEFAULT_TYPES,                    trmm_trsm),       &
]
#:set _modname = lambda name: name.replace('?','')
#:endmute

#! Per-module submodules
#:for name, supported_types, code in _COLLECT
module f77_blas_${_modname(name)}$
    use iso_fortran_env
    use iso_c_binding
    implicit none

$:f77_original(name, supported_types, code)
end module

#:endfor

! Extensions
#:if defined('MFI_EXTENSIONS')
  #:if defined('MFI_LINK_EXTERNAL')
module f77_blas_iamax
    use iso_fortran_env
    use iso_c_binding
    implicit none
$:f77_original('i?amax', DEFAULT_TYPES, iamax_iamin)
end module
module f77_blas_iamin
    use iso_fortran_env
    use iso_c_binding
    implicit none
$:f77_original('i?amin', DEFAULT_TYPES, iamax_iamin)
end module
  #:else
module f77_blas_iamax
    use iso_fortran_env
    use iso_c_binding
    implicit none
$:f77_improved('i?amax', DEFAULT_TYPES)
contains
$:f77_implement('i?amax', DEFAULT_TYPES, iamin_stub)
end module
module f77_blas_iamin
    use iso_fortran_env
    use iso_c_binding
    implicit none
$:f77_improved('i?amin', DEFAULT_TYPES)
contains
$:f77_implement('i?amin', DEFAULT_TYPES, iamin_stub)
end module
  #:endif
#:endif

!> cuBLAS v2 interfaces
#:if defined('MFI_EXTENSIONS') and defined('MFI_USE_CUBLAS')
module f77_blas_cublas
use iso_c_binding
implicit none
interface
#:for pfx in DEFAULT_TYPES
#:set pfxs = split(pfx)
#:set NAME = prefix(pfx, '?gemm')
$:cublas_gemm(NAME, pfxs)
#:set NAME = prefix(pfx, '?gemv')
$:cublas_gemv(NAME, pfxs)
#:set NAME = prefix(pfx, '?trmm')
$:cublas_trmm_trsm(NAME, pfxs)
#:set NAME = prefix(pfx, '?trsm')
$:cublas_trmm_trsm(NAME, pfxs)
#:endfor
end interface
end module
#:endif

#! Umbrella module — re-exports all submodules
module f77_blas
    use iso_fortran_env
    use iso_c_binding
#:for name, supported_types, code in _COLLECT
    use f77_blas_${_modname(name)}$
#:endfor
#:if defined('MFI_EXTENSIONS')
    use f77_blas_iamax
    use f77_blas_iamin
#:endif
#:if defined('MFI_EXTENSIONS') and defined('MFI_USE_CUBLAS')
    use f77_blas_cublas
#:endif
    implicit none

#:include "src/f77/blas/specific_interfaces.fypp"
end module
