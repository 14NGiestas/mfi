#:mute
#:include "common.fpp"
#:include "src/f77/blas/asum_nrm2.fypp"
#:include "src/f77/blas/scal.fypp"
#:include "src/f77/blas/axpy.fypp"
#:include "src/f77/blas/copy_swap.fypp"
#:include "src/f77/blas/dot_product.fypp"
#:include "src/f77/blas/sdsdot.fypp"
#:include "src/f77/blas/rot.fypp"
#:include "src/f77/blas/rotm.fypp"
#:include "src/f77/blas/rotmg.fypp"
#:include "src/f77/blas/gbmv.fypp"
#:include "src/f77/blas/gemv.fypp"
#:include "src/f77/blas/ger_gerc_geru.fypp"
#:include "src/f77/blas/hbmv_sbmv.fypp"
#:include "src/f77/blas/hemv_symv.fypp"
#:include "src/f77/blas/her.fypp"
#:include "src/f77/blas/syr.fypp"
#:include "src/f77/blas/her_syr2.fypp"
#:include "src/f77/blas/hpmv_spmv.fypp"
#:include "src/f77/blas/hpr.fypp"
#:include "src/f77/blas/spr.fypp"
#:include "src/f77/blas/hpr_spr2.fypp"
#:include "src/f77/blas/tbmv_tbsv.fypp"
#:include "src/f77/blas/tpmv_tpsv.fypp"
#:include "src/f77/blas/trmv_trsv.fypp"
#:include "src/f77/blas/gemm.fypp"
#:include "src/f77/blas/hemm_symm.fypp"
#:include "src/f77/blas/herk.fypp"
#:include "src/f77/blas/syrk.fypp"
#:include "src/f77/blas/her2k.fypp"
#:include "src/f77/blas/syr2k.fypp"
#:include "src/f77/blas/trmm_trsm.fypp"
! BLAS Level 1 - Extensions
#:include "src/f77/blas/iamax_iamin.fypp"
#:include "src/f77/blas/iamin_stub.fypp"
#:endmute
!> Improved and original F77 interfaces for blas
module f77_blas
use iso_fortran_env
implicit none

!FIXME rot, dot, rotg, nrm2: problem with functions that have TYPE /= TYPE_result

! BLAS level 1

$:f77_interface('?axpy',  DEFAULT_TYPES, axpy)
$:f77_interface('?copy',  DEFAULT_TYPES, copy_swap)
$:f77_interface('?dot',   REAL_TYPES,    dot_product)
$:f77_interface('?dotu',  COMPLEX_TYPES, dot_product)
$:f77_interface('?dotc',  COMPLEX_TYPES, dot_product)
!$:f77_interface('?rotg', DEFAULT_TYPES, rotg, result=REAL_TYPES)
$:f77_interface('?rotm',  REAL_TYPES,    rotm)
$:f77_interface('?rotmg', REAL_TYPES,    rotmg)
$:f77_interface('?swap',  DEFAULT_TYPES, copy_swap)

#! Problematic functions
#! asum has special names indicating the returns are real types
$:f77_interface('?asum',  DEFAULT_TYPES, asum_nrm2, f=MIX_REAL_COMPLEX)
#! nrm2 has the same interface as asum
$:f77_interface('?nrm2',  DEFAULT_TYPES, asum_nrm2, f=MIX_REAL_COMPLEX)
#! scal has mixed types scalars so it can multiply a real constant by a complex vector
$:f77_interface('?scal',  DEFAULT_TYPES, scal,       improved_f77=False)
$:f77_interface('?scal',  COMPLEX_TYPES, scal_mixed, improved_f77=False, f=MIX_COMPLEX_REAL)
$:f77_interface_improved('?scal', DEFAULT_TYPES + COMPLEX_REAL_TYPES)

#! ?rot has mixed types scalars so it can multiply a real constant by a complex vector
$:f77_interface('?rot',  DEFAULT_TYPES, rot,       improved_f77=False)
$:f77_interface('?rot',  COMPLEX_TYPES, rot_mixed, improved_f77=False, f=MIX_COMPLEX_REAL)
$:f77_interface_improved('?rot', DEFAULT_TYPES + COMPLEX_REAL_TYPES)

! BLAS level 2
$:f77_interface('?gbmv',  DEFAULT_TYPES, gbmv)
$:f77_interface('?gemv',  DEFAULT_TYPES, gemv)
$:f77_interface('?ger',   REAL_TYPES,    ger_gerc_geru)
$:f77_interface('?gerc',  COMPLEX_TYPES, ger_gerc_geru)
$:f77_interface('?geru',  COMPLEX_TYPES, ger_gerc_geru)
$:f77_interface('?hbmv',  COMPLEX_TYPES, hbmv_sbmv)
$:f77_interface('?hemv',  COMPLEX_TYPES, hemv_symv)
$:f77_interface('?her',   COMPLEX_TYPES, her)
$:f77_interface('?her2',  COMPLEX_TYPES, her_syr2)
$:f77_interface('?hpmv',  COMPLEX_TYPES, hpmv_spmv)
$:f77_interface('?hpr',   COMPLEX_TYPES, hpr)
$:f77_interface('?hpr2',  COMPLEX_TYPES, hpr_spr2)
$:f77_interface('?sbmv',  REAL_TYPES,    hbmv_sbmv)
$:f77_interface('?spmv',  REAL_TYPES,    hpmv_spmv)
$:f77_interface('?spr',   REAL_TYPES,    spr)
$:f77_interface('?spr2',  REAL_TYPES,    hpr_spr2)
$:f77_interface('?symv',  REAL_TYPES,    hemv_symv)
$:f77_interface('?syr',   REAL_TYPES,    syr)
$:f77_interface('?syr2',  REAL_TYPES,    her_syr2)
$:f77_interface('?tbmv',  DEFAULT_TYPES, tbmv_tbsv)
$:f77_interface('?tbsv',  DEFAULT_TYPES, tbmv_tbsv)
$:f77_interface('?tpmv',  DEFAULT_TYPES, tpmv_tpsv)
$:f77_interface('?tpsv',  DEFAULT_TYPES, tpmv_tpsv)
$:f77_interface('?trmv',  DEFAULT_TYPES, trmv_trsv)
$:f77_interface('?trsv',  DEFAULT_TYPES, trmv_trsv)

! BLAS level 3
$:f77_interface('?gemm',  DEFAULT_TYPES, gemm)
$:f77_interface('?hemm',  COMPLEX_TYPES, hemm_symm)
$:f77_interface('?herk',  COMPLEX_TYPES, herk)
$:f77_interface('?her2k', COMPLEX_TYPES, her2k)
$:f77_interface('?symm',  REAL_TYPES,    hemm_symm)
$:f77_interface('?syrk',  REAL_TYPES,    syrk)
$:f77_interface('?syr2k', REAL_TYPES,    syr2k)
$:f77_interface('?trmm',  DEFAULT_TYPES, trmm_trsm)
$:f77_interface('?trsm',  DEFAULT_TYPES, trmm_trsm)


#:include "src/f77/blas/specific_interfaces.fypp"

! Extensions
! BLAS Level 1 - Utils / Extensions
#:if defined('MFI_EXTENSIONS')
  #:if defined('MFI_LINK_EXTERNAL')
! Link with a external source
$:f77_interface('i?amax', DEFAULT_TYPES, iamax_iamin)
$:f77_interface('i?amin', DEFAULT_TYPES, iamax_iamin)
  #:else
! Implement the blas extensions in
$:f77_interface_improved('i?amax', DEFAULT_TYPES)
$:f77_interface_improved('i?amin', DEFAULT_TYPES)
contains
$:f77_implement('i?amax', DEFAULT_TYPES, iamin_stub)
$:f77_implement('i?amin', DEFAULT_TYPES, iamin_stub)
  #:endif
#:endif

end module

