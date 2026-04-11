#:include "common.fpp"
#:include "test/lapack/macros/geqrf_gerqf.fypp"
#:include "test/lapack/macros/gesvd.fypp"
#:include "test/lapack/macros/getrf.fypp"
#:include "test/lapack/macros/getri.fypp"
#:include "test/lapack/macros/getrs.fypp"
#:include "test/lapack/macros/heevr.fypp"
#:include "test/lapack/macros/hetrf.fypp"
#:include "test/lapack/macros/orgqr.fypp"
#:include "test/lapack/macros/ormqr.fypp"
#:include "test/lapack/macros/pocon.fypp"
#:include "test/lapack/macros/potrf_potri.fypp"
#:include "test/lapack/macros/sytrf.fypp"
#:include "test/lapack/macros/trtrs.fypp"

#:set _COLLECT = [                                                          &
     ('?geqrf', DEFAULT_TYPES,    geqrf_gerqf),                             &
     ('?gerqf', DEFAULT_TYPES,    geqrf_gerqf),                             &
     ('?gesvd', DEFAULT_TYPES,    gesvd),                                   &
     ('?getrf', DEFAULT_TYPES,    getrf),                                   &
     ('?getri', DEFAULT_TYPES,    getri),                                   &
     ('?getrs', DEFAULT_TYPES,    getrs),                                   &
     ('?heevr', COMPLEX_TYPES,    heevr),                                   &
     ('?hetrf', COMPLEX_TYPES,    hetrf),                                   &
     ('?orgqr', REAL_TYPES,       orgqr),                                   &
     ('?orgrq', REAL_TYPES,       orgqr),                                   &
     ('?ormqr', REAL_TYPES,       ormqr),                                   &
     ('?pocon', DEFAULT_TYPES,    pocon),                                   &
     ('?potrf', DEFAULT_TYPES,    potrf_potri),                             &
     ('?sytrf', REAL_TYPES,       sytrf),                                   &
     ('?trtrs', DEFAULT_TYPES,    trtrs),                                   &
     ('?ungqr', COMPLEX_TYPES,    orgqr),                                   &
     ('?ungrq', COMPLEX_TYPES,    orgqr),                                   &
     ('?unmqr', COMPLEX_TYPES,    ormqr),                                   &
]
#:set _modname = lambda name: name.replace('?','')

#:for name, supported_types, code in _COLLECT
 program test_${_modname(name)}$
 use iso_fortran_env
 use mfi_lapack
 #:set f77_names = ', '.join([prefix(pfx, name) for pfx in supported_types])
 use f77_lapack, only: ${f77_names}$
 implicit none
 #:for pfx in supported_types
 #:set f77 = prefix(pfx, name)
 #:set mfi = 'mfi_' + prefix('', name)
 print '(A)', "testing ${mfi}$ (CPU) against ${f77}$"
 #:endfor
 contains
 $:test_implement(name, supported_types, code)

 #:include "test/assert.inc"

 end program

#:endfor

#:for name, supported_types, code in _COLLECT
 program test_${_modname(name)}$_gpu
 use iso_fortran_env
 use mfi_lapack
 #:set f77_names = ', '.join([prefix(pfx, name) for pfx in supported_types])
 use f77_lapack, only: ${f77_names}$
 implicit none
 #:for pfx in supported_types
 #:set f77 = prefix(pfx, name)
 #:set mfi = 'mfi_' + prefix('', name)
 print '(A)', "testing ${mfi}$ (GPU) against ${f77}$"
 #:endfor
 contains
 $:test_implement_gpu(name, supported_types, code)

 #:include "test/assert.inc"

 end program

#:endfor
