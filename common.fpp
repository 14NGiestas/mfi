#:set REAL_KINDS=['s','d']
#:set COMPLEX_KINDS=['c','z']
#:set DEFAULT_KINDS=REAL_KINDS+COMPLEX_KINDS
#:set PREFIX_TO_TYPE={  &
's': 'real(REAL32)',    &
'd': 'real(REAL64)',    &
'c': 'complex(REAL32)', &
'z': 'complex(REAL64)', &
}

#:set BLAS_ROUTINES = {  &
'?rotg':  REAL_KINDS,    &
'?rotmg': REAL_KINDS,    &
'?rot':   REAL_KINDS,    &
'?rotm':  REAL_KINDS,    &
'?swap':  DEFAULT_KINDS, &
'?scal':  DEFAULT_KINDS+['cs','zd'], &
'?copy':  DEFAULT_KINDS, &
'?axpy':  DEFAULT_KINDS, &
'?dot':   REAL_KINDS+['ds','sds'], &
'?dotu':  COMPLEX_KINDS, &
'?dotc':  COMPLEX_KINDS, &
'?nrm2':  REAL_KINDS+['sc','dz'], &
'?asum':  REAL_KINDS+['sc','dz'], &
'i?amax': DEFAULT_KINDS, &
'i?amin': DEFAULT_KINDS, &
'?gemv':  DEFAULT_KINDS, &
'?gbmv':  DEFAULT_KINDS, &
'?hemv':  COMPLEX_KINDS, &
'?hbmv':  COMPLEX_KINDS, &
'?hpmv':  COMPLEX_KINDS, &
'?symv':  REAL_KINDS,    &
'?sbmv':  REAL_KINDS,    &
'?spmv':  REAL_KINDS,    &
'?trmv':  DEFAULT_KINDS, &
'?tbmv':  DEFAULT_KINDS, &
'?tpmv':  DEFAULT_KINDS, &
'?trsv':  DEFAULT_KINDS, &
'?tbsv':  DEFAULT_KINDS, &
'?tpsv':  DEFAULT_KINDS, &
'?ger':   REAL_KINDS,    &
'?geru':  COMPLEX_KINDS, &
'?gerc':  COMPLEX_KINDS, &
'?her':   COMPLEX_KINDS, &
'?hpr':   COMPLEX_KINDS, &
'?her2':  COMPLEX_KINDS, &
'?hpr2':  COMPLEX_KINDS, &
'?syr':   REAL_KINDS,    &
'?spr':   REAL_KINDS,    &
'?syr2':  REAL_KINDS,    &
'?spr2':  REAL_KINDS,    &
'?gemm':  DEFAULT_KINDS, &
'?symm':  DEFAULT_KINDS, &
'?hemm':  COMPLEX_KINDS, &
'?syrk':  DEFAULT_KINDS, &
'?herk':  COMPLEX_KINDS, &
'?syr2k': DEFAULT_KINDS, &
'?her2k': COMPLEX_KINDS, &
'?trmm':  DEFAULT_KINDS, &
'?trsm':  DEFAULT_KINDS, &
'?axpyi': DEFAULT_KINDS, &
'?doti':  REAL_KINDS,    &
'?dotci': COMPLEX_KINDS, &
'?dotui': COMPLEX_KINDS, &
'?gthr':  DEFAULT_KINDS, &
'?gthrz': DEFAULT_KINDS, &
'?sctr':  DEFAULT_KINDS, &
'?roti':  REAL_KINDS     &
}
#:set ROUTINES = BLAS_ROUTINES

#:def f77_interface(name, types)
interface f77_${name.replace('?','')}$
    #:for T in types
    module procedure ${name.replace('?',T)}$
    #:endfor
end interface
#:enddef

#:def mfi_interface(name, types)
interface ${name.replace('?','')}$
    #:for T in types
    module procedure mfi_${name.replace('?',T)}$
    #:endfor
end interface
#:enddef

#:def mfi_implement(name, args, code)
#:for T in ROUTINES[name]
pure subroutine mfi_${name.replace('?',T)}$${args}$
${code.replace('{T}', PREFIX_TO_TYPE[T])}$
end subroutine

#:endfor
#:enddef

