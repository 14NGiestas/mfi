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
'?scal':  DEFAULT_KINDS, &
'?copy':  DEFAULT_KINDS, &
'?axpy':  DEFAULT_KINDS, &
'?dot':   REAL_KINDS,    &
'?dotu':  COMPLEX_KINDS, &
'?dotc':  COMPLEX_KINDS, &
'?nrm2':  REAL_KINDS,    &
'?asum':  REAL_KINDS,    &
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

#:def optargs(dtype, *args)
#:for variable in args
    ${dtype}$, intent(in), optional :: ${variable}$
    ${dtype}$ :: local_${variable}$
#:endfor
#:enddef

#:def localvars(dtype, *args)
#:for variable in args
    ${dtype}$ :: ${variable}$
#:endfor
#:enddef

#:def defaults(**kwargs)
#:for variable, default in kwargs.items()
    if (present(${variable}$)) then
        local_${variable}$ = ${variable}$
    else
        local_${variable}$ = ${default}$
    end if
#:endfor
#:enddef

#:def args(dtype, intent, *args)
#:for variable in args
    ${dtype}$, intent(${intent}$), :: ${variable}$
#:endfor
#:enddef


#:def mfi_implement(name, code)
#:for P in BLAS_ROUTINES[name]
#:set T=PREFIX_TO_TYPE[P]
#:set N=name.replace('?',P)
$:code.replace('SNAME', N) &
      .replace('DTYPE', T)
#:endfor
#:enddef
