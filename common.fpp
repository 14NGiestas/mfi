#:mute

#:set REAL    = 'real(wp)'
#:set COMPLEX = 'complex(wp)'
#:set PREFIX = { &
    's': { 'type': 'real(wp)',    'wp': 'REAL32'}, &
    'd': { 'type': 'real(wp)',    'wp': 'REAL64'}, &
    'c': { 'type': 'complex(wp)', 'wp': 'REAL32'}, &
    'z': { 'type': 'complex(wp)', 'wp': 'REAL64'}, &
}

#:set ERROR = lambda pfx: { 'type': f'error: {pfx}', 'wp' : f'error: {pfx}' }

#:set mix    = lambda l, r: list(lp + rp for lp, rp in zip(l,r))
#:set split  = lambda pfx: list(pfx) if len(pfx) > 1 else pfx
#:set get_types = lambda pfxs: (pfxs[0], pfxs[0] if len(pfxs) == 1 else pfxs[1])
#:set get    = lambda pfx,what: PREFIX.get(pfx).get(what)
#:set prefix = lambda pfx, name: name.replace('?',pfx)
#:set kind   = lambda pfx: get(pfx,'wp')
#:set type   = lambda pfx: get(pfx,'type').replace('wp',kind(pfx))
#:set real   = lambda pfx: REAL.replace('wp',kind(pfx))
#:set complex= lambda pfx: COMPLEX.replace('wp',kind(pfx))

#:set SINGLE_TYPES  = ['s','c']
#:set DOUBLE_TYPES  = ['d','z']
#:set REAL_TYPES    = ['s','d']
#:set COMPLEX_TYPES = ['c','z']
#:set DEFAULT_TYPES = REAL_TYPES + COMPLEX_TYPES

#:set MIX_REAL_COMPLEX  = mix(REAL_TYPES,COMPLEX_TYPES)
#:set MIX_COMPLEX_REAL  = mix(COMPLEX_TYPES,REAL_TYPES)
#:set MIX_SINGLE_DOUBLE = mix(SINGLE_TYPES,DOUBLE_TYPES)
#:set MIX_DOUBLE_SINGLE = mix(DOUBLE_TYPES,SINGLE_TYPES)

${type('s')}$ :: variable
${get('s','wp')}$
${prefix('s','?gemm')}$

${mix(REAL_TYPES,COMPLEX_TYPES)}$
${mix(COMPLEX_TYPES,REAL_TYPES)}$
${mix(SINGLE_TYPES,DOUBLE_TYPES)}$
${mix(DOUBLE_TYPES,SINGLE_TYPES)}$

${list(map(split, mix(REAL_TYPES, COMPLEX_TYPES)))}$
${list(map(split, mix(COMPLEX_TYPES, REAL_TYPES)))}$
${list(map(split, mix(SINGLE_TYPES,DOUBLE_TYPES)))}$
${list(map(split, mix(DOUBLE_TYPES,SINGLE_TYPES)))}$

#:def timeit(message, code)
block
real :: t1, t2
call cpu_time(t1)
$:code
call cpu_time(t2)
print '(A," (",G0,"s)")', ${message}$, t2-t1
end block
#:enddef

#:def random_number(type, name, shape='')
#:if type.startswith('complex')
    $:random_complex(type, name,shape)
#:else
    call random_number(${name}$)
#:endif
#:enddef

#:def random_complex(type, name, shape='')
#:set REAL = type.replace('complex','real')
block
    ${REAL}$ :: re${shape}$
    ${REAL}$ :: im${shape}$
    call random_number(im)
    call random_number(re)
    ${name}$ = cmplx(re,im)
end block
#:enddef

#! Handles parameters (usage: working precision)
#:def parameter(dtype, **kwargs)
#:for variable, value in kwargs.items()
    ${dtype}$, parameter :: ${variable}$ = ${value}$
#:endfor
#:enddef

#! Handles importing and setting precision constants in interfaces
#:def imports(pfxs)
#:set wps = set(list(map(kind, pfxs)))
#:if len(wps) > 1
    import :: ${', '.join(wps)}$
#:else
    import :: ${''.join(wps)}$
#:endif
#:enddef

#! Handles the input/output arguments
#:def args(dtype, intent, *args)
#:for variable in args
    ${dtype}$, intent(${intent}$) :: ${variable}$
#:endfor
#:enddef

#! Defines a optional variable, creating local corresponding variable by default
#:def optional(dtype, intent, *args)
#:for variable in args
    ${dtype}$, intent(${intent}$), optional :: ${variable}$
    ${dtype}$ :: local_${variable}$
#:endfor
#:enddef

#! Handles default values of a optional variable
#:def defaults(**kwargs)
#:for variable, default in kwargs.items()
    if (present(${variable}$)) then
        local_${variable}$ = ${variable}$
    else
        local_${variable}$ = ${default}$
    end if
#:endfor
#:enddef

#! Handles a value of "variable" depending on "condition"
#:def optval(condition, variable, true_value, false_value)
    if (${condition}$) then
        ${variable}$ = ${true_value}$
    else
        ${variable}$ = ${false_value}$
    end if
#:enddef

#:def interface(functions, procedure='procedure', name='')
interface ${name}$
    #:for function_name in functions
    ${procedure}$ :: ${function_name}$
    #:endfor
end interface
#:enddef

#! Interfaces for the original f77 routines
#! code must implement a routine interface
#:def f77_original(generic_name, prefixes, code)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set f77 = 'f77_' + prefix('',generic_name)
!> ${generic_name}$ supports ${', '.join(prefixes)}$.
!> See also: [[${mfi}$]], [[${f77}$]].
interface
#:for pfx in prefixes
#:set name = prefix(pfx,generic_name)
#:set pfxs = list(map(split,pfx))
$:code(name,pfxs)
#:endfor
end interface
#:enddef

#! Define a common interface with the original f77 interfaces
#! So you can call the original function without the prefix
#:def f77_improved(generic_name, prefixes)
#:set functions = map(lambda pfx: prefix(pfx,generic_name), prefixes)
$:interface(functions, name=f"f77_{prefix('',generic_name)}")
#:enddef

#! In case of missing functions / extensions you can pass a code
#! in which case must provide the routine implementation
#! Must be called inside a contains block
#:def f77_implement(generic_name, prefixes, code)
#:for pfx in prefixes
#:set name = prefix(pfx,generic_name)
#:set pfxs = list(map(split,pfx))
$:code(name,pfxs)
#:endfor
#:enddef

#:def mfi_interface(generic_name, prefixes)
#:set functions = map(lambda pfx: 'mfi_' + prefix(pfx,generic_name), prefixes)
$:interface(functions, &
            procedure='module procedure', &
            name=f"mfi_{prefix('',generic_name)}")
#:enddef

#! Implements the modern interface in code
#! for each supported prefix combination
#! Must be called inside a contains block
#:def mfi_implement(generic_name, prefixes, code)
#:for pfx in prefixes
#:set mfi_name  = 'mfi_' + prefix(pfx,generic_name)
#:set f77_name  =          prefix(pfx,generic_name)
#:set pfxs      = list(map(split,pfx))
$:code(mfi_name,f77_name,pfxs)
#:endfor
#:enddef


#! Implements the test for all interfaces
#! and each supported prefix combination
#! Must be called inside a contains block
#:def test_implement(generic_name, prefixes, code)
#:for pfx in prefixes
#:set f77  =          prefix(pfx,generic_name)
#:set f90 = 'f77_' + prefix('',generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set pfxs    = list(map(split,pfx))
$:code(f77,f90,mfi,pfxs)
#:endfor
#:enddef

#:def test_run(generic_name, prefixes)
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
@:timeit("testing ${mfi}$ against ${f77}$", { call test_${f77}$ })
#:endfor
#:enddef

#:def rot_f77(name,pfxs)
#:set A = pfxs[0]
#:set B = A if len(pfxs) == 1 else pfxs[1]
!> ${name.upper()}$ applies a plane rotation.
pure subroutine ${name}$(n, x, incx, y, incy, c, s)
    $:imports(pfxs)
@:args(${type(A)}$, in, x(*), y(*))
@:args(${real(A)}$, in, c)
@:args(${type(B)}$, in, s)
    integer,     intent(in) :: n, incx, incy
end subroutine
#:enddef

#:def rot_mfi(mfi_name,f77_name,pfxs)
#:set A = pfxs[0]
#:set B = A if len(pfxs) == 1 else pfxs[1]
!> Given two vectors x and y,
!> each vector element of these vectors is replaced as follows:
!>```fortran
#:if type(A) == real(A)
!> xi = c*xi + s*yi
!> yi = c*yi - s*xi
#:elif type(A) == complex(A)
!> xi = c*xi + s*yi
!> yi = c*yi - conj(s)*xi
#:endif
!>```
pure subroutine ${mfi_name}$(x, y, c, s, incx, incy)
@:args(${type(A)}$, inout, x(:), y(:))
@:args(${real(A)}$, in, c)
@:args(${type(B)}$, in, s)
@:optional(integer, in, incx, incy)
    integer :: n
@:defaults(incx=1, incy=1)
    n = size(x)
    call ${f77_name}$(n,x,local_incx,y,local_incy,c,s)
end subroutine
#:enddef

$:f77_original('?rot', DEFAULT_TYPES + mix(COMPLEX_TYPES,REAL_TYPES), rot_f77)
$:mfi_interface('?rot', DEFAULT_TYPES + mix(COMPLEX_TYPES,REAL_TYPES))
$:mfi_implement('?rot', DEFAULT_TYPES + mix(COMPLEX_TYPES,REAL_TYPES), rot_mfi)


#:endmute
