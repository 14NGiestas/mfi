#:mute

#:set REAL    = 'real(wp)'
#:set COMPLEX = 'complex(wp)'
#:set PREFIX = { &
    's': { 'type': 'real(wp)',    'wp': 'REAL32', 'c_kind': 'c_float'  }, &
    'd': { 'type': 'real(wp)',    'wp': 'REAL64', 'c_kind': 'c_double' }, &
    'c': { 'type': 'complex(wp)', 'wp': 'REAL32', 'c_kind': 'c_float'  }, &
    'z': { 'type': 'complex(wp)', 'wp': 'REAL64', 'c_kind': 'c_double' }, &
}

#:set ERROR = lambda pfx: { 'type': f'error: {pfx}', 'wp' : f'error: {pfx}' }

#:set mix    = lambda l, r: list(lp + rp for lp, rp in zip(l,r))
#:set split  = lambda pfx: list(pfx) if len(pfx) > 1 else pfx
#:set get_types = lambda pfxs: (pfxs[0], pfxs[0] if len(pfxs) == 1 else pfxs[1])
#:set get    = lambda pfx,what: PREFIX.get(pfx).get(what)
#:set prefix = lambda pfx, name: name.replace('?',pfx)
#:set kind   = lambda pfx: get(pfx,'wp')
#:set c_kind = lambda pfx: get(pfx,'c_kind')
#:set type   = lambda pfx: get(pfx,'type').replace('wp',kind(pfx))
#:set c_type = lambda pfx: get(pfx,'type').replace('wp',c_kind(pfx))
#:set real      = lambda pfx: REAL.replace('wp',kind(pfx))
#:set complex   = lambda pfx: COMPLEX.replace('wp',kind(pfx))

#:set functions = lambda gen_name, pfxs: map(lambda pfx: prefix(pfx,gen_name), pfxs)

#:set SINGLE_TYPES  = ['s','c']
#:set DOUBLE_TYPES  = ['d','z']
#:set REAL_TYPES    = ['s','d']
#:set COMPLEX_TYPES = ['c','z']
#:set DEFAULT_TYPES = REAL_TYPES + COMPLEX_TYPES

#:set MIX_REAL_COMPLEX  = mix(REAL_TYPES,COMPLEX_TYPES)
#:set MIX_COMPLEX_REAL  = mix(COMPLEX_TYPES,REAL_TYPES)
#:set MIX_SINGLE_DOUBLE = mix(SINGLE_TYPES,DOUBLE_TYPES)
#:set MIX_DOUBLE_SINGLE = mix(DOUBLE_TYPES,SINGLE_TYPES)

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
    ${name}$ = cmplx(re,im, kind=${type.replace('complex(', '').replace(')', '')}$)
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

#:def interface(functions, procedure='procedure', name='', prefix='')
interface ${prefix}$${name}$
    #:for function_name in functions
    ${procedure}$ :: ${function_name}$
    #:endfor
end interface
#:enddef

#! Interfaces for the original f77 routines
#! code must implement a routine interface
#:def f77_original(generic_name, prefixes, code)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set f90 =          prefix('',generic_name)
#:set f77 = [prefix(pfx,generic_name) for pfx in prefixes]
!> Generic old style interface for ${prefix('',generic_name).upper()}$.
!> Supports ${', '.join(prefixes)}$.
!> See also: [[${mfi}$]], ${'[[' + ']], [['.join(f77) + ']]'}$.
interface f77_${prefix('',generic_name)}$
#:for pfx in prefixes
#:set name = prefix(pfx,generic_name)
#:set pfxs = list(map(split,pfx))
!> Original interface for ${name.upper()}$
!> See also: [[${mfi}$]], [[${f90}$]].
$:code(name,pfxs)
#:endfor
end interface
#:enddef


#! Define a common interface with the original f77 interfaces
#! So you can call the original function without the prefix
#:def f77_improved(generic_name, prefixes)
$:interface(functions(generic_name, prefixes), name=f"f77_{prefix('',generic_name)}")
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
#:set f77 = ['f77_' + prefix('',generic_name) + ':' + prefix(pfx,generic_name) for pfx in prefixes]
!> Generic modern interface for ${prefix('',generic_name).upper()}$.
!> Supports ${', '.join(prefixes)}$.
!> See also:
!> ${'[[' + ']], [['.join(f77) + ']]'}$.
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
#:set f77_name  = 'f77_' + prefix('',generic_name)
#:set pfxs      = list(map(split,pfx))
#:set fun       = prefix('',generic_name)
!> Modern interface for [[f77_${fun}$:${f77_name}$]].
!> See also: [[mfi_${fun}$]], [[f77_${fun}$]].
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

#! GPU variant — calls the original macro with keyword args + suffix
#:def test_implement_gpu(generic_name, prefixes, code)
#:for pfx in prefixes
#:set f77  =          prefix(pfx,generic_name)
#:set f90 = 'f77_' + prefix('',generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set pfxs    = list(map(split,pfx))
$:code(f77=f77, f90=f90, mfi=mfi, pfxs=pfxs, suffix='_gpu')
#:endfor
#:enddef

#:def test_run(generic_name, prefixes)
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
@:timeit("testing ${mfi}$ (CPU) against ${f77}$", { call ${f77}$ })
#:endfor
#:enddef

#:def test_run_lapack(generic_name, prefixes)
#:set f77_names = ', '.join([prefix(pfx, generic_name) for pfx in prefixes])
use f77_lapack, only: ${f77_names}$
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
@:timeit("testing ${mfi}$ (CPU) against ${f77}$", { call ${f77}$ })
#:endfor
#:enddef

#:def test_run_gpu(generic_name, prefixes)
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
@:timeit("testing ${mfi}$ (GPU) against ${f77}$", { call ${f77 + '_gpu'}$ })
#:endfor
#:enddef

#:def test_run_lapack_gpu(generic_name, prefixes)
#:set f77_names = ', '.join([prefix(pfx, generic_name) + '_gpu' for pfx in prefixes])
use f77_lapack, only: ${f77_names}$
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
@:timeit("testing ${mfi}$ (GPU) against ${f77}$", { call ${f77 + '_gpu'}$ })
#:endfor
#:enddef

#! Run code on GPU mode (only when compiled with MFI_CUBLAS).
#! If MFI_CUBLAS is not defined, expands to just the code (no-op).
#! Usage: @:mfi_gpu({ call mfi_gemm(A, B, C) })
#! Requires: use mfi_blas (for mfi_force_gpu/cpu when cublas is available)
#:def mfi_gpu(code)
#if defined(MFI_CUBLAS)
    call mfi_force_gpu()
#endif
    $:code
#if defined(MFI_CUBLAS)
    call mfi_force_cpu()
#endif
#:enddef

#! Run code on CPU mode — explicit marker for readability.
#! Always expands to just the code (CPU is the default).
#:def mfi_cpu(code)
    $:code
#:enddef

#! Initialize random seed from MFI_TEST_SEED env var (or default=42).
#! Usage: $:test_seed()
#:def test_seed()
block
    integer, parameter :: seed_size = 8
    integer :: seed_arr(seed_size)
    integer :: env_seed
    integer :: seed_stat
    integer :: ii
    character(64) :: env_val
    call get_environment_variable('MFI_TEST_SEED', value=env_val, status=seed_stat)
    if (seed_stat == 0 .and. len_trim(env_val) > 0) then
        read(env_val, '(I10)', iostat=seed_stat) env_seed
    end if
    if (seed_stat /= 0) env_seed = 42
    do ii = 0, seed_size - 1
        seed_arr(ii + 1) = mod(env_seed * (ii + 1), 2147483647)
    end do
    call random_seed(put=seed_arr)
end block
#:enddef

#! Assert two 2D arrays are close within type-appropriate tolerance.
#! Usage: $:assert_close('A', 'A_rf', 'label', wp)
#:def assert_close(actual, expected, label, wp)
#:set knd = kind(wp)
#:if wp in ['s','d']
call assert(maxval(abs(${actual}$ - ${expected}$)) < sqrt(epsilon(1.0_${knd}$)), "${label}$: mismatch")
#:else
call assert(maxval(abs(${actual}$ - ${expected}$)) < 2.0 * sqrt(epsilon(1.0_${knd}$)), "${label}$: mismatch")
#:endif
#:enddef

#! Assert two 1D arrays are close within type-appropriate tolerance.
#! Usage: $:assert_close_1d('x', 'x_rf', 'label', wp)
#:def assert_close_1d(actual, expected, label, wp)
#:set knd = kind(wp)
#:if wp in ['s','d']
call assert(maxval(abs(${actual}$ - ${expected}$)) < sqrt(epsilon(1.0_${knd}$)), "${label}$: mismatch")
#:else
call assert(maxval(abs(${actual}$ - ${expected}$)) < 2.0 * sqrt(epsilon(1.0_${knd}$)), "${label}$: mismatch")
#:endif
#:enddef

#:endmute
