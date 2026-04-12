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
    integer :: t_n, t_i, t_stat
    real :: t_t1, t_t2, t_sum, t_sum2, t_tmin, t_tmax, t_sig
    real, allocatable :: t_dt(:)
    character(16) :: t_mu, t_ms, t_mn, t_mx
    character(32) :: t_env
    call get_environment_variable('MFI_TEST_SAMPLES', value=t_env, status=t_stat)
    if (t_stat == 0 .and. len_trim(t_env) > 0) then
        read(t_env, '(I10)', iostat=t_stat) t_n
    else
        t_n = 3
    end if
    if (t_n < 1) t_n = 1
    allocate(t_dt(t_n))
    t_tmin = huge(1.0)
    t_tmax = -huge(1.0)
    t_sum  = 0.0
    t_sum2 = 0.0
    do t_i = 1, t_n
        call cpu_time(t_t1)
        $:code
        call cpu_time(t_t2)
        t_dt(t_i) = t_t2 - t_t1
        t_tmin = min(t_tmin, t_dt(t_i))
        t_tmax = max(t_tmax, t_dt(t_i))
        t_sum  = t_sum  + t_dt(t_i)
        t_sum2 = t_sum2 + t_dt(t_i)**2
    end do
    deallocate(t_dt)
    t_sig = sqrt(max(t_sum2/t_n - (t_sum/t_n)**2, 0.0))
    call fmt_time(t_sum/t_n, t_mu)
    call fmt_time(t_sig, t_ms)
    call fmt_time(t_tmin, t_mn)
    call fmt_time(t_tmax, t_mx)
    print '(A,"  μ=",A16," σ=",A16," min=",A16," max=",A16,"  (",I0," runs)")', &
        ${message + ' ' * max(0, 42 - len(str(message).replace(chr(27)+'[', '<<')))}$, t_mu, t_ms, t_mn, t_mx, t_n
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

#! ANSI color codes (Python chr() for Fypp-time evaluation)
#:set _GREEN    = chr(27)+'[32m'
#:set _BLUE     = chr(27)+'[34m'
#:set _RED      = chr(27)+'[31m'
#:set _YELLOW   = chr(27)+'[33m'
#:set _RESET    = chr(27)+'[0m'

#! Label padding — adds trailing spaces to align stats columns.
#:set _pad = lambda s: s + ' ' * max(0, 40 - len(s.replace(chr(27)+'[', '<<')))

#:def test_run(generic_name, prefixes)
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set pfx_c = _RED if pfx == 's' else (_GREEN if pfx == 'd' else (_BLUE if pfx == 'c' else _YELLOW))
#:set f77_s = pfx_c + f77 + _RESET
@:timeit("testing ${pfx_c}$${pfx}$${_RESET}$ ${mfi}$ (${_BLUE}$CPU${_RESET}$) against ${f77_s}$", { call test_${f77}$ })
#:endfor
#:enddef

#:def test_run_lapack(generic_name, prefixes)
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set pfx_c = _RED if pfx == 's' else (_GREEN if pfx == 'd' else (_BLUE if pfx == 'c' else _YELLOW))
#:set f77_s = pfx_c + f77 + _RESET
@:timeit("testing ${pfx_c}$${pfx}$${_RESET}$ ${mfi}$ (${_BLUE}$CPU${_RESET}$) against ${f77_s}$", { call test_${f77}$ })
#:endfor
#:enddef

#:def test_run_lapack_gpu(generic_name, prefixes)
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set pfx_c = _RED if pfx == 's' else (_GREEN if pfx == 'd' else (_BLUE if pfx == 'c' else _YELLOW))
#:set f77_s = pfx_c + f77 + _RESET
@:timeit("testing ${pfx_c}$${pfx}$${_RESET}$ ${mfi}$ (${_GREEN}$GPU${_RESET}$) against ${f77_s}$", { call test_${f77}$_gpu })
#:endfor
#:enddef

#:def test_run_gpu(generic_name, prefixes)
#:for pfx in prefixes
#:set f77 =          prefix(pfx,generic_name)
#:set mfi = 'mfi_' + prefix('',generic_name)
#:set pfx_c = _RED if pfx == 's' else (_GREEN if pfx == 'd' else (_BLUE if pfx == 'c' else _YELLOW))
#:set f77_s = pfx_c + f77 + _RESET
@:timeit("testing ${pfx_c}$${pfx}$${_RESET}$ ${mfi}$ (${_GREEN}$GPU${_RESET}$) against ${f77_s}$", { call test_${f77}$_gpu })
#:endfor
#:enddef

#! GPU test wrapper — activates GPU, runs code, resets to CPU.
#! Only expands to GPU switches when suffix is non-empty and MFI_CUBLAS is defined.
#! Usage: @:on_gpu({ call mfi_gemm(A, B, C) }, suffix)
#:def on_gpu(code, suffix='')
#:if suffix
#if defined(MFI_CUBLAS)
    call mfi_force_gpu()
#endif
    $:code
#if defined(MFI_CUBLAS)
    call mfi_force_cpu()
#endif
#:else
    $:code
#:endif
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

#:def fmt_time_fn()

  subroutine fmt_time(t, out)
      real, intent(in) :: t
      character(*), intent(out) :: out
      if (t < 1.0e-3) then
          write(out, '(F12.3,"µs")') t * 1.0e6
      else if (t < 1.0) then
          write(out, '(F12.3,"ms")') t * 1.0e3
      else
          write(out, '(F12.3,"s ")') t
      end if
  end subroutine fmt_time
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

#! Test problem size — reads MFI_TEST_ELEMENTS env var at Fypp time (compile-time default).
#! For runtime sizing, use $:test_get_size_fn() in the contains block and call
#!   call test_get_1d(500000, N)   ! 1D array length from env var or default
#!   call test_get_2d(500000, N)   ! NxN from env var or default (N=sqrt(elements))
#! Default: 500K elements (N=707 for 2D matrices).
#:set _total_elements = int(os.environ.get('MFI_TEST_ELEMENTS', '500000'))
#:set N1D = lambda n=_total_elements: n
#:set N2D = lambda n=_total_elements: int(n ** 0.5)

#! Runtime test size — use inside subroutines with allocatable arrays.
#! Usage:
#!   call test_get_1d(10000000, N)   ! sets N from env var or default
#!   call test_get_2d(10000000, N)   ! sets N=sqrt(elements) from env var or default
#:def test_get_size_fn()

  subroutine test_get_1d(default, n)
      integer, intent(in) :: default
      integer, intent(out) :: n
      integer :: tgs_stat
      character(32) :: tgs_env
      call get_environment_variable('MFI_TEST_ELEMENTS', value=tgs_env, status=tgs_stat)
      if (tgs_stat == 0 .and. len_trim(tgs_env) > 0) then
          read(tgs_env, '(I10)', iostat=tgs_stat) n
      else
          n = default
      end if
      if (n < 1) n = default
  end subroutine test_get_1d

  subroutine test_get_2d(default, n)
      integer, intent(in) :: default
      integer, intent(out) :: n
      integer :: tgs_stat, tgs_val
      character(32) :: tgs_env
      call get_environment_variable('MFI_TEST_ELEMENTS', value=tgs_env, status=tgs_stat)
      if (tgs_stat == 0 .and. len_trim(tgs_env) > 0) then
          read(tgs_env, '(I10)', iostat=tgs_stat) tgs_val
          n = int(sqrt(real(tgs_val)))
      else
          n = int(sqrt(real(default)))
      end if
      if (n < 4) n = 4
  end subroutine test_get_2d
#:enddef

#! Runtime test size functions — read MFI_TEST_ELEMENTS env var at runtime.
#! Like test_seed() but for array dimensions. Use with allocatable arrays.
#! Usage: $:test_size_fn() in the contains block, then:
#!        integer :: N = test_size_1d(2000)   ! 1D vector length
#!        integer :: N = test_size_2d(2000)   ! sqrt(elements) for NxN matrix
#:def test_size_fn()

  function test_size_1d(default) result(n)
      integer, intent(in) :: default
      integer :: n
      integer :: env_stat
      character(32) :: env_val
      call get_environment_variable('MFI_TEST_ELEMENTS', value=env_val, status=env_stat)
      if (env_stat == 0 .and. len_trim(env_val) > 0) then
          read(env_val, '(I10)', iostat=env_stat) n
      end if
      if (env_stat /= 0) n = default
  end function test_size_1d

  function test_size_2d(default) result(n)
      integer, intent(in) :: default
      integer :: n
      integer :: env_stat
      character(32) :: env_val
      call get_environment_variable('MFI_TEST_ELEMENTS', value=env_val, status=env_stat)
      if (env_stat == 0 .and. len_trim(env_val) > 0) then
          read(env_val, '(I10)', iostat=env_stat) n
          n = int(sqrt(real(n)))
      end if
      if (env_stat /= 0) n = int(sqrt(real(default)))
  end function test_size_2d

  function test_samples(default) result(ns)
      integer, intent(in) :: default
      integer :: ns
      integer :: env_stat
      character(32) :: env_val
      call get_environment_variable('MFI_TEST_SAMPLES', value=env_val, status=env_stat)
      if (env_stat == 0 .and. len_trim(env_val) > 0) then
          read(env_val, '(I10)', iostat=env_stat) ns
      end if
      if (env_stat /= 0) ns = default
  end function test_samples
#:enddef

#! Assert two 2D arrays are close within type-appropriate tolerance.
#! Usage: @:assert_close(A, A_rf, 'label', wp)
#!   @:assert_close(A, A_rf, 'label', wp, 10.0)  — custom tolerance multiplier
#:def assert_close(actual, expected, label, wp, tol=None)
#:set knd = kind(wp)
#:set eps = 'epsilon(1.0_' + knd + ')'
#:if wp in ['s','d']
#:if tol is None
    call assert(maxval(abs(${actual}$ - ${expected}$)) < sqrt(${eps}$), ${label}$ // ": mismatch")
#:else
    call assert(maxval(abs(${actual}$ - ${expected}$)) < ${tol}$ * sqrt(${eps}$), ${label}$ // ": mismatch")
#:endif
#:else
#:if tol is None
    call assert(maxval(abs(${actual}$ - ${expected}$)) < 2.0 * sqrt(${eps}$), ${label}$ // ": mismatch")
#:else
    call assert(maxval(abs(${actual}$ - ${expected}$)) < ${tol}$ * 2.0 * sqrt(${eps}$), ${label}$ // ": mismatch")
#:endif
#:endif
#:enddef

#! Assert two 1D arrays are close within type-appropriate tolerance.
#! Usage: @:assert_close_1d(x, x_rf, 'label', wp)
#!   @:assert_close_1d(x, x_rf, 'label', wp, 10.0)  — custom tolerance multiplier
#:def assert_close_1d(actual, expected, label, wp, tol=None)
#:set knd = kind(wp)
#:set eps = 'epsilon(1.0_' + knd + ')'
#:if wp in ['s','d']
#:if tol is None
    call assert(maxval(abs(${actual}$ - ${expected}$)) < sqrt(${eps}$), ${label}$ // ": mismatch")
#:else
    call assert(maxval(abs(${actual}$ - ${expected}$)) < ${tol}$ * sqrt(${eps}$), ${label}$ // ": mismatch")
#:endif
#:else
#:if tol is None
    call assert(maxval(abs(${actual}$ - ${expected}$)) < 2.0 * sqrt(${eps}$), ${label}$ // ": mismatch")
#:else
    call assert(maxval(abs(${actual}$ - ${expected}$)) < ${tol}$ * 2.0 * sqrt(${eps}$), ${label}$ // ": mismatch")
#:endif
#:endif
#:enddef

#:endmute
