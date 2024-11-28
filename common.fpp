#:mute
#:set REAL_TYPE='real(wp)'
#:set COMPLEX_TYPE='complex(wp)'
#:set REAL_TYPES=['s','d']
#:set COMPLEX_TYPES=['c','z']
#:set DEFAULT_TYPES=REAL_TYPES+COMPLEX_TYPES

#:set PREFIX_TO_TYPE={   &
    's':   REAL_TYPE,    &
    'd':   REAL_TYPE,    &
    'c':   COMPLEX_TYPE, &
    'z':   COMPLEX_TYPE, &
}

#:set PREFIX_TO_KIND={&
    's':   'REAL32',  &
    'd':   'REAL64',  &
    'c':   'REAL32',  &
    'z':   'REAL64',  &
}

#:set TYPE_AND_KIND_TO_PREFIX = { &
    'real(REAL32)':    's',       &
    'real(REAL64)':    'd',       &
    'complex(REAL32)': 'c',       &
    'complex(REAL64)': 'z',       &
}

#! Defines a optional variable, creating local corresponding variable by default
#:def optional(dtype, intent, *args)
#:for variable in args
    ${dtype}$, intent(${intent}$), optional :: ${variable}$
    ${dtype}$ :: local_${variable}$
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

#! Handles default values of the optional
#:def defaults(**kwargs)
#:for variable, default in kwargs.items()
    if (present(${variable}$)) then
        local_${variable}$ = ${variable}$
    else
        local_${variable}$ = ${default}$
    end if
#:endfor
#:enddef

#! Handles the input/output arguments
#:def args(dtype, intent, *args)
#:for variable in args
    ${dtype}$, intent(${intent}$) :: ${variable}$
#:endfor
#:enddef

#! Handles parameters (usage: working precision)
#:def parameter(dtype, **kwargs)
#:for variable, value in kwargs.items()
    ${dtype}$, parameter :: ${variable}$ = ${value}$
#:endfor
#:enddef

#! Handles the implementation of the modern interface to each supported type and kind
#:def mfi_implement(name, supports, code, f=lambda x: x)
#:for PREFIX in supports
#:set MFI_NAME = "mfi_" + name.replace('?',f(PREFIX))
#:set F77_NAME = name.replace('?',f(PREFIX))
#:set TYPE = PREFIX_TO_TYPE.get(PREFIX,None)
#:set KIND = PREFIX_TO_KIND.get(PREFIX,None)
$:code(MFI_NAME,F77_NAME,TYPE,KIND,PREFIX)
#:endfor
#:enddef

#! Define mfi interfaces to implemented routines
#:def mfi_interface(name, types, f=lambda x: x)
interface mfi_${name.replace('?','')}$
    #:for T in types
    module procedure mfi_${name.replace('?',f(T))}$
    #:endfor
end interface
#:enddef

#! Define f77 interfaces to implemented routines
#:def f77_interface_improved(name, types, f=lambda x: x)
interface f77_${name.replace('?','')}$
    #:for T in types
    procedure :: ${name.replace('?',f(T))}$
    #:endfor
end interface
#:enddef

#! Define a f77 interfaces to the external blas/lapack library
#:def f77_interface(name, supports, code, f=lambda x: x, improved_f77=True)

interface
#:for PREFIX in supports
#:set NAME = name.replace('?',f(PREFIX))
#:set TYPE = PREFIX_TO_TYPE.get(PREFIX,None)
#:set KIND = PREFIX_TO_KIND.get(PREFIX,None)
$:code(NAME,TYPE,KIND,PREFIX)
#:endfor
end interface

#:if improved_f77
$:f77_interface_improved(name, supports, f=f)
#:endif

#:enddef

#! Implements a f77 function / extension
#:def f77_implement(name, supports, code)
#:for PREFIX in supports
#:set NAME = name.replace('?',PREFIX)
#:set TYPE = PREFIX_TO_TYPE.get(PREFIX,None)
#:set KIND = PREFIX_TO_KIND.get(PREFIX,None)
$:code(NAME,TYPE,KIND)
#:endfor
#:enddef

#! Implements a test
#:def test_implement(name, supports, code, f=lambda x: x)
#:for PREFIX in supports
#:set ORIGINAL = name.replace('?',f(PREFIX))
#:set IMPROVED = "f77_" + name.replace('?','')
#:set MODERN   = "mfi_" + name.replace('?','')
#:set TYPE = PREFIX_TO_TYPE.get(PREFIX,None)
#:set KIND = PREFIX_TO_KIND.get(PREFIX,None)
$:code(ORIGINAL,IMPROVED,MODERN,TYPE,KIND,PREFIX)
#:endfor
#:enddef

#! Call the subroutine test
#:def test_run(name, supports, f=lambda x: x)
#:for PREFIX in supports
#:set ORIGINAL = name.replace('?',f(PREFIX))
print*, 'calling ${ORIGINAL}$'
call test_${ORIGINAL}$
#:endfor
#:enddef

#:def timeit(message, code)
block
real :: t1, t2
call cpu_time(t1)
$:code
call cpu_time(t2)
print '(A,G0)', ${message}$, t2-t1
end block
#:enddef
#:endmute
