#:set REAL_TYPES=['s','d']
#:set COMPLEX_TYPES=['c','z']
#:set DEFAULT_TYPES=REAL_TYPES+COMPLEX_TYPES

#:set PREFIX_TO_TYPE={    &
    's':   'real(wp)',    &
    'd':   'real(wp)',    &
    'c':   'complex(wp)', &
    'z':   'complex(wp)', &
}

#:set PREFIX_TO_KIND={&
    's':   'REAL32',  &
    'd':   'REAL64',  &
    'c':   'REAL32',  &
    'z':   'REAL64',  &
}

#! Defines a optional variable, creating local corresponding variable by default
#:def optional(dtype, *args)
#:for variable in args
    ${dtype}$, intent(in), optional :: ${variable}$
    ${dtype}$ :: local_${variable}$
#:endfor
#:enddef

#! Define a local variable
#:def localvars(dtype, *args)
#:for variable in args
    ${dtype}$ :: ${variable}$
#:endfor
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
#:def mfi_implement(name, supports, code)
#:for PREFIX in supports
#:set NAME = name.replace('?',PREFIX)
#:set TYPE = PREFIX_TO_TYPE.get(PREFIX,None)
#:set KIND = PREFIX_TO_KIND.get(PREFIX,None)
$:code(NAME,TYPE,KIND)
#:endfor
#:enddef

#:def mfi_interface(name, types)
interface ${name.replace('?','')}$
    #:for T in types
    module procedure mfi_${name.replace('?',T)}$
    #:endfor
end interface
#:enddef
