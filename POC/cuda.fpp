#:mute
#:def copyin(*varlist)
#:for var in varlist
    call cublas_alloc(size(${var}$),wp,device_${var}$)
#:endfor
#:for var in varlist
    call cublas_set_matrix(size(${var}$,1),size(${var}$,2),wp,${var}$,max(1,size(${var}$,1)),device_${var}$,max(1,size(${var}$,1)))
#:endfor
#:enddef

#:def delete(*varlist)
#:for var in varlist
    call cublas_free(device_${var}$)
#:endfor
#:enddef

#:def copyout(*varlist)
#:for var in varlist
    call cublas_get_matrix(size(${var}$,1),size(${var}$,2),wp,device_${var}$,max(1,size(${var}$,1)),${var}$,max(1,size(${var}$,1)))
#:endfor
#:enddef

#:def gpu_interface(name, supports, code)
interface gpu_${name.replace('?','')}$
#:for PREFIX in supports
#:set NAME = f"cublas{name.replace('?',PREFIX.upper())}"
#:set TYPE = PREFIX_TO_TYPE.get(PREFIX,None)
#:set KIND = PREFIX_TO_C_KIND.get(PREFIX,None)
$:code(NAME,TYPE,KIND)
#:endfor
end interface
#:enddef

#:endmute
