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
#:endmute
