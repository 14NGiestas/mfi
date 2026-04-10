#:mute
#:include "common.fpp"
#:def mfi_extensions_interfaces()
! BLAS level 1 - Utils / Extensions
$:mfi_interface('i?amax', DEFAULT_TYPES)
$:mfi_interface('i?amin', DEFAULT_TYPES)

! Global variables for execution mode control - always defined
integer, save :: MFI_USE_CUBLAS = 0
integer, save :: MFI_USE_CUBLAS_PREV_STATE = 0
#:enddef

#:def mfi_extensions_implement()
! BLAS level 1 - Utils / Extensions
$:mfi_implement('i?amax', DEFAULT_TYPES, iamin_iamax)
$:mfi_implement('i?amin', DEFAULT_TYPES, iamin_iamax)

!> Initialize execution mode from environment variables
subroutine mfi_execution_init()
    integer :: info
    character(len=10) :: env_mfi_use_cublas
    integer :: env_len_cublas

    ! Save previous states
    MFI_USE_CUBLAS_PREV_STATE = MFI_USE_CUBLAS

    ! Initialize to defaults
    MFI_USE_CUBLAS = 0

    ! Only check environment if extensions are enabled
    #:if defined('MFI_USE_CUBLAS')
    ! Check CUBLAS setting
    call get_environment_variable("MFI_USE_CUBLAS",env_mfi_use_cublas,env_len_cublas)
    if (env_len_cublas > 0) then
        read(env_mfi_use_cublas(1:env_len_cublas),*,iostat=info) MFI_USE_CUBLAS
    end if
    ! Initialize cuBLAS handle if GPU mode is enabled
    if (MFI_USE_CUBLAS == 1) then
        call mfi_cublas_handle_ensure()
    end if
    #:else
    ! Issue warning when used without proper compilation
    print *, 'WARNING: mfi_execution_init() called but compiled without CUBLAS support'
    #:endif
end subroutine

!> Ensure cuBLAS handle is created (called internally)
#:if defined('MFI_USE_CUBLAS')
subroutine mfi_cublas_handle_ensure()
    integer(c_int) :: stat
    if (.not. c_associated(mfi_cublas_handle)) then
        stat = cublasCreate(mfi_cublas_handle)
        if (stat /= 0) error stop 'cublasCreate_v2 failed - check CUDA driver version'
        stat = cublasSetPointerMode(mfi_cublas_handle, CUBLAS_POINTER_MODE_HOST)
        if (stat /= 0) error stop 'cublasSetPointerMode_v2 failed'
    end if
end subroutine
#:endif

!> Sets execution mode to use GPU with CUBLAS (only effective when compiled with support)
subroutine mfi_force_gpu()
    MFI_USE_CUBLAS_PREV_STATE = MFI_USE_CUBLAS
    #:if defined('MFI_USE_CUBLAS')
    MFI_USE_CUBLAS = 1
    call mfi_cublas_handle_ensure()
    #:else
    print *, 'WARNING: mfi_force_gpu() called but compiled without CUBLAS support'
    ! Don't actually enable GPU mode
    #:endif
end subroutine

!> Finalize cuBLAS resources
#:if defined('MFI_USE_CUBLAS')
subroutine mfi_cublas_finalize()
    integer(c_int) :: stat
    if (c_associated(mfi_cublas_handle)) then
        stat = cublasDestroy(mfi_cublas_handle)
        mfi_cublas_handle = c_null_ptr
    end if
end subroutine
#:endif

!> Report cuBLAS error (called from pure wrappers)
#:if defined('MFI_USE_CUBLAS')
pure subroutine mfi_cublas_error(stat, name)
    integer(c_int), value, intent(in) :: stat
    character(*), intent(in) :: name
    character(len=120) :: msg
    write(msg, '(A, I0, A)') 'cuBLAS error: ', stat, ' (' // trim(name) // ')'
    error stop msg
end subroutine
#:endif

!> Sets execution mode to use traditional CPU
subroutine mfi_force_cpu()
    MFI_USE_CUBLAS_PREV_STATE = MFI_USE_CUBLAS
    MFI_USE_CUBLAS = 0
end subroutine

!> Restores the previous execution mode state
subroutine mfi_execution_restore()
    MFI_USE_CUBLAS = MFI_USE_CUBLAS_PREV_STATE
end subroutine

!> Returns current execution mode as string
function mfi_get_execution_mode() result(mode_str)
    character(len=20) :: mode_str
    #:if defined('MFI_USE_CUBLAS')
    if (MFI_USE_CUBLAS == 1) then
        mode_str = "GPU-CUBLAS"
    else
        mode_str = "CPU-DEFAULT"
    end if
    #:else
    mode_str = "CPU-ONLY"
    #:endif
end function
#:enddef
#:endmute

