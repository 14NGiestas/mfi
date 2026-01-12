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
    #:else
    ! Issue warning when used without proper compilation
    print *, 'WARNING: mfi_execution_init() called but compiled without CUBLAS support'
    #:endif
end subroutine

!> Sets execution mode to use GPU with CUBLAS (only effective when compiled with support)
subroutine mfi_force_gpu()
    MFI_USE_CUBLAS_PREV_STATE = MFI_USE_CUBLAS
    #:if defined('MFI_USE_CUBLAS')
    MFI_USE_CUBLAS = 1
    #:else
    print *, 'WARNING: mfi_force_gpu() called but compiled without CUBLAS support'
    ! Don't actually enable GPU mode
    #:endif
end subroutine

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

