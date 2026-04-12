#:mute
#:include "common.fpp"
#:def mfi_extensions_interfaces()
! BLAS level 1 - Utils / Extensions
$:mfi_interface('i?amax', DEFAULT_TYPES)
$:mfi_interface('i?amin', DEFAULT_TYPES)

!> Public API — always available (stubs when no cuBLAS)
public :: mfi_iamax, mfi_isamax, mfi_idamax, mfi_icamax, mfi_izamax
public :: mfi_iamin, mfi_isamin, mfi_idamin, mfi_icamin, mfi_izamin
public :: mfi_force_gpu, mfi_force_cpu
public :: mfi_cublas_finalize, mfi_cublas_is_active

!> Internal state — not exposed to users
integer, save :: MFI_USE_CUBLAS = 0
logical, save :: mfi_cublas_env_checked = .false.
logical, save :: mfi_cublas_global_initialized = .false.
integer, save :: mfi_cublas_handle_count = 0
private :: MFI_USE_CUBLAS
private :: mfi_cublas_env_checked
private :: mfi_cublas_global_initialized
private :: mfi_cublas_handle_count
#:enddef

#:def mfi_extensions_implement()
! BLAS level 1 - Utils / Extensions
$:mfi_implement('i?amax', DEFAULT_TYPES, iamin_iamax)
$:mfi_implement('i?amin', DEFAULT_TYPES, iamin_iamax)

!> Check environment and initialize cuBLAS on first GPU call (lazy init)
#if defined(MFI_CUBLAS)
subroutine mfi_cublas_lazy_init()
    integer :: info
    integer(c_int) :: stat
    character(len=16) :: env_val
    integer :: env_len

    if (mfi_cublas_global_initialized) return

    if (.not. mfi_cublas_env_checked) then
        mfi_cublas_env_checked = .true.
        call get_environment_variable("MFI_USE_CUBLAS", env_val, env_len)
        if (env_len > 0) then
            read(env_val(1:env_len), *, iostat=info) MFI_USE_CUBLAS
        end if
    end if

    if (MFI_USE_CUBLAS == 1 .and. mfi_cublas_handle_count == 0) then
        call get_environment_variable("OMP_NUM_THREADS", env_val, env_len)
        if (env_len > 0) then
            read(env_val(1:env_len), *, iostat=info) mfi_cublas_handle_count
            if (info /= 0 .or. mfi_cublas_handle_count < 1) mfi_cublas_handle_count = 1
        else
            mfi_cublas_handle_count = 1
        end if
        call mfi_cublas_set_count(mfi_cublas_handle_count)
        call mfi_cublas_init_handles(stat)
        if (stat /= 0) error stop 'mfi_cublas_init_handles failed - check CUDA driver version'
    end if

    mfi_cublas_global_initialized = .true.
end subroutine
#endif

!> Returns .true. if GPU execution is active (pure, checks flag only)
pure logical function mfi_cublas_is_active() result(active)
    active = (MFI_USE_CUBLAS == 1)
end function

!> Get the cuBLAS handle for the current thread (pure, delegates to C wrapper)
pure function mfi_cublas_handle_get() result(handle)
#if defined(MFI_CUBLAS)
    use mfi_blas_cublas, only: mfi_cublas_get_thread_handle
    type(c_ptr) :: handle
    integer(c_int) :: stat
    call mfi_cublas_get_thread_handle(handle, stat)
    if (stat /= 0) error stop 'mfi: thread ID out of range. '// &
               'Set OMP_NUM_THREADS before running with MFI_USE_CUBLAS=1.'
#else
    type(c_ptr) :: handle
    handle = c_null_ptr
#endif
end function

!> Switch to GPU mode (no-op when compiled without cuBLAS)
#if defined(MFI_CUBLAS)
subroutine mfi_force_gpu()
    MFI_USE_CUBLAS = 1
    mfi_cublas_global_initialized = .false.
    mfi_cublas_env_checked = .true.
    call mfi_cublas_lazy_init()
end subroutine
#else
subroutine mfi_force_gpu()
end subroutine
#endif

!> Switch to CPU mode (no-op when compiled without cuBLAS)
#if defined(MFI_CUBLAS)
subroutine mfi_force_cpu()
    MFI_USE_CUBLAS = 0
    mfi_cublas_global_initialized = .false.
    mfi_cublas_env_checked = .true.
end subroutine
#else
subroutine mfi_force_cpu()
    ! No-op: already on CPU
end subroutine
#endif

!> Finalize all cuBLAS handles (no-op when compiled without cuBLAS)
#if defined(MFI_CUBLAS)
subroutine mfi_cublas_finalize()
    use mfi_blas_cublas, only: mfi_cublas_finalize_all
    integer(c_int) :: stat
    call mfi_cublas_finalize_all(stat)
    mfi_cublas_handle_count = 0
    MFI_USE_CUBLAS = 0
    mfi_cublas_global_initialized = .false.
end subroutine
#else
subroutine mfi_cublas_finalize()
    ! No-op: no cuBLAS resources to release
end subroutine
#endif

!> Report cuBLAS error (called from pure wrappers)
#if defined(MFI_CUBLAS)
pure subroutine mfi_cublas_error(stat, name)
    integer(c_int), value, intent(in) :: stat
    character(*), intent(in) :: name
    character(len=120) :: msg
    character(len=40) :: stat_str

    select case(stat)
        case (CUBLAS_STATUS_SUCCESS)
            stat_str = 'CUBLAS_STATUS_SUCCESS'
        case (CUBLAS_STATUS_NOT_INITIALIZED)
            stat_str = 'CUBLAS_STATUS_NOT_INITIALIZED'
        case (CUBLAS_STATUS_ALLOC_FAILED)
            stat_str = 'CUBLAS_STATUS_ALLOC_FAILED'
        case (CUBLAS_STATUS_INVALID_VALUE)
            stat_str = 'CUBLAS_STATUS_INVALID_VALUE'
        case (CUBLAS_STATUS_ARCH_MISMATCH)
            stat_str = 'CUBLAS_STATUS_ARCH_MISMATCH'
        case (CUBLAS_STATUS_MAPPING_ERROR)
            stat_str = 'CUBLAS_STATUS_MAPPING_ERROR'
        case (CUBLAS_STATUS_EXECUTION_FAILED)
            stat_str = 'CUBLAS_STATUS_EXECUTION_FAILED'
        case (CUBLAS_STATUS_INTERNAL_ERROR)
            stat_str = 'CUBLAS_STATUS_INTERNAL_ERROR'
        case (CUBLAS_STATUS_NOT_SUPPORTED)
            stat_str = 'CUBLAS_STATUS_NOT_SUPPORTED'
        case (CUBLAS_STATUS_LICENSE_ERROR)
            stat_str = 'CUBLAS_STATUS_LICENSE_ERROR'
        case default
            stat_str = 'UNKNOWN_CUBLAS_ERROR'
    end select

    write(msg, '(A, I0, A, A, A)') 'cuBLAS error: ', stat, ' [', trim(stat_str), '] (' // trim(name) // ')'
    error stop msg
end subroutine
#endif
#:enddef
#:endmute
