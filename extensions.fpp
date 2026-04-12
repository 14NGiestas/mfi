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
#:enddef

#:def mfi_extensions_implement()
! BLAS level 1 - Utils / Extensions
$:mfi_implement('i?amax', DEFAULT_TYPES, iamin_iamax)
$:mfi_implement('i?amin', DEFAULT_TYPES, iamin_iamax)

!> Returns .true. if GPU execution is active (pure, delegates to C)
pure logical function mfi_cublas_is_active() result(active)
#if defined(MFI_CUBLAS)
    use mfi_blas_cublas, only: mfi_cublas_lazy_init_c
    call mfi_cublas_lazy_init_c()
    active = mfi_cublas_is_active_c() /= 0
#else
    active = .false.
#endif
end function

!> Get the cuBLAS handle for the current thread (pure, delegates to C)
pure function mfi_cublas_handle_get() result(handle)
#if defined(MFI_CUBLAS)
    use mfi_blas_cublas, only: mfi_cublas_handle_get_c
    type(c_ptr) :: handle
    handle = mfi_cublas_handle_get_c()
    if (.not. c_associated(handle)) &
        error stop 'mfi: thread ID out of range. '// &
                   'Set OMP_NUM_THREADS before running with MFI_USE_CUBLAS=1.'
#else
    type(c_ptr) :: handle
    handle = c_null_ptr
#endif
end function

!> Switch to GPU mode (no-op when compiled without cuBLAS)
#if defined(MFI_CUBLAS)
subroutine mfi_force_gpu()
    use mfi_blas_cublas, only: mfi_cublas_force_gpu_c
    call mfi_cublas_force_gpu_c()
end subroutine
#else
subroutine mfi_force_gpu()
end subroutine
#endif

!> Switch to CPU mode (no-op when compiled without cuBLAS)
#if defined(MFI_CUBLAS)
subroutine mfi_force_cpu()
    use mfi_blas_cublas, only: mfi_cublas_force_cpu_c
    call mfi_cublas_force_cpu_c()
end subroutine
#else
subroutine mfi_force_cpu()
end subroutine
#endif

!> Finalize all cuBLAS handles (no-op when compiled without cuBLAS)
#if defined(MFI_CUBLAS)
subroutine mfi_cublas_finalize()
    use mfi_blas_cublas, only: mfi_cublas_finalize_all
    integer(c_int) :: stat
    call mfi_cublas_finalize_all(stat)
end subroutine
#else
subroutine mfi_cublas_finalize()
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
        case (CUBLAS_STATUS_SUCCESS);          stat_str = 'CUBLAS_STATUS_SUCCESS'
        case (CUBLAS_STATUS_NOT_INITIALIZED);  stat_str = 'CUBLAS_STATUS_NOT_INITIALIZED'
        case (CUBLAS_STATUS_ALLOC_FAILED);     stat_str = 'CUBLAS_STATUS_ALLOC_FAILED'
        case (CUBLAS_STATUS_INVALID_VALUE);    stat_str = 'CUBLAS_STATUS_INVALID_VALUE'
        case (CUBLAS_STATUS_ARCH_MISMATCH);    stat_str = 'CUBLAS_STATUS_ARCH_MISMATCH'
        case (CUBLAS_STATUS_MAPPING_ERROR);    stat_str = 'CUBLAS_STATUS_MAPPING_ERROR'
        case (CUBLAS_STATUS_EXECUTION_FAILED); stat_str = 'CUBLAS_STATUS_EXECUTION_FAILED'
        case (CUBLAS_STATUS_INTERNAL_ERROR);   stat_str = 'CUBLAS_STATUS_INTERNAL_ERROR'
        case (CUBLAS_STATUS_NOT_SUPPORTED);    stat_str = 'CUBLAS_STATUS_NOT_SUPPORTED'
        case (CUBLAS_STATUS_LICENSE_ERROR);    stat_str = 'CUBLAS_STATUS_LICENSE_ERROR'
        case default;                          stat_str = 'UNKNOWN_CUBLAS_ERROR'
    end select

    write(msg, '(A, I0, A, A, A)') 'cuBLAS error: ', stat, ' [', trim(stat_str), '] (' // trim(name) // ')'
    error stop msg
end subroutine
#endif
#:enddef
#:endmute
