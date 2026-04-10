module mfi_blas_extensions
#if defined(MFI_EXTENSIONS)
    use iso_fortran_env
    use f77_blas
#if defined(MFI_CUBLAS)
    use iso_c_binding
    use mfi_blas_cublas
#endif
    implicit none

! BLAS level 1 - Utils / Extensions
!> Generic modern interface for IAMAX.
!> Supports s, d, c, z.
!> See also:
!> [[f77_iamax:isamax]], [[f77_iamax:idamax]], [[f77_iamax:icamax]], [[f77_iamax:izamax]].
interface mfi_iamax
    module procedure :: mfi_isamax
    module procedure :: mfi_idamax
    module procedure :: mfi_icamax
    module procedure :: mfi_izamax
end interface
!> Generic modern interface for IAMIN.
!> Supports s, d, c, z.
!> See also:
!> [[f77_iamin:isamin]], [[f77_iamin:idamin]], [[f77_iamin:icamin]], [[f77_iamin:izamin]].
interface mfi_iamin
    module procedure :: mfi_isamin
    module procedure :: mfi_idamin
    module procedure :: mfi_icamin
    module procedure :: mfi_izamin
end interface

! Global variables for execution mode control - always defined
integer, save :: MFI_USE_CUBLAS = 0
integer, save :: MFI_USE_CUBLAS_PREV_STATE = 0

contains

! BLAS level 1 - Utils / Extensions
!> Modern interface for [[f77_iamax:f77_iamax]].
!> See also: [[mfi_iamax]], [[f77_iamax]].
pure function mfi_isamax(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_isamax
    real(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_isamax = f77_iamax(n,x,local_incx)
end function
!> Modern interface for [[f77_iamax:f77_iamax]].
!> See also: [[mfi_iamax]], [[f77_iamax]].
pure function mfi_idamax(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_idamax
    real(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_idamax = f77_iamax(n,x,local_incx)
end function
!> Modern interface for [[f77_iamax:f77_iamax]].
!> See also: [[mfi_iamax]], [[f77_iamax]].
pure function mfi_icamax(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_icamax
    complex(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_icamax = f77_iamax(n,x,local_incx)
end function
!> Modern interface for [[f77_iamax:f77_iamax]].
!> See also: [[mfi_iamax]], [[f77_iamax]].
pure function mfi_izamax(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_izamax
    complex(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_izamax = f77_iamax(n,x,local_incx)
end function
!> Modern interface for [[f77_iamin:f77_iamin]].
!> See also: [[mfi_iamin]], [[f77_iamin]].
pure function mfi_isamin(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_isamin
    real(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_isamin = f77_iamin(n,x,local_incx)
end function
!> Modern interface for [[f77_iamin:f77_iamin]].
!> See also: [[mfi_iamin]], [[f77_iamin]].
pure function mfi_idamin(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_idamin
    real(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_idamin = f77_iamin(n,x,local_incx)
end function
!> Modern interface for [[f77_iamin:f77_iamin]].
!> See also: [[mfi_iamin]], [[f77_iamin]].
pure function mfi_icamin(x, incx)
    integer, parameter :: wp = REAL32
    integer :: mfi_icamin
    complex(REAL32), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_icamin = f77_iamin(n,x,local_incx)
end function
!> Modern interface for [[f77_iamin:f77_iamin]].
!> See also: [[mfi_iamin]], [[f77_iamin]].
pure function mfi_izamin(x, incx)
    integer, parameter :: wp = REAL64
    integer :: mfi_izamin
    complex(REAL64), intent(in) :: x(:)
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer :: n
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    n = size(x)
    mfi_izamin = f77_iamin(n,x,local_incx)
end function

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
#if defined(MFI_CUBLAS)
    ! Check CUBLAS setting
    call get_environment_variable("MFI_USE_CUBLAS",env_mfi_use_cublas,env_len_cublas)
    if (env_len_cublas > 0) then
        read(env_mfi_use_cublas(1:env_len_cublas),*,iostat=info) MFI_USE_CUBLAS
    end if
    ! Initialize cuBLAS handle if GPU mode is enabled
    if (MFI_USE_CUBLAS == 1) then
        call mfi_cublas_handle_ensure()
    end if
#else
    ! Issue warning when used without proper compilation
    print *, 'WARNING: mfi_execution_init() called but compiled without CUBLAS support'
#endif
end subroutine

!> Ensure cuBLAS handle is created (called internally)
#if defined(MFI_CUBLAS)
subroutine mfi_cublas_handle_ensure()
    integer(c_int) :: stat
    if (.not. c_associated(mfi_cublas_handle)) then
        call cublasCreate(mfi_cublas_handle, stat)
        if (stat /= 0) error stop 'cublasCreate_v2 failed - check CUDA driver version'
        call cublasSetPointerMode(mfi_cublas_handle, CUBLAS_POINTER_MODE_HOST, stat)
        if (stat /= 0) error stop 'cublasSetPointerMode_v2 failed'
    end if
end subroutine
#endif

!> Sets execution mode to use GPU with CUBLAS (only effective when compiled with support)
subroutine mfi_force_gpu()
    MFI_USE_CUBLAS_PREV_STATE = MFI_USE_CUBLAS
#if defined(MFI_CUBLAS)
    MFI_USE_CUBLAS = 1
    call mfi_cublas_handle_ensure()
#else
    print *, 'WARNING: mfi_force_gpu() called but compiled without CUBLAS support'
    ! Don't actually enable GPU mode
#endif
end subroutine

!> Finalize cuBLAS resources
#if defined(MFI_CUBLAS)
subroutine mfi_cublas_finalize()
    integer(c_int) :: stat
    if (c_associated(mfi_cublas_handle)) then
        call cublasDestroy(mfi_cublas_handle, stat)
        mfi_cublas_handle = c_null_ptr
    end if
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
#if defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        mode_str = "GPU-CUBLAS"
    else
        mode_str = "CPU-DEFAULT"
    end if
#else
    mode_str = "CPU-ONLY"
#endif
end function
#endif
end module

