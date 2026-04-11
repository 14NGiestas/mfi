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

!> Check environment and initialize cuBLAS on first GPU call (lazy init)
#if defined(MFI_CUBLAS)
subroutine mfi_cublas_lazy_init()
    integer :: info
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
        call mfi_cublas_preallocate_handles()
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
    use mfi_blas_cublas, only: mfi_cublas_handles
    type(c_ptr) :: handle
    call mfi_cublas_get_thread_handle(mfi_cublas_handles, mfi_cublas_handle_count, handle)
#else
    type(c_ptr) :: handle
    handle = c_null_ptr
#endif
end function

!> Pre-allocate all cuBLAS handles (called during lazy init, serial)
#if defined(MFI_CUBLAS)
subroutine mfi_cublas_preallocate_handles()
    integer(c_int) :: stat
    integer :: i
    if (mfi_cublas_handle_count <= 0) return
    do i = 1, mfi_cublas_handle_count
        if (.not. c_associated(mfi_cublas_handles(i))) then
            call cublasCreate(mfi_cublas_handles(i), stat)
            if (stat /= 0) error stop 'cublasCreate_v2 failed - check CUDA driver version'
            call cublasSetPointerMode(mfi_cublas_handles(i), CUBLAS_POINTER_MODE_HOST, stat)
            if (stat /= 0) error stop 'cublasSetPointerMode_v2 failed'
        end if
    end do
end subroutine
#endif

!> Switch to GPU mode (no-op when compiled without cuBLAS)
#if defined(MFI_CUBLAS)
subroutine mfi_force_gpu()
    MFI_USE_CUBLAS = 1
    mfi_cublas_global_initialized = .false.
    mfi_cublas_env_checked = .true.
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
    integer(c_int) :: stat
    integer :: i
    do i = 1, mfi_cublas_handle_count
        if (c_associated(mfi_cublas_handles(i))) then
            call cublasDestroy(mfi_cublas_handles(i), stat)
            mfi_cublas_handles(i) = c_null_ptr
        end if
    end do
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
#endif
end module

