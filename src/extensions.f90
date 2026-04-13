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
public :: mfi_cublas_set_threads, mfi_cublas_handle_count

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
    type(c_ptr) :: h
    integer(c_int) :: stat
    call mfi_cublas_handle_get_c(h, stat)
    handle = h
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

!> Set number of cuBLAS handles for multi-threaded GPU use (no-op when no cuBLAS)
#if defined(MFI_CUBLAS)
subroutine mfi_cublas_set_threads(n)
    use mfi_blas_cublas, only: mfi_cublas_set_threads_c
    integer, intent(in) :: n
    call mfi_cublas_set_threads_c(n)
end subroutine
#else
subroutine mfi_cublas_set_threads(n)
    integer, intent(in) :: n
end subroutine
#endif

!> Return current cuBLAS handle count (0 when no cuBLAS)
#if defined(MFI_CUBLAS)
function mfi_cublas_handle_count() result(n)
    use mfi_blas_cublas, only: mfi_cublas_handle_count_c
    integer :: n
    n = mfi_cublas_handle_count_c()
end function
#else
function mfi_cublas_handle_count() result(n)
    integer :: n
    n = 0
end function
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
#endif
end module

