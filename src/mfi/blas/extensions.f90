module mfi_blas_extensions
    use iso_fortran_env
    use f77_blas
    use iso_c_binding
    use mfi_blas_cublas
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
    ! Check CUBLAS setting
    call get_environment_variable("MFI_USE_CUBLAS",env_mfi_use_cublas,env_len_cublas)
    if (env_len_cublas > 0) then
        read(env_mfi_use_cublas(1:env_len_cublas),*,iostat=info) MFI_USE_CUBLAS
    end if
    ! Initialize cuBLAS handle if GPU mode is enabled
    if (MFI_USE_CUBLAS == 1) then
        call mfi_cublas_handle_ensure()
    end if
end subroutine

!> Ensure cuBLAS handle is created (called internally)
subroutine mfi_cublas_handle_ensure()
    integer(c_int) :: stat
    if (.not. c_associated(mfi_cublas_handle)) then
        call cublasCreate(mfi_cublas_handle, stat)
        if (stat /= 0) error stop 'cublasCreate_v2 failed - check CUDA driver version'
        call cublasSetPointerMode(mfi_cublas_handle, CUBLAS_POINTER_MODE_HOST, stat)
        if (stat /= 0) error stop 'cublasSetPointerMode_v2 failed'
    end if
end subroutine

!> Sets execution mode to use GPU with CUBLAS (only effective when compiled with support)
subroutine mfi_force_gpu()
    MFI_USE_CUBLAS_PREV_STATE = MFI_USE_CUBLAS
    MFI_USE_CUBLAS = 1
    call mfi_cublas_handle_ensure()
end subroutine

!> Finalize cuBLAS resources
subroutine mfi_cublas_finalize()
    integer(c_int) :: stat
    if (c_associated(mfi_cublas_handle)) then
        call cublasDestroy(mfi_cublas_handle, stat)
        mfi_cublas_handle = c_null_ptr
    end if
end subroutine

!> Report cuBLAS error (called from pure wrappers)
pure subroutine mfi_cublas_error(stat, name)
    integer(c_int), value, intent(in) :: stat
    character(*), intent(in) :: name
    character(len=120) :: msg
    write(msg, '(A, I0, A)') 'cuBLAS error: ', stat, ' (' // trim(name) // ')'
    error stop msg
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
    if (MFI_USE_CUBLAS == 1) then
        mode_str = "GPU-CUBLAS"
    else
        mode_str = "CPU-DEFAULT"
    end if
end function
end module

