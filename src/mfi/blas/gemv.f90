module mfi_blas_gemv
    use iso_fortran_env
    use f77_blas
#if defined(MFI_CUBLAS)
    use iso_c_binding
    use mfi_blas_cublas
#endif
#if defined(MFI_EXTENSIONS)
    use mfi_blas_extensions
#endif
    implicit none

!> Generic modern interface for GEMV.
!> Supports s, d, c, z.
!> See also:
!> [[f77_gemv:sgemv]], [[f77_gemv:dgemv]], [[f77_gemv:cgemv]], [[f77_gemv:zgemv]].
interface mfi_gemv
    module procedure :: mfi_sgemv
    module procedure :: mfi_dgemv
    module procedure :: mfi_cgemv
    module procedure :: mfi_zgemv
end interface

contains

!> Modern interface for [[f77_gemv:f77_gemv]].
!> See also: [[mfi_gemv]], [[f77_gemv]].
pure subroutine mfi_sgemv(a, x, y, trans, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in), target :: a(:,:), x(:)
    real(REAL32), intent(inout), target :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    if (present(incy)) then
        local_incy = incy
    else
        local_incy = 1
    end if
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_x
        type(c_ptr) :: device_y
        integer(c_int) :: op, cublas_stat
        real(REAL32), target :: alpha_target, beta_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_x, &
                              int(size(x) * storage_size(x)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating x'
        call cuda_malloc(device_y, &
                              int(size(y) * storage_size(y)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating y'
        if (local_trans == 'N' .or. local_trans == 'n') then
            op = CUBLAS_OP_N
        else if (local_trans == 'T' .or. local_trans == 't') then
            op = CUBLAS_OP_T
        else
            op = CUBLAS_OP_C
        end if
        alpha_target = local_alpha
        beta_target = local_beta
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_x, c_loc(x), &
                        int(size(x) * storage_size(x)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_y, c_loc(y), &
                        int(size(y) * storage_size(y)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasSgemv(mfi_cublas_handle, op, &
                 int(m,c_int), int(n,c_int), &
                 c_loc(alpha_target), device_a, int(lda,c_int), &
                 device_x, int(local_incx,c_int), &
                 c_loc(beta_target), device_y, int(local_incy,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'Sgemv')
        call cudaMemcpy(c_loc(y), device_y, &
                        int(size(y) * storage_size(y)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_x, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating x'
        call cuda_free(device_y, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating y'
        end block
        return
    end if
#endif
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gemv:f77_gemv]].
!> See also: [[mfi_gemv]], [[f77_gemv]].
pure subroutine mfi_dgemv(a, x, y, trans, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in), target :: a(:,:), x(:)
    real(REAL64), intent(inout), target :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    if (present(incy)) then
        local_incy = incy
    else
        local_incy = 1
    end if
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_x
        type(c_ptr) :: device_y
        integer(c_int) :: op, cublas_stat
        real(REAL64), target :: alpha_target, beta_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_x, &
                              int(size(x) * storage_size(x)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating x'
        call cuda_malloc(device_y, &
                              int(size(y) * storage_size(y)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating y'
        if (local_trans == 'N' .or. local_trans == 'n') then
            op = CUBLAS_OP_N
        else if (local_trans == 'T' .or. local_trans == 't') then
            op = CUBLAS_OP_T
        else
            op = CUBLAS_OP_C
        end if
        alpha_target = local_alpha
        beta_target = local_beta
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_x, c_loc(x), &
                        int(size(x) * storage_size(x)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_y, c_loc(y), &
                        int(size(y) * storage_size(y)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasDgemv(mfi_cublas_handle, op, &
                 int(m,c_int), int(n,c_int), &
                 c_loc(alpha_target), device_a, int(lda,c_int), &
                 device_x, int(local_incx,c_int), &
                 c_loc(beta_target), device_y, int(local_incy,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'Dgemv')
        call cudaMemcpy(c_loc(y), device_y, &
                        int(size(y) * storage_size(y)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_x, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating x'
        call cuda_free(device_y, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating y'
        end block
        return
    end if
#endif
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gemv:f77_gemv]].
!> See also: [[mfi_gemv]], [[f77_gemv]].
pure subroutine mfi_cgemv(a, x, y, trans, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in), target :: a(:,:), x(:)
    complex(REAL32), intent(inout), target :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    if (present(incy)) then
        local_incy = incy
    else
        local_incy = 1
    end if
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_x
        type(c_ptr) :: device_y
        integer(c_int) :: op, cublas_stat
        complex(REAL32), target :: alpha_target, beta_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_x, &
                              int(size(x) * storage_size(x)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating x'
        call cuda_malloc(device_y, &
                              int(size(y) * storage_size(y)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating y'
        if (local_trans == 'N' .or. local_trans == 'n') then
            op = CUBLAS_OP_N
        else if (local_trans == 'T' .or. local_trans == 't') then
            op = CUBLAS_OP_T
        else
            op = CUBLAS_OP_C
        end if
        alpha_target = local_alpha
        beta_target = local_beta
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_x, c_loc(x), &
                        int(size(x) * storage_size(x)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_y, c_loc(y), &
                        int(size(y) * storage_size(y)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasCgemv(mfi_cublas_handle, op, &
                 int(m,c_int), int(n,c_int), &
                 c_loc(alpha_target), device_a, int(lda,c_int), &
                 device_x, int(local_incx,c_int), &
                 c_loc(beta_target), device_y, int(local_incy,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'Cgemv')
        call cudaMemcpy(c_loc(y), device_y, &
                        int(size(y) * storage_size(y)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_x, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating x'
        call cuda_free(device_y, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating y'
        end block
        return
    end if
#endif
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
!> Modern interface for [[f77_gemv:f77_gemv]].
!> See also: [[mfi_gemv]], [[f77_gemv]].
pure subroutine mfi_zgemv(a, x, y, trans, alpha, beta, incx, incy)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in), target :: a(:,:), x(:)
    complex(REAL64), intent(inout), target :: y(:)
    character, intent(in), optional :: trans
    character :: local_trans
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
    integer, intent(in), optional :: incx
    integer :: local_incx
    integer, intent(in), optional :: incy
    integer :: local_incy
    integer :: m, n, lda
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    if (present(beta)) then
        local_beta = beta
    else
        local_beta = 0.0_wp
    end if
    if (present(incx)) then
        local_incx = incx
    else
        local_incx = 1
    end if
    if (present(incy)) then
        local_incy = incy
    else
        local_incy = 1
    end if
    m = size(a,1)
    n = size(a,2)
    lda = max(1,m)
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_x
        type(c_ptr) :: device_y
        integer(c_int) :: op, cublas_stat
        complex(REAL64), target :: alpha_target, beta_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_x, &
                              int(size(x) * storage_size(x)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating x'
        call cuda_malloc(device_y, &
                              int(size(y) * storage_size(y)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating y'
        if (local_trans == 'N' .or. local_trans == 'n') then
            op = CUBLAS_OP_N
        else if (local_trans == 'T' .or. local_trans == 't') then
            op = CUBLAS_OP_T
        else
            op = CUBLAS_OP_C
        end if
        alpha_target = local_alpha
        beta_target = local_beta
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_x, c_loc(x), &
                        int(size(x) * storage_size(x)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_y, c_loc(y), &
                        int(size(y) * storage_size(y)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasZgemv(mfi_cublas_handle, op, &
                 int(m,c_int), int(n,c_int), &
                 c_loc(alpha_target), device_a, int(lda,c_int), &
                 device_x, int(local_incx,c_int), &
                 c_loc(beta_target), device_y, int(local_incy,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'Zgemv')
        call cudaMemcpy(c_loc(y), device_y, &
                        int(size(y) * storage_size(y)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_x, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating x'
        call cuda_free(device_y, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating y'
        end block
        return
    end if
#endif
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,local_incx,local_beta,y,local_incy)
end subroutine
end module

