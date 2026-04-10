module mfi_blas_gemm
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

!> Generic modern interface for GEMM.
!> Supports s, d, c, z.
!> See also:
!> [[f77_gemm:sgemm]], [[f77_gemm:dgemm]], [[f77_gemm:cgemm]], [[f77_gemm:zgemm]].
interface mfi_gemm
    module procedure :: mfi_sgemm
    module procedure :: mfi_dgemm
    module procedure :: mfi_cgemm
    module procedure :: mfi_zgemm
end interface

contains

!> Modern interface for [[f77_gemm:f77_gemm]].
!> See also: [[mfi_gemm]], [[f77_gemm]].
pure subroutine mfi_sgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in), target :: a(:,:), b(:,:)
    real(REAL32), intent(inout), target :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    real(REAL32), intent(in), optional :: beta
    real(REAL32) :: local_beta
    integer :: m, n, k, lda, ldb, ldc
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(transb)) then
        local_transb = transb
    else
        local_transb = 'N'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_b
        type(c_ptr) :: device_c
        integer(c_int) :: op_a, op_b, cublas_stat
        real(REAL32), target :: alpha_target, beta_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_b, &
                              int(size(b) * storage_size(b)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating b'
        call cuda_malloc(device_c, &
                              int(size(c) * storage_size(c)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating c'
        if (local_transa == 'N' .or. local_transa == 'n') then
            op_a = CUBLAS_OP_N
        else if (local_transa == 'T' .or. local_transa == 't') then
            op_a = CUBLAS_OP_T
        else
            op_a = CUBLAS_OP_C
        end if
        if (local_transb == 'N' .or. local_transb == 'n') then
            op_b = CUBLAS_OP_N
        else if (local_transb == 'T' .or. local_transb == 't') then
            op_b = CUBLAS_OP_T
        else
            op_b = CUBLAS_OP_C
        end if
        alpha_target = local_alpha
        beta_target = local_beta
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_b, c_loc(b), &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_c, c_loc(c), &
                        int(size(c) * storage_size(c)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasSgemm(mfi_cublas_handle, op_a, op_b, &
                 int(m,c_int), int(n,c_int), int(k,c_int), &
                 c_loc(alpha_target), device_a, int(lda,c_int), &
                 device_b, int(ldb,c_int), c_loc(beta_target), device_c, int(ldc,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'Sgemm')
        call cudaMemcpy(c_loc(c), device_c, &
                        int(size(c) * storage_size(c)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_b, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating b'
        call cuda_free(device_c, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating c'
        end block
        return
    end if
#endif
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
!> Modern interface for [[f77_gemm:f77_gemm]].
!> See also: [[mfi_gemm]], [[f77_gemm]].
pure subroutine mfi_dgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in), target :: a(:,:), b(:,:)
    real(REAL64), intent(inout), target :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    real(REAL64), intent(in), optional :: beta
    real(REAL64) :: local_beta
    integer :: m, n, k, lda, ldb, ldc
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(transb)) then
        local_transb = transb
    else
        local_transb = 'N'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_b
        type(c_ptr) :: device_c
        integer(c_int) :: op_a, op_b, cublas_stat
        real(REAL64), target :: alpha_target, beta_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_b, &
                              int(size(b) * storage_size(b)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating b'
        call cuda_malloc(device_c, &
                              int(size(c) * storage_size(c)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating c'
        if (local_transa == 'N' .or. local_transa == 'n') then
            op_a = CUBLAS_OP_N
        else if (local_transa == 'T' .or. local_transa == 't') then
            op_a = CUBLAS_OP_T
        else
            op_a = CUBLAS_OP_C
        end if
        if (local_transb == 'N' .or. local_transb == 'n') then
            op_b = CUBLAS_OP_N
        else if (local_transb == 'T' .or. local_transb == 't') then
            op_b = CUBLAS_OP_T
        else
            op_b = CUBLAS_OP_C
        end if
        alpha_target = local_alpha
        beta_target = local_beta
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_b, c_loc(b), &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_c, c_loc(c), &
                        int(size(c) * storage_size(c)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasDgemm(mfi_cublas_handle, op_a, op_b, &
                 int(m,c_int), int(n,c_int), int(k,c_int), &
                 c_loc(alpha_target), device_a, int(lda,c_int), &
                 device_b, int(ldb,c_int), c_loc(beta_target), device_c, int(ldc,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'Dgemm')
        call cudaMemcpy(c_loc(c), device_c, &
                        int(size(c) * storage_size(c)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_b, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating b'
        call cuda_free(device_c, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating c'
        end block
        return
    end if
#endif
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
!> Modern interface for [[f77_gemm:f77_gemm]].
!> See also: [[mfi_gemm]], [[f77_gemm]].
pure subroutine mfi_cgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in), target :: a(:,:), b(:,:)
    complex(REAL32), intent(inout), target :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    complex(REAL32), intent(in), optional :: beta
    complex(REAL32) :: local_beta
    integer :: m, n, k, lda, ldb, ldc
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(transb)) then
        local_transb = transb
    else
        local_transb = 'N'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_b
        type(c_ptr) :: device_c
        integer(c_int) :: op_a, op_b, cublas_stat
        complex(REAL32), target :: alpha_target, beta_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_b, &
                              int(size(b) * storage_size(b)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating b'
        call cuda_malloc(device_c, &
                              int(size(c) * storage_size(c)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating c'
        if (local_transa == 'N' .or. local_transa == 'n') then
            op_a = CUBLAS_OP_N
        else if (local_transa == 'T' .or. local_transa == 't') then
            op_a = CUBLAS_OP_T
        else
            op_a = CUBLAS_OP_C
        end if
        if (local_transb == 'N' .or. local_transb == 'n') then
            op_b = CUBLAS_OP_N
        else if (local_transb == 'T' .or. local_transb == 't') then
            op_b = CUBLAS_OP_T
        else
            op_b = CUBLAS_OP_C
        end if
        alpha_target = local_alpha
        beta_target = local_beta
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_b, c_loc(b), &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_c, c_loc(c), &
                        int(size(c) * storage_size(c)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasCgemm(mfi_cublas_handle, op_a, op_b, &
                 int(m,c_int), int(n,c_int), int(k,c_int), &
                 c_loc(alpha_target), device_a, int(lda,c_int), &
                 device_b, int(ldb,c_int), c_loc(beta_target), device_c, int(ldc,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'Cgemm')
        call cudaMemcpy(c_loc(c), device_c, &
                        int(size(c) * storage_size(c)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_b, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating b'
        call cuda_free(device_c, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating c'
        end block
        return
    end if
#endif
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
!> Modern interface for [[f77_gemm:f77_gemm]].
!> See also: [[mfi_gemm]], [[f77_gemm]].
pure subroutine mfi_zgemm(a, b, c, transa, transb, alpha, beta)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in), target :: a(:,:), b(:,:)
    complex(REAL64), intent(inout), target :: c(:,:)
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: transb
    character :: local_transb
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    complex(REAL64), intent(in), optional :: beta
    complex(REAL64) :: local_beta
    integer :: m, n, k, lda, ldb, ldc
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(transb)) then
        local_transb = transb
    else
        local_transb = 'N'
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
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    if (local_transa == 'N' .or. local_transa == 'n') then
        k = size(a,2)
    else
        k = size(a,1)
    end if
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_b
        type(c_ptr) :: device_c
        integer(c_int) :: op_a, op_b, cublas_stat
        complex(REAL64), target :: alpha_target, beta_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_b, &
                              int(size(b) * storage_size(b)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating b'
        call cuda_malloc(device_c, &
                              int(size(c) * storage_size(c)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating c'
        if (local_transa == 'N' .or. local_transa == 'n') then
            op_a = CUBLAS_OP_N
        else if (local_transa == 'T' .or. local_transa == 't') then
            op_a = CUBLAS_OP_T
        else
            op_a = CUBLAS_OP_C
        end if
        if (local_transb == 'N' .or. local_transb == 'n') then
            op_b = CUBLAS_OP_N
        else if (local_transb == 'T' .or. local_transb == 't') then
            op_b = CUBLAS_OP_T
        else
            op_b = CUBLAS_OP_C
        end if
        alpha_target = local_alpha
        beta_target = local_beta
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_b, c_loc(b), &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_c, c_loc(c), &
                        int(size(c) * storage_size(c)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasZgemm(mfi_cublas_handle, op_a, op_b, &
                 int(m,c_int), int(n,c_int), int(k,c_int), &
                 c_loc(alpha_target), device_a, int(lda,c_int), &
                 device_b, int(ldb,c_int), c_loc(beta_target), device_c, int(ldc,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'Zgemm')
        call cudaMemcpy(c_loc(c), device_c, &
                        int(size(c) * storage_size(c)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_b, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating b'
        call cuda_free(device_c, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating c'
        end block
        return
    end if
#endif
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
end subroutine
end module

