module mfi_blas_trmm
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

!> Generic modern interface for TRMM.
!> Supports s, d, c, z.
!> See also:
!> [[f77_trmm:strmm]], [[f77_trmm:dtrmm]], [[f77_trmm:ctrmm]], [[f77_trmm:ztrmm]].
interface mfi_trmm
    module procedure :: mfi_strmm
    module procedure :: mfi_dtrmm
    module procedure :: mfi_ctrmm
    module procedure :: mfi_ztrmm
end interface

contains

!> Modern interface for [[f77_trmm:f77_trmm]].
!> See also: [[mfi_trmm]], [[f77_trmm]].
pure subroutine mfi_strmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in), target :: a(:,:)
    real(REAL32), intent(inout), target :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(REAL32), intent(in), optional :: alpha
    real(REAL32) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_b
        integer(c_int) :: cu_side, cu_uplo, cu_transa, cu_diag, cublas_stat
        real(REAL32), target :: alpha_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_b, &
                              int(size(b) * storage_size(b)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating b'
        if (local_side == 'L' .or. local_side == 'l') then
            cu_side = CUBLAS_SIDE_LEFT
        else
            cu_side = CUBLAS_SIDE_RIGHT
        end if
        if (local_uplo == 'U' .or. local_uplo == 'u') then
            cu_uplo = CUBLAS_TRSM_FILL_UPPER
        else
            cu_uplo = CUBLAS_TRSM_FILL_LOWER
        end if
        if (local_transa == 'N' .or. local_transa == 'n') then
            cu_transa = CUBLAS_OP_N
        else if (local_transa == 'T' .or. local_transa == 't') then
            cu_transa = CUBLAS_OP_T
        else
            cu_transa = CUBLAS_OP_C
        end if
        if (local_diag == 'U' .or. local_diag == 'u') then
            cu_diag = CUBLAS_DIAG_UNIT
        else
            cu_diag = CUBLAS_DIAG_NON_UNIT
        end if
        alpha_target = local_alpha
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_b, c_loc(b), &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasStrmm(mfi_cublas_handle, cu_side, cu_uplo, cu_transa, cu_diag, &
                 int(m,c_int), int(n,c_int), c_loc(alpha_target), &
                 device_a, int(lda,c_int), device_b, int(ldb,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'f77_trmm')
        call cudaMemcpy(c_loc(b), device_b, &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_b, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating b'
        end block
        return
    end if
#endif
    call f77_trmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
!> Modern interface for [[f77_trmm:f77_trmm]].
!> See also: [[mfi_trmm]], [[f77_trmm]].
pure subroutine mfi_dtrmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in), target :: a(:,:)
    real(REAL64), intent(inout), target :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    real(REAL64), intent(in), optional :: alpha
    real(REAL64) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_b
        integer(c_int) :: cu_side, cu_uplo, cu_transa, cu_diag, cublas_stat
        real(REAL64), target :: alpha_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_b, &
                              int(size(b) * storage_size(b)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating b'
        if (local_side == 'L' .or. local_side == 'l') then
            cu_side = CUBLAS_SIDE_LEFT
        else
            cu_side = CUBLAS_SIDE_RIGHT
        end if
        if (local_uplo == 'U' .or. local_uplo == 'u') then
            cu_uplo = CUBLAS_TRSM_FILL_UPPER
        else
            cu_uplo = CUBLAS_TRSM_FILL_LOWER
        end if
        if (local_transa == 'N' .or. local_transa == 'n') then
            cu_transa = CUBLAS_OP_N
        else if (local_transa == 'T' .or. local_transa == 't') then
            cu_transa = CUBLAS_OP_T
        else
            cu_transa = CUBLAS_OP_C
        end if
        if (local_diag == 'U' .or. local_diag == 'u') then
            cu_diag = CUBLAS_DIAG_UNIT
        else
            cu_diag = CUBLAS_DIAG_NON_UNIT
        end if
        alpha_target = local_alpha
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_b, c_loc(b), &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasDtrmm(mfi_cublas_handle, cu_side, cu_uplo, cu_transa, cu_diag, &
                 int(m,c_int), int(n,c_int), c_loc(alpha_target), &
                 device_a, int(lda,c_int), device_b, int(ldb,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'f77_trmm')
        call cudaMemcpy(c_loc(b), device_b, &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_b, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating b'
        end block
        return
    end if
#endif
    call f77_trmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
!> Modern interface for [[f77_trmm:f77_trmm]].
!> See also: [[mfi_trmm]], [[f77_trmm]].
pure subroutine mfi_ctrmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in), target :: a(:,:)
    complex(REAL32), intent(inout), target :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(REAL32), intent(in), optional :: alpha
    complex(REAL32) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_b
        integer(c_int) :: cu_side, cu_uplo, cu_transa, cu_diag, cublas_stat
        complex(REAL32), target :: alpha_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_b, &
                              int(size(b) * storage_size(b)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating b'
        if (local_side == 'L' .or. local_side == 'l') then
            cu_side = CUBLAS_SIDE_LEFT
        else
            cu_side = CUBLAS_SIDE_RIGHT
        end if
        if (local_uplo == 'U' .or. local_uplo == 'u') then
            cu_uplo = CUBLAS_TRSM_FILL_UPPER
        else
            cu_uplo = CUBLAS_TRSM_FILL_LOWER
        end if
        if (local_transa == 'N' .or. local_transa == 'n') then
            cu_transa = CUBLAS_OP_N
        else if (local_transa == 'T' .or. local_transa == 't') then
            cu_transa = CUBLAS_OP_T
        else
            cu_transa = CUBLAS_OP_C
        end if
        if (local_diag == 'U' .or. local_diag == 'u') then
            cu_diag = CUBLAS_DIAG_UNIT
        else
            cu_diag = CUBLAS_DIAG_NON_UNIT
        end if
        alpha_target = local_alpha
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_b, c_loc(b), &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasCtrmm(mfi_cublas_handle, cu_side, cu_uplo, cu_transa, cu_diag, &
                 int(m,c_int), int(n,c_int), c_loc(alpha_target), &
                 device_a, int(lda,c_int), device_b, int(ldb,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'f77_trmm')
        call cudaMemcpy(c_loc(b), device_b, &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_b, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating b'
        end block
        return
    end if
#endif
    call f77_trmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
!> Modern interface for [[f77_trmm:f77_trmm]].
!> See also: [[mfi_trmm]], [[f77_trmm]].
pure subroutine mfi_ztrmm(a, b, side, uplo, transa, diag, alpha)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in), target :: a(:,:)
    complex(REAL64), intent(inout), target :: b(:,:)
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: transa
    character :: local_transa
    character, intent(in), optional :: diag
    character :: local_diag
    complex(REAL64), intent(in), optional :: alpha
    complex(REAL64) :: local_alpha
    integer :: m, n, lda, ldb
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(transa)) then
        local_transa = transa
    else
        local_transa = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    if (present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1.0_wp
    end if
    m = size(b,1)
    n = size(b,2)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
#if defined(MFI_EXTENSIONS) && defined(MFI_CUBLAS)
    if (MFI_USE_CUBLAS == 1) then
        block
    integer(c_int) :: cuda_allocation_status
        type(c_ptr) :: device_a
        type(c_ptr) :: device_b
        integer(c_int) :: cu_side, cu_uplo, cu_transa, cu_diag, cublas_stat
        complex(REAL64), target :: alpha_target
        call cuda_malloc(device_a, &
                              int(size(a) * storage_size(a)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating a'
        call cuda_malloc(device_b, &
                              int(size(b) * storage_size(b)/8, c_size_t), &
                              cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaMalloc failed allocating b'
        if (local_side == 'L' .or. local_side == 'l') then
            cu_side = CUBLAS_SIDE_LEFT
        else
            cu_side = CUBLAS_SIDE_RIGHT
        end if
        if (local_uplo == 'U' .or. local_uplo == 'u') then
            cu_uplo = CUBLAS_TRSM_FILL_UPPER
        else
            cu_uplo = CUBLAS_TRSM_FILL_LOWER
        end if
        if (local_transa == 'N' .or. local_transa == 'n') then
            cu_transa = CUBLAS_OP_N
        else if (local_transa == 'T' .or. local_transa == 't') then
            cu_transa = CUBLAS_OP_T
        else
            cu_transa = CUBLAS_OP_C
        end if
        if (local_diag == 'U' .or. local_diag == 'u') then
            cu_diag = CUBLAS_DIAG_UNIT
        else
            cu_diag = CUBLAS_DIAG_NON_UNIT
        end if
        alpha_target = local_alpha
        call cudaMemcpy(device_a, c_loc(a), &
                        int(size(a) * storage_size(a)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        call cudaMemcpy(device_b, c_loc(b), &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyHostToDevice)
        cublas_stat = cublasZtrmm(mfi_cublas_handle, cu_side, cu_uplo, cu_transa, cu_diag, &
                 int(m,c_int), int(n,c_int), c_loc(alpha_target), &
                 device_a, int(lda,c_int), device_b, int(ldb,c_int))
        if (cublas_stat /= 0) call mfi_cublas_error(cublas_stat, 'f77_trmm')
        call cudaMemcpy(c_loc(b), device_b, &
                        int(size(b) * storage_size(b)/8, c_size_t), &
                        cudaMemcpyDeviceToHost)
        call cuda_free(device_a, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating a'
        call cuda_free(device_b, cuda_allocation_status)
        if (cuda_allocation_status /= 0) error stop 'cudaFree failed deallocating b'
        end block
        return
    end if
#endif
    call f77_trmm(local_side,local_uplo,local_transa,local_diag,m,n,local_alpha,a,lda,b,ldb)
end subroutine
end module

