program trsm_diag
    use iso_fortran_env
    use iso_c_binding
    use mfi_blas
    implicit none

    integer, parameter :: wp = REAL32
    real(wp) :: A(3,3), B(3,1), B_ref(3,1), B_gpu(3,1)
    real(wp) :: alpha
    type(c_ptr) :: d_a, d_b
    integer(c_int) :: stat
    real(wp), target :: alpha_t
    integer :: i

    ! Simple known upper triangular system: A * x = alpha * b
    ! A = [2 1 0; 0 3 1; 0 0 4], b = [1; 2; 3], alpha = 1
    ! Expected: x = A\b = [0.0833; 0.4167; 0.75]
    A = 0.0_wp
    A(1,1) = 2.0_wp; A(1,2) = 1.0_wp
    A(2,2) = 3.0_wp; A(2,3) = 1.0_wp
    A(3,3) = 4.0_wp
    B = reshape([1.0_wp, 2.0_wp, 3.0_wp], [3,1])
    alpha = 1.0_wp

    call mfi_execution_init()
    print *, "Mode:", trim(mfi_get_execution_mode())

    print *, "A (upper triangular):"
    do i = 1, 3
        print '(3F8.3)', A(i,:)
    end do
    print *, "B =", B(:,1)

    ! Reference: CPU BLAS  (side=L, uplo=U, trans=N, diag=N)
    B_ref = B
    call f77_trsm('L', 'U', 'N', 'N', 3, 1, alpha, A, 3, B_ref, 3)
    print *, "CPU TRSM result:", B_ref(:,1)

    if (MFI_USE_CUBLAS == 1) then
        ! cuBLAS direct call with enum 0/1
        stat = cuda_malloc(d_a, int(size(A)*storage_size(A)/8, c_size_t))
        if (stat /= 0) print *, "cuda_malloc A failed:", stat
        stat = cuda_malloc(d_b, int(size(B)*storage_size(B)/8, c_size_t))
        if (stat /= 0) print *, "cuda_malloc B failed:", stat

        B_gpu = B
        call cudaMemcpy(d_a, c_loc(A), int(size(A)*storage_size(A)/8, c_size_t), cudaMemcpyHostToDevice)
        call cudaMemcpy(d_b, c_loc(B_gpu), int(size(B_gpu)*storage_size(B_gpu)/8, c_size_t), cudaMemcpyHostToDevice)

        alpha_t = alpha
        ! side=0(L), uplo=0(U), trans=0(N), diag=0(NON_UNIT)
        stat = cublasStrsm(mfi_cublas_handle, 0_c_int, 0_c_int, 0_c_int, 0_c_int, &
                 3_c_int, 1_c_int, c_loc(alpha_t), d_a, 3_c_int, d_b, 3_c_int)
        print *, "cublasStrsm(0,0,0,0) stat =", stat

        call cudaMemcpy(c_loc(B_gpu), d_b, int(size(B_gpu)*storage_size(B_gpu)/8, c_size_t), cudaMemcpyDeviceToHost)
        print *, "cuBLAS (0,0,0,0) result:", B_gpu(:,1)
        if (all(abs(B_gpu - B_ref) < 1e-5_wp)) then
            print *, "  => MATCH"
        else
            print *, "  => DIFF max =", maxval(abs(B_gpu - B_ref))
        end if

        ! Try side=1 (RIGHT) with same other params
        B_gpu = B
        call cudaMemcpy(d_b, c_loc(B_gpu), int(size(B_gpu)*storage_size(B_gpu)/8, c_size_t), cudaMemcpyHostToDevice)
        stat = cublasStrsm(mfi_cublas_handle, 1_c_int, 0_c_int, 0_c_int, 0_c_int, &
                 3_c_int, 1_c_int, c_loc(alpha_t), d_a, 3_c_int, d_b, 3_c_int)
        print *, "cublasStrsm(1,0,0,0) stat =", stat
        call cudaMemcpy(c_loc(B_gpu), d_b, int(size(B_gpu)*storage_size(B_gpu)/8, c_size_t), cudaMemcpyDeviceToHost)
        print *, "cuBLAS (1,0,0,0) result:", B_gpu(:,1)

        ! Try uplo=1 (LOWER)
        B_gpu = B
        call cudaMemcpy(d_b, c_loc(B_gpu), int(size(B_gpu)*storage_size(B_gpu)/8, c_size_t), cudaMemcpyHostToDevice)
        stat = cublasStrsm(mfi_cublas_handle, 0_c_int, 1_c_int, 0_c_int, 0_c_int, &
                 3_c_int, 1_c_int, c_loc(alpha_t), d_a, 3_c_int, d_b, 3_c_int)
        print *, "cublasStrsm(0,1,0,0) stat =", stat
        call cudaMemcpy(c_loc(B_gpu), d_b, int(size(B_gpu)*storage_size(B_gpu)/8, c_size_t), cudaMemcpyDeviceToHost)
        print *, "cuBLAS (0,1,0,0) result:", B_gpu(:,1)

        stat = cuda_free(d_a)
        stat = cuda_free(d_b)
    end if

end program
