program trsm_diag2
    use iso_fortran_env
    use iso_c_binding
    use mfi_blas
    implicit none

    integer, parameter :: wp = REAL32
    real(wp), target :: A(3,3), B(3,1), B_ref(3,1), B_gpu(3,1)
    real(wp) :: alpha
    type(c_ptr) :: d_a, d_b
    integer(c_int) :: stat
    real(wp), target :: alpha_t
    integer :: s, u, d
    integer(c_int), parameter :: sides(2) = [0_c_int, 1_c_int]
    integer(c_int), parameter :: uplos(2) = [0_c_int, 1_c_int]
    integer(c_int), parameter :: diags(2) = [0_c_int, 1_c_int]
    character, parameter :: sname(2) = ['L', 'R']
    character, parameter :: uname(2) = ['U', 'L']
    character, parameter :: dname(2) = ['N', 'U']

    ! Simple upper triangular system
    A = 0.0_wp
    A(1,1) = 2.0_wp; A(1,2) = 1.0_wp
    A(2,2) = 3.0_wp; A(2,3) = 1.0_wp
    A(3,3) = 4.0_wp
    B = reshape([1.0_wp, 2.0_wp, 3.0_wp], [3,1])
    alpha = 1.0_wp

    call mfi_execution_init()
    print *, "Mode:", trim(mfi_get_execution_mode())

    do s = 1, 2
    do u = 1, 2
    do d = 1, 2
        ! CPU reference
        B_ref = B
        call f77_trsm(sname(s), uname(u), 'N', dname(d), 3, 1, alpha, A, 3, B_ref, 3)

        ! cuBLAS
        stat = cuda_malloc(d_a, int(size(A)*storage_size(A)/8, c_size_t))
        stat = cuda_malloc(d_b, int(size(B)*storage_size(B)/8, c_size_t))
        B_gpu = B
        call cudaMemcpy(d_a, c_loc(A), int(size(A)*storage_size(A)/8, c_size_t), cudaMemcpyHostToDevice)
        call cudaMemcpy(d_b, c_loc(B_gpu), int(size(B_gpu)*storage_size(B_gpu)/8, c_size_t), cudaMemcpyHostToDevice)
        alpha_t = alpha
        stat = cublasStrsm(mfi_cublas_handle, sides(s), uplos(u), 0_c_int, diags(d), &
                 3_c_int, 1_c_int, c_loc(alpha_t), d_a, 3_c_int, d_b, 3_c_int)
        call cudaMemcpy(c_loc(B_gpu), d_b, int(size(B_gpu)*storage_size(B_gpu)/8, c_size_t), cudaMemcpyDeviceToHost)
        stat = cuda_free(d_a)
        stat = cuda_free(d_b)

        if (all(abs(B_gpu - B_ref) < 1e-5_wp)) then
            print '("  side=",A1," uplo=",A1," diag=",A1," -> MATCH (stat=",I0,")")', &
                sname(s), uname(u), dname(d), stat
        else
            print '("  side=",A1," uplo=",A1," diag=",A1," -> DIFF max=",F12.6," (stat=",I0,")")', &
                sname(s), uname(u), dname(d), maxval(abs(B_gpu - B_ref)), stat
            print '("    CPU: ",3F12.6)', B_ref(:,1)
            print '("    GPU: ",3F12.6)', B_gpu(:,1)
        end if
    end do
    end do
    end do

end program
