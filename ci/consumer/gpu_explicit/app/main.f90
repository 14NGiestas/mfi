program main
    use mfi_blas, only: mfi_gemm, mfi_cublas_is_active, mfi_force_gpu, mfi_force_cpu
    use iso_fortran_env
    implicit none
    real(real64) :: A(4,4), B(4,4), C(4,4)

    call random_number(A)
    call random_number(B)
    C = 0.0_real64

    ! 1. Default state (no env var should be set here)
    if (mfi_cublas_is_active()) then
        print *, "FAIL: GPU mode active by default!"
        error stop 1
    end if

    ! 2. Explicit GPU activation
    call mfi_force_gpu()
    if (.not. mfi_cublas_is_active()) then
        print *, "FAIL: GPU mode not active after mfi_force_gpu()!"
        error stop 1
    end if
    
    call mfi_gemm(A, B, C)
    print *, "PASS: mfi_gemm ran via explicit cuBLAS."

    ! 3. Explicit CPU fallback
    call mfi_force_cpu()
    if (mfi_cublas_is_active()) then
        print *, "FAIL: GPU mode still active after mfi_force_cpu()!"
        error stop 1
    end if
    
    call mfi_gemm(A, B, C)
    print *, "PASS: mfi_gemm ran via explicit CPU fallback."
end program main