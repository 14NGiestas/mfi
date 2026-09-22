program main
    use mfi_blas, only: mfi_gemm, mfi_hip_is_active, mfi_force_gpu, mfi_force_cpu
    use iso_fortran_env
    implicit none
    real(real64) :: A(4,4), B(4,4), C(4,4), C_ref(4,4)
    integer :: i

    call random_number(A)
    call random_number(B)
    C = 0.0_real64

    ! 1. Default state (no env var should be set here)
    if (mfi_hip_is_active()) then
        print *, "FAIL: GPU mode active by default!"
        error stop 1
    end if

    call mfi_gemm(A, B, C)
    C_ref = C

    ! 2. Explicit GPU activation
    call mfi_force_gpu()
    if (.not. mfi_hip_is_active()) then
        print *, "FAIL: GPU mode not active after mfi_force_gpu()!"
        error stop 1
    end if

    C = 0.0_real64
    call mfi_gemm(A, B, C)
    do i = 1, 4
        if (maxval(abs(C(i,:) - C_ref(i,:))) > sqrt(epsilon(1.0_real64))) then
            print *, "FAIL: rocBLAS gemm result differs from the CPU result!"
            error stop 1
        end if
    end do
    print *, "PASS: mfi_gemm ran via native rocBLAS."

    ! 3. Explicit CPU fallback
    call mfi_force_cpu()
    if (mfi_hip_is_active()) then
        print *, "FAIL: GPU mode still active after mfi_force_cpu()!"
        error stop 1
    end if

    C = 0.0_real64
    call mfi_gemm(A, B, C)
    do i = 1, 4
        if (maxval(abs(C(i,:) - C_ref(i,:))) > sqrt(epsilon(1.0_real64))) then
            print *, "FAIL: CPU fallback result differs from the rocBLAS result!"
            error stop 1
        end if
    end do
    print *, "PASS: mfi_gemm ran via the CPU fallback."
end program main
