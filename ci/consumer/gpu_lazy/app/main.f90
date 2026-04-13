program main
    use mfi_blas, only: mfi_gemm, mfi_cublas_is_active
    use iso_fortran_env
    implicit none
    real(real64) :: A(4,4), B(4,4), C(4,4)

    print *, "UNIQUE_STRING_FROM_MAIN_12345"
    call random_number(A)
    call random_number(B)
    C = 0.0_real64

    ! This call will trigger lazy init because MFI_USE_CUBLAS=1 in CI
    call mfi_gemm(A, B, C)
    
    if (.not. mfi_cublas_is_active()) then
        print *, "FAIL: GPU mode not active after lazy init!"
        error stop 1
    else
        print *, "PASS: GPU mode active (auto-detected via env var)"
    end if

    print *, "PASS: mfi_gemm ran successfully via cuBLAS lazy init."
end program main
