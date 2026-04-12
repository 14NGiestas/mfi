program main
    use mfi_blas, only: mfi_gemm, mfi_cublas_is_active
    use iso_fortran_env
    implicit none
    real(real64) :: A(4,4), B(4,4), C(4,4)

    ! Check if env var activation works (no explicit force_gpu call)
    if (mfi_cublas_is_active()) then
        print *, "GPU mode active (auto-detected via env var)"
    else
        print *, "CPU mode active"
    end if

    call random_number(A)
    call random_number(B)
    C = 0.0_real64
    call mfi_gemm(A, B, C)
    print *, "Consumer build successful, mfi_gemm ran."
end program main
