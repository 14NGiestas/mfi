program main
    use mfi_blas, only: mfi_gemm, mfi_cublas_is_active
    use iso_fortran_env
    implicit none
    real(real64) :: A(4,4), B(4,4), C(4,4)

    if (mfi_cublas_is_active()) then
        print *, "FAIL: GPU mode active, but cpu-basic should be CPU-only!"
        error stop 1
    else
        print *, "PASS: CPU mode active"
    end if

    call random_number(A)
    call random_number(B)
    C = 0.0_real64
    call mfi_gemm(A, B, C)
    print *, "PASS: mfi_gemm ran successfully on CPU."
end program main