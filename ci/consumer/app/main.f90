program main
    use mfi_blas, only: mfi_gemm
    use iso_fortran_env
    implicit none
    real(real64) :: A(4,4), B(4,4), C(4,4)
    call random_number(A)
    call random_number(B)
    C = 0.0_real64
    call mfi_gemm(A, B, C)
    print *, "Consumer build successful, mfi_gemm ran."
end program main
