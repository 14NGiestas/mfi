program test_cublas_gemv
    use iso_fortran_env
    use iso_c_binding
    use mfi_blas, only: mfi_execution_init, mfi_sgemv, mfi_get_execution_mode
    implicit none
    real(real32) :: M(4,4), X(4), Y(4), Y_ref(4)
    real(real32) :: alpha, beta
    character :: trans
    integer :: i, j
    
    call mfi_execution_init()
    print *, "Execution mode:", trim(mfi_get_execution_mode())
    
    call random_number(M)
    call random_number(X)
    Y = 0.0
    Y_ref = 0.0
    alpha = 2.0
    beta = 0.5
    
    trans = 'N'
    do i = 1, 4
        do j = 1, 4
            Y_ref(i) = Y_ref(i) + alpha * M(i,j) * X(j)
        end do
        Y_ref(i) = Y_ref(i) + beta * Y(i)
    end do
    
    Y = 0.0
    call mfi_sgemv(M, X, Y, alpha=alpha, beta=beta, trans=trans)
    
    print *, "CPU reference Y:", Y_ref
    print *, "MFI result Y:   ", Y
    print *, "Difference:      ", abs(Y - Y_ref)
    
    if (all(abs(Y - Y_ref) < 1.0e-5)) then
        print *, "SUCCESS"
    else
        print *, "FAILURE"
    end if
end program
