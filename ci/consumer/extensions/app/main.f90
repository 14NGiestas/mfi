program main
    use mfi_blas, only: mfi_iamax, mfi_cublas_is_active, mfi_cublas_finalize
    use iso_fortran_env
    implicit none
    real(real64) :: x(5)
    integer :: imax

    ! 1. Call extension before GPU is active (runs on CPU)
    x = [1.0_real64, 3.0_real64, 5.0_real64, 2.0_real64, 4.0_real64]
    imax = mfi_iamax(x)
    if (imax /= 3) then
        print *, "FAIL: mfi_iamax returned", imax, "expected 3"
        error stop 1
    end if
    print *, "PASS: mfi_iamax ran correctly on CPU."
    
    ! 2. Test cublas finalize (should be safe to call even if not active)
    call mfi_cublas_finalize()
    print *, "PASS: mfi_cublas_finalize ran successfully."

end program main