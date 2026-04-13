program main
    use mfi_lapack, only: mfi_getrf
    use iso_fortran_env
    implicit none
    real(real64) :: A(4,4)
    integer :: ipiv(4)
    integer :: info

    call random_number(A)
    
    ! Perform LU factorization
    call mfi_getrf(A, ipiv, info)
    
    if (info /= 0) then
        print *, "FAIL: mfi_getrf returned info =", info
        error stop 1
    end if
    
    print *, "PASS: mfi_getrf ran successfully."
end program main