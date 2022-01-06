
program test_mfi_lapack
    use iso_fortran_env
    use mfi_lapack
    use f77_lapack
    implicit none
    integer, parameter :: N = 2000
    real(REAL64) :: S(N,N)
    integer :: i, j, info

    call test_gesvd
    call test_potrf
    call test_potri

contains

    subroutine symmetric
        call random_number(S)
        do j=1,N
            do i=j+1,N
                S(i,j) = S(j,i)
            end do
        end do
    end subroutine

    subroutine positive_definite
        call symmetric
        S = matmul(S,S)
    end subroutine

    subroutine test_gesvd
        real(REAL64) :: A(3,5)
        real(REAL64) :: S(3), ES(3)
        ES = [real(REAL64) :: 14.0828, 10.1124, 3.92600]
        A(1,:) = [real(REAL64) :: -2,  1, -5, 11, 2]
        A(2,:) = [real(REAL64) ::  2, -2,  0,  1, 4]
        A(3,:) = [real(REAL64) ::  6, -8,  0,  6, 0]
block
real :: t1, t2
call cpu_time(t1)
 call mfi_gesvd(A,S) 
call cpu_time(t2)
print '(A,G0)', "mfi_gesvd: ", t2-t1
end block
        call assert(all(abs(S-ES) < 1e-4))
    end subroutine

    subroutine test_potrf
        call positive_definite
block
real :: t1, t2
call cpu_time(t1)
 call f77_potrf('U',N,S,N,info) 
call cpu_time(t2)
print '(A,G0)', "f77_potrf: ", t2-t1
end block
        call assert(info == 0)

        call positive_definite
block
real :: t1, t2
call cpu_time(t1)
 call mfi_potrf(S,info)         
call cpu_time(t2)
print '(A,G0)', "mfi_potrf: ", t2-t1
end block
        call assert(info == 0)
    end subroutine

    subroutine test_potri
        call positive_definite
        call f77_potrf('U',N,S,N,info)
block
real :: t1, t2
call cpu_time(t1)
 call f77_potri('U',N,S,N,info) 
call cpu_time(t2)
print '(A,G0)', "f77_potri: ", t2-t1
end block
        call assert(info == 0)

        call positive_definite
        call mfi_potrf(S,info)
block
real :: t1, t2
call cpu_time(t1)
 call mfi_potri(S,info)         
call cpu_time(t2)
print '(A,G0)', "mfi_potri: ", t2-t1
end block
        call assert(info == 0)
    end subroutine

    pure subroutine assert(test)
        logical, intent(in) :: test
        if (.not. test) then
            error stop 'assertion failed'
        end if
    end subroutine

    logical pure elemental function is_almost_equal(x, y)
        real(REAL64), intent(in) :: x, y
        is_almost_equal = abs(x-y) < 10**6*epsilon(x)
    end function

end program
