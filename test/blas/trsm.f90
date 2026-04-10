

program test_trsm
use iso_fortran_env
implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_strsm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trsm against strsm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dtrsm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trsm against dtrsm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_ctrsm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trsm against ctrsm", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_ztrsm 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trsm against ztrsm", t2-t1
end block
contains
subroutine test_strsm
    use f77_blas, only: strsm, f77_trsm
    use mfi_blas, only: mfi_trsm
    use mfi_blas

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 10, N = 10
    real(REAL32) :: A(M,M), B(M,N), B_in(M,N), B_rf(M,N)
    real(REAL32) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag

    character, parameter :: sides(*) = ['L', 'R']
    character, parameter :: uplos(*) = ['U', 'L']
    character, parameter :: transas(*) = ['N', 'T', 'C']
    character, parameter :: diags(*) = ['N', 'U']

    call mfi_execution_init
    print*, "Testing mfi_trsm with strsm, mode: ", mfi_get_execution_mode()

    ! Initialize A as a triangular matrix to ensure it is non-singular
    A = 0.0_wp
    do i_side = 1, M
        A(i_side, i_side) = 1.0_wp + real(i_side, wp)
        do i_uplo = 1, i_side - 1
            A(i_uplo, i_side) = real(i_uplo + i_side, wp)
        end do
    end do
    
    call random_number(B)
    call random_number(alpha)

    do i_side = 1, size(sides)
    do i_uplo = 1, size(uplos)
    do i_transa = 1, size(transas)
    do i_diag = 1, size(diags)
        side = sides(i_side)
        uplo = uplos(i_uplo)
        transa = transas(i_transa)
        diag = diags(i_diag)

        ! Reference run
        B_rf = B
        call strsm(side, uplo, transa, diag, M, N, alpha, A, M, B_rf, M)

        ! MFI run
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        block
        real(wp) :: max_diff, scale
        max_diff = maxval(abs(B_in - B_rf))
        scale = max(1.0_wp, maxval(abs(B_rf)))
        if (.not. max_diff < sqrt(epsilon(1.0_wp)) * scale) then
            print *, "FAIL: side=", side, " uplo=", uplo, " transa=", transa, " diag=", diag
            print *, "max diff =", max_diff
            print *, "B_rf(1,1) =", B_rf(1,1)
            print *, "B_in(1,1) =", B_in(1,1)
            error stop "mfi diff"
        end if
        end block

    end do
    end do
    end do
    end do
end subroutine
subroutine test_dtrsm
    use f77_blas, only: dtrsm, f77_trsm
    use mfi_blas, only: mfi_trsm
    use mfi_blas

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 10, N = 10
    real(REAL64) :: A(M,M), B(M,N), B_in(M,N), B_rf(M,N)
    real(REAL64) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag

    character, parameter :: sides(*) = ['L', 'R']
    character, parameter :: uplos(*) = ['U', 'L']
    character, parameter :: transas(*) = ['N', 'T', 'C']
    character, parameter :: diags(*) = ['N', 'U']

    call mfi_execution_init
    print*, "Testing mfi_trsm with dtrsm, mode: ", mfi_get_execution_mode()

    ! Initialize A as a triangular matrix to ensure it is non-singular
    A = 0.0_wp
    do i_side = 1, M
        A(i_side, i_side) = 1.0_wp + real(i_side, wp)
        do i_uplo = 1, i_side - 1
            A(i_uplo, i_side) = real(i_uplo + i_side, wp)
        end do
    end do
    
    call random_number(B)
    call random_number(alpha)

    do i_side = 1, size(sides)
    do i_uplo = 1, size(uplos)
    do i_transa = 1, size(transas)
    do i_diag = 1, size(diags)
        side = sides(i_side)
        uplo = uplos(i_uplo)
        transa = transas(i_transa)
        diag = diags(i_diag)

        ! Reference run
        B_rf = B
        call dtrsm(side, uplo, transa, diag, M, N, alpha, A, M, B_rf, M)

        ! MFI run
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        block
        real(wp) :: max_diff, scale
        max_diff = maxval(abs(B_in - B_rf))
        scale = max(1.0_wp, maxval(abs(B_rf)))
        if (.not. max_diff < sqrt(epsilon(1.0_wp)) * scale) then
            print *, "FAIL: side=", side, " uplo=", uplo, " transa=", transa, " diag=", diag
            print *, "max diff =", max_diff
            print *, "B_rf(1,1) =", B_rf(1,1)
            print *, "B_in(1,1) =", B_in(1,1)
            error stop "mfi diff"
        end if
        end block

    end do
    end do
    end do
    end do
end subroutine
subroutine test_ctrsm
    use f77_blas, only: ctrsm, f77_trsm
    use mfi_blas, only: mfi_trsm
    use mfi_blas

    integer, parameter :: wp = REAL32
    integer, parameter :: M = 10, N = 10
    complex(REAL32) :: A(M,M), B(M,N), B_in(M,N), B_rf(M,N)
    complex(REAL32) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag

    character, parameter :: sides(*) = ['L', 'R']
    character, parameter :: uplos(*) = ['U', 'L']
    character, parameter :: transas(*) = ['N', 'T', 'C']
    character, parameter :: diags(*) = ['N', 'U']

    call mfi_execution_init
    print*, "Testing mfi_trsm with ctrsm, mode: ", mfi_get_execution_mode()

    ! Initialize A as a triangular matrix to ensure it is non-singular
    A = 0.0_wp
    do i_side = 1, M
        A(i_side, i_side) = 1.0_wp + real(i_side, wp)
        do i_uplo = 1, i_side - 1
            A(i_uplo, i_side) = real(i_uplo + i_side, wp)
        end do
    end do
    
block
    real(REAL32) :: re(M,N)
    real(REAL32) :: im(M,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im, kind=REAL32)
end block
block
    real(REAL32) :: re
    real(REAL32) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im, kind=REAL32)
end block

    do i_side = 1, size(sides)
    do i_uplo = 1, size(uplos)
    do i_transa = 1, size(transas)
    do i_diag = 1, size(diags)
        side = sides(i_side)
        uplo = uplos(i_uplo)
        transa = transas(i_transa)
        diag = diags(i_diag)

        ! Reference run
        B_rf = B
        call ctrsm(side, uplo, transa, diag, M, N, alpha, A, M, B_rf, M)

        ! MFI run
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        block
        real(wp) :: max_diff, scale
        max_diff = maxval(abs(B_in - B_rf))
        scale = max(1.0_wp, maxval(abs(B_rf)))
        if (.not. max_diff < sqrt(epsilon(1.0_wp)) * scale) then
            print *, "FAIL: side=", side, " uplo=", uplo, " transa=", transa, " diag=", diag
            print *, "max diff =", max_diff
            print *, "B_rf(1,1) =", B_rf(1,1)
            print *, "B_in(1,1) =", B_in(1,1)
            error stop "mfi diff"
        end if
        end block

    end do
    end do
    end do
    end do
end subroutine
subroutine test_ztrsm
    use f77_blas, only: ztrsm, f77_trsm
    use mfi_blas, only: mfi_trsm
    use mfi_blas

    integer, parameter :: wp = REAL64
    integer, parameter :: M = 10, N = 10
    complex(REAL64) :: A(M,M), B(M,N), B_in(M,N), B_rf(M,N)
    complex(REAL64) :: alpha
    character :: side, uplo, transa, diag
    integer :: i_side, i_uplo, i_transa, i_diag

    character, parameter :: sides(*) = ['L', 'R']
    character, parameter :: uplos(*) = ['U', 'L']
    character, parameter :: transas(*) = ['N', 'T', 'C']
    character, parameter :: diags(*) = ['N', 'U']

    call mfi_execution_init
    print*, "Testing mfi_trsm with ztrsm, mode: ", mfi_get_execution_mode()

    ! Initialize A as a triangular matrix to ensure it is non-singular
    A = 0.0_wp
    do i_side = 1, M
        A(i_side, i_side) = 1.0_wp + real(i_side, wp)
        do i_uplo = 1, i_side - 1
            A(i_uplo, i_side) = real(i_uplo + i_side, wp)
        end do
    end do
    
block
    real(REAL64) :: re(M,N)
    real(REAL64) :: im(M,N)
    call random_number(im)
    call random_number(re)
    B = cmplx(re,im, kind=REAL64)
end block
block
    real(REAL64) :: re
    real(REAL64) :: im
    call random_number(im)
    call random_number(re)
    alpha = cmplx(re,im, kind=REAL64)
end block

    do i_side = 1, size(sides)
    do i_uplo = 1, size(uplos)
    do i_transa = 1, size(transas)
    do i_diag = 1, size(diags)
        side = sides(i_side)
        uplo = uplos(i_uplo)
        transa = transas(i_transa)
        diag = diags(i_diag)

        ! Reference run
        B_rf = B
        call ztrsm(side, uplo, transa, diag, M, N, alpha, A, M, B_rf, M)

        ! MFI run
        B_in = B
        call mfi_trsm(A, B_in, side=side, uplo=uplo, transa=transa, diag=diag, alpha=alpha)
        block
        real(wp) :: max_diff, scale
        max_diff = maxval(abs(B_in - B_rf))
        scale = max(1.0_wp, maxval(abs(B_rf)))
        if (.not. max_diff < sqrt(epsilon(1.0_wp)) * scale) then
            print *, "FAIL: side=", side, " uplo=", uplo, " transa=", transa, " diag=", diag
            print *, "max diff =", max_diff
            print *, "B_rf(1,1) =", B_rf(1,1)
            print *, "B_in(1,1) =", B_in(1,1)
            error stop "mfi diff"
        end if
        end block

    end do
    end do
    end do
    end do
end subroutine

subroutine assert(test, msg, info)
    logical, intent(in) :: test
    character(*), intent(in) :: msg
    integer, intent(in), optional :: info
    character(1024) :: buffer

    if (.not. test) then
        if (present(info)) then
            write(buffer, *) 'Error ', info, ': ', msg
        else
            write(buffer, *) 'Error: ', msg
        end if
        error stop trim(buffer)
    end if
end subroutine

end program
