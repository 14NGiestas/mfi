
program test_trtrs
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_strtrs 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trtrs against strtrs", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dtrtrs 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trtrs against dtrtrs", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_ctrtrs 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trtrs against ctrtrs", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_ztrtrs 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_trtrs against ztrtrs", t2-t1
end block
contains

subroutine test_strtrs
    use f77_lapack, only: strtrs, f77_trtrs
    use mfi_lapack, only: mfi_trtrs, mfi_strtrs

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3, NRHS = 2
    real(REAL32) :: A(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: info, info_rf, info_mfi

    ! Create a triangular matrix A
    A = reshape([4.0_wp, 0.0_wp, 0.0_wp, &  ! column 1
                 2.0_wp, 3.0_wp, 0.0_wp, &  ! column 2
                 1.0_wp, 2.0_wp, 5.0_wp], [N,N]) ! column 3

    ! Create B matrix such that A*X = [1,2; 2,1; 1,0] is the solution
    ! So B should be A * [1,2; 2,1; 1,0] = [4*1+0*2+0*1, 4*2+0*1+0*0; 2*1+3*2+0*1, 2*2+3*1+0*0; 1*1+2*2+5*1, 1*2+2*1+5*0]
    ! = [4,8; 8,7; 10,4]
    B = reshape([4.0_wp, 8.0_wp, 10.0_wp, &  ! column 1
                 8.0_wp, 7.0_wp, 4.0_wp], [N,NRHS]) ! column 2

    ! Test f77 interface for trtrs
    B_in = B
    call strtrs('U', 'N', 'N', N, NRHS, A, N, B_in, N, info)
    B_rf = B_in
    info_rf = info
    
    call assert(info == 0, "f77_strtrs failed", info)

    ! Test mfi interface (short form)
    B_in = B
    call mfi_strtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < 1e-10), &
                "different results for mfi_strtrs")

    ! Test mfi interface (full form)
    B_in = B
    call mfi_trtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < 1e-10), &
                "different results for mfi_trtrs")

end subroutine
subroutine test_dtrtrs
    use f77_lapack, only: dtrtrs, f77_trtrs
    use mfi_lapack, only: mfi_trtrs, mfi_dtrtrs

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3, NRHS = 2
    real(REAL64) :: A(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: info, info_rf, info_mfi

    ! Create a triangular matrix A
    A = reshape([4.0_wp, 0.0_wp, 0.0_wp, &  ! column 1
                 2.0_wp, 3.0_wp, 0.0_wp, &  ! column 2
                 1.0_wp, 2.0_wp, 5.0_wp], [N,N]) ! column 3

    ! Create B matrix such that A*X = [1,2; 2,1; 1,0] is the solution
    ! So B should be A * [1,2; 2,1; 1,0] = [4*1+0*2+0*1, 4*2+0*1+0*0; 2*1+3*2+0*1, 2*2+3*1+0*0; 1*1+2*2+5*1, 1*2+2*1+5*0]
    ! = [4,8; 8,7; 10,4]
    B = reshape([4.0_wp, 8.0_wp, 10.0_wp, &  ! column 1
                 8.0_wp, 7.0_wp, 4.0_wp], [N,NRHS]) ! column 2

    ! Test f77 interface for trtrs
    B_in = B
    call dtrtrs('U', 'N', 'N', N, NRHS, A, N, B_in, N, info)
    B_rf = B_in
    info_rf = info
    
    call assert(info == 0, "f77_dtrtrs failed", info)

    ! Test mfi interface (short form)
    B_in = B
    call mfi_dtrtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < 1e-10), &
                "different results for mfi_dtrtrs")

    ! Test mfi interface (full form)
    B_in = B
    call mfi_trtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < 1e-10), &
                "different results for mfi_trtrs")

end subroutine
subroutine test_ctrtrs
    use f77_lapack, only: ctrtrs, f77_trtrs
    use mfi_lapack, only: mfi_trtrs, mfi_ctrtrs

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3, NRHS = 2
    complex(REAL32) :: A(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: info, info_rf, info_mfi

    ! Create a triangular matrix A
    A = reshape([4.0_wp, 0.0_wp, 0.0_wp, &  ! column 1
                 2.0_wp, 3.0_wp, 0.0_wp, &  ! column 2
                 1.0_wp, 2.0_wp, 5.0_wp], [N,N]) ! column 3

    ! Create B matrix such that A*X = [1,2; 2,1; 1,0] is the solution
    ! So B should be A * [1,2; 2,1; 1,0] = [4*1+0*2+0*1, 4*2+0*1+0*0; 2*1+3*2+0*1, 2*2+3*1+0*0; 1*1+2*2+5*1, 1*2+2*1+5*0]
    ! = [4,8; 8,7; 10,4]
    B = reshape([4.0_wp, 8.0_wp, 10.0_wp, &  ! column 1
                 8.0_wp, 7.0_wp, 4.0_wp], [N,NRHS]) ! column 2

    ! Test f77 interface for trtrs
    B_in = B
    call ctrtrs('U', 'N', 'N', N, NRHS, A, N, B_in, N, info)
    B_rf = B_in
    info_rf = info
    
    call assert(info == 0, "f77_ctrtrs failed", info)

    ! Test mfi interface (short form)
    B_in = B
    call mfi_ctrtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < 1e-10), &
                "different results for mfi_ctrtrs")

    ! Test mfi interface (full form)
    B_in = B
    call mfi_trtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < 1e-10), &
                "different results for mfi_trtrs")

end subroutine
subroutine test_ztrtrs
    use f77_lapack, only: ztrtrs, f77_trtrs
    use mfi_lapack, only: mfi_trtrs, mfi_ztrtrs

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3, NRHS = 2
    complex(REAL64) :: A(N,N), B(N,NRHS), B_in(N,NRHS), B_rf(N,NRHS)
    integer :: info, info_rf, info_mfi

    ! Create a triangular matrix A
    A = reshape([4.0_wp, 0.0_wp, 0.0_wp, &  ! column 1
                 2.0_wp, 3.0_wp, 0.0_wp, &  ! column 2
                 1.0_wp, 2.0_wp, 5.0_wp], [N,N]) ! column 3

    ! Create B matrix such that A*X = [1,2; 2,1; 1,0] is the solution
    ! So B should be A * [1,2; 2,1; 1,0] = [4*1+0*2+0*1, 4*2+0*1+0*0; 2*1+3*2+0*1, 2*2+3*1+0*0; 1*1+2*2+5*1, 1*2+2*1+5*0]
    ! = [4,8; 8,7; 10,4]
    B = reshape([4.0_wp, 8.0_wp, 10.0_wp, &  ! column 1
                 8.0_wp, 7.0_wp, 4.0_wp], [N,NRHS]) ! column 2

    ! Test f77 interface for trtrs
    B_in = B
    call ztrtrs('U', 'N', 'N', N, NRHS, A, N, B_in, N, info)
    B_rf = B_in
    info_rf = info
    
    call assert(info == 0, "f77_ztrtrs failed", info)

    ! Test mfi interface (short form)
    B_in = B
    call mfi_ztrtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < 1e-10), &
                "different results for mfi_ztrtrs")

    ! Test mfi interface (full form)
    B_in = B
    call mfi_trtrs(A, B_in, info=info_mfi)
    call assert(info_mfi == info_rf .and. all(abs(B_in - B_rf) < 1e-10), &
                "different results for mfi_trtrs")

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

    subroutine report_test_result(test_name, success)
        character(*), intent(in) :: test_name
        logical, intent(in) :: success

        if (success) then
            write(*, '(A, ": ", A)') trim(test_name), 'PASSED'
        else
            write(*, '(A, ": ", A)') trim(test_name), 'FAILED'
        end if
    end subroutine

end program
