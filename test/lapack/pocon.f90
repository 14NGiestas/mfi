
program test_pocon
    use iso_fortran_env
    implicit none
block
real :: t1, t2
call cpu_time(t1)
 call test_spocon 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_pocon against spocon", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_dpocon 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_pocon against dpocon", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_cpocon 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_pocon against cpocon", t2-t1
end block
block
real :: t1, t2
call cpu_time(t1)
 call test_zpocon 
call cpu_time(t2)
print '(A," (",G0,"s)")', "testing mfi_pocon against zpocon", t2-t1
end block
contains

subroutine test_spocon
    use f77_lapack, only: spocon, f77_pocon
    use mfi_lapack, only: mfi_pocon, mfi_spocon

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    real(REAL32) :: A(N,N)
    real(REAL32) :: anorm, rcond, rcond_rf, rcond_mfi
    integer :: info, info_rf, info_mfi
    integer :: i, j
    real(REAL32) :: col_sum
    real(REAL32), allocatable :: work(:)
    integer, allocatable :: iwork(:)

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Calculate the 1-norm of A (maximum absolute column sum)
    anorm = 0.0_wp
    do j = 1, N
        col_sum = 0.0_wp
        do i = 1, N
            col_sum = col_sum + abs(A(i,j))
        end do
        if (col_sum > anorm) anorm = col_sum
    end do

    ! Test f77 interface for pocon
    allocate(work(3*N))  ! Real pocon functions use 3*N for work array
    allocate(iwork(N))
    call spocon('U', N, A, N, anorm, rcond, work, iwork, info)
    if (allocated(work)) deallocate(work)
    if (allocated(iwork)) deallocate(iwork)
    rcond_rf = rcond
    info_rf = info

    ! Only continue if the condition number computation was successful
    if (info == 0) then
        ! Test mfi interface (short form)
        call mfi_spocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-10, &
                    "different results for mfi_spocon")

        ! Test mfi interface (full form)
        call mfi_pocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-10, &
                    "different results for mfi_pocon")
    end if

end subroutine
subroutine test_dpocon
    use f77_lapack, only: dpocon, f77_pocon
    use mfi_lapack, only: mfi_pocon, mfi_dpocon

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    real(REAL64) :: A(N,N)
    real(REAL64) :: anorm, rcond, rcond_rf, rcond_mfi
    integer :: info, info_rf, info_mfi
    integer :: i, j
    real(REAL64) :: col_sum
    real(REAL64), allocatable :: work(:)
    integer, allocatable :: iwork(:)

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Calculate the 1-norm of A (maximum absolute column sum)
    anorm = 0.0_wp
    do j = 1, N
        col_sum = 0.0_wp
        do i = 1, N
            col_sum = col_sum + abs(A(i,j))
        end do
        if (col_sum > anorm) anorm = col_sum
    end do

    ! Test f77 interface for pocon
    allocate(work(3*N))  ! Real pocon functions use 3*N for work array
    allocate(iwork(N))
    call dpocon('U', N, A, N, anorm, rcond, work, iwork, info)
    if (allocated(work)) deallocate(work)
    if (allocated(iwork)) deallocate(iwork)
    rcond_rf = rcond
    info_rf = info

    ! Only continue if the condition number computation was successful
    if (info == 0) then
        ! Test mfi interface (short form)
        call mfi_dpocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-10, &
                    "different results for mfi_dpocon")

        ! Test mfi interface (full form)
        call mfi_pocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-10, &
                    "different results for mfi_pocon")
    end if

end subroutine
subroutine test_cpocon
    use f77_lapack, only: cpocon, f77_pocon
    use mfi_lapack, only: mfi_pocon, mfi_cpocon

    integer, parameter :: wp = REAL32
    integer, parameter :: N = 3
    complex(REAL32) :: A(N,N)
    real(REAL32) :: anorm, rcond, rcond_rf, rcond_mfi
    integer :: info, info_rf, info_mfi
    integer :: i, j
    real(REAL32) :: col_sum
    complex(REAL32), allocatable :: work(:)
    real(REAL32), allocatable :: rwork(:)

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Calculate the 1-norm of A (maximum absolute column sum)
    anorm = 0.0_wp
    do j = 1, N
        col_sum = 0.0_wp
        do i = 1, N
            col_sum = col_sum + abs(A(i,j))
        end do
        if (col_sum > anorm) anorm = col_sum
    end do

    ! Test f77 interface for pocon
    allocate(work(2*N))
    allocate(rwork(N))
    call cpocon('U', N, A, N, anorm, rcond, work, rwork, info)
    if (allocated(work)) deallocate(work)
    if (allocated(rwork)) deallocate(rwork)
    rcond_rf = rcond
    info_rf = info

    ! Only continue if the condition number computation was successful
    if (info == 0) then
        ! Test mfi interface (short form)
        call mfi_cpocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-10, &
                    "different results for mfi_cpocon")

        ! Test mfi interface (full form)
        call mfi_pocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-10, &
                    "different results for mfi_pocon")
    end if

end subroutine
subroutine test_zpocon
    use f77_lapack, only: zpocon, f77_pocon
    use mfi_lapack, only: mfi_pocon, mfi_zpocon

    integer, parameter :: wp = REAL64
    integer, parameter :: N = 3
    complex(REAL64) :: A(N,N)
    real(REAL64) :: anorm, rcond, rcond_rf, rcond_mfi
    integer :: info, info_rf, info_mfi
    integer :: i, j
    real(REAL64) :: col_sum
    complex(REAL64), allocatable :: work(:)
    real(REAL64), allocatable :: rwork(:)

    ! Create a positive definite matrix
    A = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
                 1.0_wp, 4.0_wp, 1.0_wp, &
                 1.0_wp, 1.0_wp, 4.0_wp], [N,N])

    ! Calculate the 1-norm of A (maximum absolute column sum)
    anorm = 0.0_wp
    do j = 1, N
        col_sum = 0.0_wp
        do i = 1, N
            col_sum = col_sum + abs(A(i,j))
        end do
        if (col_sum > anorm) anorm = col_sum
    end do

    ! Test f77 interface for pocon
    allocate(work(2*N))
    allocate(rwork(N))
    call zpocon('U', N, A, N, anorm, rcond, work, rwork, info)
    if (allocated(work)) deallocate(work)
    if (allocated(rwork)) deallocate(rwork)
    rcond_rf = rcond
    info_rf = info

    ! Only continue if the condition number computation was successful
    if (info == 0) then
        ! Test mfi interface (short form)
        call mfi_zpocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-10, &
                    "different results for mfi_zpocon")

        ! Test mfi interface (full form)
        call mfi_pocon(A, anorm, rcond_mfi, info=info_mfi)
        call assert(info_mfi == info_rf .and. abs(rcond_mfi - rcond_rf) < 1e-10, &
                    "different results for mfi_pocon")
    end if

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
