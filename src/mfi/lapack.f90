module mfi_lapack
use iso_fortran_env
use f77_lapack
use f77_lapack, only: mfi_lartg => f77_lartg
implicit none

interface mfi_geqrf
    module procedure mfi_sgeqrf
    module procedure mfi_dgeqrf
    module procedure mfi_cgeqrf
    module procedure mfi_zgeqrf
end interface
interface mfi_gerqf
    module procedure mfi_sgerqf
    module procedure mfi_dgerqf
    module procedure mfi_cgerqf
    module procedure mfi_zgerqf
end interface
interface mfi_getrf
    module procedure mfi_sgetrf
    module procedure mfi_dgetrf
    module procedure mfi_cgetrf
    module procedure mfi_zgetrf
end interface
interface mfi_getri
    module procedure mfi_sgetri
    module procedure mfi_dgetri
    module procedure mfi_cgetri
    module procedure mfi_zgetri
end interface
interface mfi_getrs
    module procedure mfi_sgetrs
    module procedure mfi_dgetrs
    module procedure mfi_cgetrs
    module procedure mfi_zgetrs
end interface
interface mfi_hetrf
    module procedure mfi_chetrf
    module procedure mfi_zhetrf
end interface
interface mfi_hegv
    module procedure mfi_chegv
    module procedure mfi_zhegv
end interface
interface mfi_heevd
    module procedure mfi_cheevd
    module procedure mfi_zheevd
end interface
interface mfi_gesvd
    module procedure mfi_sgesvd
    module procedure mfi_dgesvd
    module procedure mfi_cgesvd
    module procedure mfi_zgesvd
end interface
interface mfi_potrf
    module procedure mfi_spotrf
    module procedure mfi_dpotrf
    module procedure mfi_cpotrf
    module procedure mfi_zpotrf
end interface
interface mfi_potri
    module procedure mfi_spotri
    module procedure mfi_dpotri
    module procedure mfi_cpotri
    module procedure mfi_zpotri
end interface
interface mfi_potrs
    module procedure mfi_spotrs
    module procedure mfi_dpotrs
    module procedure mfi_cpotrs
    module procedure mfi_zpotrs
end interface

contains

pure subroutine mfi_sgeqrf(a, tau, info)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(wp), pointer :: local_tau(:), work(:)
    real(wp), target  :: s_work(1)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(tau)) then
        local_tau => tau
    else
        allocate(local_tau(min(m,n)), stat=allocation_status)
    end if
    ! Retrieve work array size
    lwork = -1
    call sgeqrf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call sgeqrf(m,n,a,lda,local_tau,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(tau)) then
        deallocate(local_tau, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sgeqrf', -local_info)
    end if
end subroutine
pure subroutine mfi_dgeqrf(a, tau, info)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(wp), pointer :: local_tau(:), work(:)
    real(wp), target  :: s_work(1)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(tau)) then
        local_tau => tau
    else
        allocate(local_tau(min(m,n)), stat=allocation_status)
    end if
    ! Retrieve work array size
    lwork = -1
    call dgeqrf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call dgeqrf(m,n,a,lda,local_tau,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(tau)) then
        deallocate(local_tau, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dgeqrf', -local_info)
    end if
end subroutine
pure subroutine mfi_cgeqrf(a, tau, info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(wp), pointer :: local_tau(:), work(:)
    complex(wp), target  :: s_work(1)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(tau)) then
        local_tau => tau
    else
        allocate(local_tau(min(m,n)), stat=allocation_status)
    end if
    ! Retrieve work array size
    lwork = -1
    call cgeqrf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call cgeqrf(m,n,a,lda,local_tau,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(tau)) then
        deallocate(local_tau, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cgeqrf', -local_info)
    end if
end subroutine
pure subroutine mfi_zgeqrf(a, tau, info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(wp), pointer :: local_tau(:), work(:)
    complex(wp), target  :: s_work(1)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(tau)) then
        local_tau => tau
    else
        allocate(local_tau(min(m,n)), stat=allocation_status)
    end if
    ! Retrieve work array size
    lwork = -1
    call zgeqrf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call zgeqrf(m,n,a,lda,local_tau,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(tau)) then
        deallocate(local_tau, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zgeqrf', -local_info)
    end if
end subroutine
pure subroutine mfi_sgerqf(a, tau, info)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(wp), pointer :: local_tau(:), work(:)
    real(wp), target  :: s_work(1)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(tau)) then
        local_tau => tau
    else
        allocate(local_tau(min(m,n)), stat=allocation_status)
    end if
    ! Retrieve work array size
    lwork = -1
    call sgerqf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call sgerqf(m,n,a,lda,local_tau,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(tau)) then
        deallocate(local_tau, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sgerqf', -local_info)
    end if
end subroutine
pure subroutine mfi_dgerqf(a, tau, info)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(wp), pointer :: local_tau(:), work(:)
    real(wp), target  :: s_work(1)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(tau)) then
        local_tau => tau
    else
        allocate(local_tau(min(m,n)), stat=allocation_status)
    end if
    ! Retrieve work array size
    lwork = -1
    call dgerqf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call dgerqf(m,n,a,lda,local_tau,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(tau)) then
        deallocate(local_tau, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dgerqf', -local_info)
    end if
end subroutine
pure subroutine mfi_cgerqf(a, tau, info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(wp), pointer :: local_tau(:), work(:)
    complex(wp), target  :: s_work(1)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(tau)) then
        local_tau => tau
    else
        allocate(local_tau(min(m,n)), stat=allocation_status)
    end if
    ! Retrieve work array size
    lwork = -1
    call cgerqf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call cgerqf(m,n,a,lda,local_tau,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(tau)) then
        deallocate(local_tau, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cgerqf', -local_info)
    end if
end subroutine
pure subroutine mfi_zgerqf(a, tau, info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(wp), pointer :: local_tau(:), work(:)
    complex(wp), target  :: s_work(1)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(tau)) then
        local_tau => tau
    else
        allocate(local_tau(min(m,n)), stat=allocation_status)
    end if
    ! Retrieve work array size
    lwork = -1
    call zgerqf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call zgerqf(m,n,a,lda,local_tau,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(tau)) then
        deallocate(local_tau, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zgerqf', -local_info)
    end if
end subroutine
pure subroutine mfi_sgetrf(a, ipiv, info)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(:,:)
    integer, intent(out), optional, target :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, allocation_status, deallocation_status
    integer, pointer :: local_ipiv(:)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(ipiv)) then
        local_ipiv => ipiv
    else
        allocate(local_ipiv(min(m,n)), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call sgetrf(m,n,a,lda,local_ipiv,local_info)
    else
        local_info = -1000
    end if
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sgetrf', -local_info)
    end if
end subroutine
pure subroutine mfi_dgetrf(a, ipiv, info)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(:,:)
    integer, intent(out), optional, target :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, allocation_status, deallocation_status
    integer, pointer :: local_ipiv(:)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(ipiv)) then
        local_ipiv => ipiv
    else
        allocate(local_ipiv(min(m,n)), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call dgetrf(m,n,a,lda,local_ipiv,local_info)
    else
        local_info = -1000
    end if
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dgetrf', -local_info)
    end if
end subroutine
pure subroutine mfi_cgetrf(a, ipiv, info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    integer, intent(out), optional, target :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, allocation_status, deallocation_status
    integer, pointer :: local_ipiv(:)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(ipiv)) then
        local_ipiv => ipiv
    else
        allocate(local_ipiv(min(m,n)), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call cgetrf(m,n,a,lda,local_ipiv,local_info)
    else
        local_info = -1000
    end if
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cgetrf', -local_info)
    end if
end subroutine
pure subroutine mfi_zgetrf(a, ipiv, info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    integer, intent(out), optional, target :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, allocation_status, deallocation_status
    integer, pointer :: local_ipiv(:)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    allocation_status = 0
    if (present(ipiv)) then
        local_ipiv => ipiv
    else
        allocate(local_ipiv(min(m,n)), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call zgetrf(m,n,a,lda,local_ipiv,local_info)
    else
        local_info = -1000
    end if
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zgetrf', -local_info)
    end if
end subroutine
pure subroutine mfi_sgetri(a, ipiv, info)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(:,:)
    integer, intent(in) :: ipiv(:)
    real(wp), pointer :: work(:)
    real(wp) :: s_work(1)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, lwork, allocation_status, deallocation_status
    lda = max(1,size(a,1))
    n = size(a,2)
    lwork = -1
    call sgetri(n,a,lda,ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call sgetri(n,a,lda,ipiv,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sgetri',-local_info)
    end if
end subroutine
pure subroutine mfi_dgetri(a, ipiv, info)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(:,:)
    integer, intent(in) :: ipiv(:)
    real(wp), pointer :: work(:)
    real(wp) :: s_work(1)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, lwork, allocation_status, deallocation_status
    lda = max(1,size(a,1))
    n = size(a,2)
    lwork = -1
    call dgetri(n,a,lda,ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call dgetri(n,a,lda,ipiv,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dgetri',-local_info)
    end if
end subroutine
pure subroutine mfi_cgetri(a, ipiv, info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    integer, intent(in) :: ipiv(:)
    complex(wp), pointer :: work(:)
    complex(wp) :: s_work(1)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, lwork, allocation_status, deallocation_status
    lda = max(1,size(a,1))
    n = size(a,2)
    lwork = -1
    call cgetri(n,a,lda,ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call cgetri(n,a,lda,ipiv,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cgetri',-local_info)
    end if
end subroutine
pure subroutine mfi_zgetri(a, ipiv, info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    integer, intent(in) :: ipiv(:)
    complex(wp), pointer :: work(:)
    complex(wp) :: s_work(1)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, lwork, allocation_status, deallocation_status
    lda = max(1,size(a,1))
    n = size(a,2)
    lwork = -1
    call zgetri(n,a,lda,ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call zgetri(n,a,lda,ipiv,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zgetri',-local_info)
    end if
end subroutine
pure subroutine mfi_sgetrs(a,ipiv,b,trans,info)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(inout) :: b(:,:)
    integer, intent(in) :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    character, intent(in), optional :: trans
    character :: local_trans
    integer :: n, nrhs, lda, ldb
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call sgetrs(local_trans,n,nrhs,a,lda,ipiv,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sgetrs',-local_info)
    end if
end subroutine
pure subroutine mfi_dgetrs(a,ipiv,b,trans,info)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(inout) :: b(:,:)
    integer, intent(in) :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    character, intent(in), optional :: trans
    character :: local_trans
    integer :: n, nrhs, lda, ldb
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call dgetrs(local_trans,n,nrhs,a,lda,ipiv,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dgetrs',-local_info)
    end if
end subroutine
pure subroutine mfi_cgetrs(a,ipiv,b,trans,info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    integer, intent(in) :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    character, intent(in), optional :: trans
    character :: local_trans
    integer :: n, nrhs, lda, ldb
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call cgetrs(local_trans,n,nrhs,a,lda,ipiv,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cgetrs',-local_info)
    end if
end subroutine
pure subroutine mfi_zgetrs(a,ipiv,b,trans,info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    integer, intent(in) :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    character, intent(in), optional :: trans
    character :: local_trans
    integer :: n, nrhs, lda, ldb
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call zgetrs(local_trans,n,nrhs,a,lda,ipiv,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zgetrs',-local_info)
    end if
end subroutine
pure subroutine mfi_chetrf(a, uplo, ipiv, info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    integer, intent(out), optional, target :: ipiv(:)
    integer, pointer :: local_ipiv(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer   :: n, lda, lwork, allocation_status, deallocation_status
    complex(wp), target :: s_work(1)
    complex(wp), pointer :: work(:)
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    allocation_status = 0
    if (present(ipiv)) then
        local_ipiv => ipiv
    else
        allocate(local_ipiv(n), stat=allocation_status)
    end if
    lwork = -1
    call chetrf(local_uplo,n,a,lda,local_ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (.not. present(ipiv)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('chetrf',-local_info)
    end if
end subroutine
pure subroutine mfi_zhetrf(a, uplo, ipiv, info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    integer, intent(out), optional, target :: ipiv(:)
    integer, pointer :: local_ipiv(:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer   :: n, lda, lwork, allocation_status, deallocation_status
    complex(wp), target :: s_work(1)
    complex(wp), pointer :: work(:)
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    allocation_status = 0
    if (present(ipiv)) then
        local_ipiv => ipiv
    else
        allocate(local_ipiv(n), stat=allocation_status)
    end if
    lwork = -1
    call zhetrf(local_uplo,n,a,lda,local_ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (.not. present(ipiv)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zhetrf',-local_info)
    end if
end subroutine
pure subroutine mfi_chegv(a, b, w, itype, jobz, uplo, info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    real(wp), intent(out) :: w(:)
    integer, intent(in), optional :: itype
    integer :: local_itype
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    complex(wp),      pointer :: work(:)
    real(wp), pointer :: rwork(:)
    complex(wp) :: s_work(1)
    integer :: n, lda, ldb, lwork, allocation_status, deallocation_status
    if (present(itype)) then
        local_itype = itype
    else
        local_itype = 1
    end if
    if (present(jobz)) then
        local_jobz = jobz
    else
        local_jobz = 'N'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    allocation_status = 0
    allocate(rwork(max(1,3*N-2)), stat=allocation_status)
    lwork = -1
    call chegv(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,s_work,lwork,rwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call chegv(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,work,lwork,rwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    deallocate(rwork, stat=deallocation_status)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('chegv', -local_info)
    end if
end subroutine
pure subroutine mfi_zhegv(a, b, w, itype, jobz, uplo, info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    real(wp), intent(out) :: w(:)
    integer, intent(in), optional :: itype
    integer :: local_itype
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    complex(wp),      pointer :: work(:)
    real(wp), pointer :: rwork(:)
    complex(wp) :: s_work(1)
    integer :: n, lda, ldb, lwork, allocation_status, deallocation_status
    if (present(itype)) then
        local_itype = itype
    else
        local_itype = 1
    end if
    if (present(jobz)) then
        local_jobz = jobz
    else
        local_jobz = 'N'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    allocation_status = 0
    allocate(rwork(max(1,3*N-2)), stat=allocation_status)
    lwork = -1
    call zhegv(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,s_work,lwork,rwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call zhegv(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,work,lwork,rwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    deallocate(rwork, stat=deallocation_status)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zhegv', -local_info)
    end if
end subroutine
pure subroutine mfi_cheevd(a, w, jobz, uplo, info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    real(wp), intent(out) :: w(:)
    integer, intent(out), optional :: info
    integer :: local_info
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp),      pointer :: work(:)
    real(wp), pointer :: rwork(:)
    integer,       pointer :: iwork(:)
    complex(wp)      :: s_work(1)
    real(wp) :: s_rwork(1)
    integer       :: s_iwork(1)
    integer :: n, lda, lwork, lrwork, liwork, allocation_status, deallocation_status
    if (present(jobz)) then
        local_jobz = jobz
    else
        local_jobz = 'N'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    lwork  = -1
    lrwork = -1
    liwork = -1

    call cheevd(local_jobz,local_uplo,n,a,lda,w, &
                      s_work,lwork,s_rwork,lrwork,s_iwork,liwork,local_info)
    if (local_info /= 0) goto 404
    lwork  = int(s_work(1))
    lrwork = int(s_rwork(1))
    liwork = int(s_iwork(1))

    allocate(iwork(liwork), stat=allocation_status)

    if (allocation_status == 0) then
        allocate(rwork(lrwork), stat=allocation_status)
        allocate(work(lwork),   stat=allocation_status)
        call cheevd(local_jobz,local_uplo,n,a,lda,w, &
                      work,lwork,rwork,lrwork,iwork,liwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(iwork, stat=deallocation_status)
    deallocate(rwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cheevd', -local_info)
    end if
end subroutine
pure subroutine mfi_zheevd(a, w, jobz, uplo, info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    real(wp), intent(out) :: w(:)
    integer, intent(out), optional :: info
    integer :: local_info
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(wp),      pointer :: work(:)
    real(wp), pointer :: rwork(:)
    integer,       pointer :: iwork(:)
    complex(wp)      :: s_work(1)
    real(wp) :: s_rwork(1)
    integer       :: s_iwork(1)
    integer :: n, lda, lwork, lrwork, liwork, allocation_status, deallocation_status
    if (present(jobz)) then
        local_jobz = jobz
    else
        local_jobz = 'N'
    end if
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    lwork  = -1
    lrwork = -1
    liwork = -1

    call zheevd(local_jobz,local_uplo,n,a,lda,w, &
                      s_work,lwork,s_rwork,lrwork,s_iwork,liwork,local_info)
    if (local_info /= 0) goto 404
    lwork  = int(s_work(1))
    lrwork = int(s_rwork(1))
    liwork = int(s_iwork(1))

    allocate(iwork(liwork), stat=allocation_status)

    if (allocation_status == 0) then
        allocate(rwork(lrwork), stat=allocation_status)
        allocate(work(lwork),   stat=allocation_status)
        call zheevd(local_jobz,local_uplo,n,a,lda,w, &
                      work,lwork,rwork,lrwork,iwork,liwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(iwork, stat=deallocation_status)
    deallocate(rwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zheevd', -local_info)
    end if
end subroutine
pure subroutine mfi_sgesvd(a, s, u, vt, ww, job, info)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(out) :: s(:)
    real(wp),      intent(out), optional, target :: u(:,:), vt(:,:)
    real(wp), intent(out), optional, target :: ww(:)
    character, intent(in), optional :: job
    character :: local_job
    integer, intent(out), optional :: info
    integer :: local_info
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    real(wp),      target  :: s_work(1), l_a2(1,1)
    real(wp),      pointer :: local_u(:,:), local_vt(:,:), work(:)
    if (present(job)) then
        local_job = job
    else
        local_job = 'N'
    end if
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    if (present(u)) then
        ldu = max(1,size(u,1))
    else
        ldu = 1
    end if
    if (present(vt)) then
        ldvt = max(1,size(vt,1))
    else
        ldvt = 1
    end if
    if (present(u)) then
        if (size(u,2) == m) then
            jobu = 'A'
        else
            jobu = 'S'
        end if
        local_u => u
    else
        if (local_job == 'u' .or. local_job == 'U') then
            jobu = 'O'
        else
            jobu = 'N'
        end if
        local_u => l_a2
    end if
    if (present(vt)) then
        if (size(vt,1) == n) then
            jobvt = 'A'
        else
            jobvt = 'S'
        end if
        local_vt => vt
    else
        if (local_job == 'v' .or. local_job == 'V') then
            jobvt = 'O'
        else
            jobvt = 'N'
        end if
        local_vt => l_a2
    end if
    allocation_status = 0
    lwork = -1
    call sgesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,local_info)
    if (local_info /= 0) then
        goto 404
    end if
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call sgesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,local_info)
    else
        local_info = -1000
    end if

    if (present(ww)) then
        ww = real(work(2:min(m,n)-1))
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sgesvd', -local_info)
    end if
end subroutine
pure subroutine mfi_dgesvd(a, s, u, vt, ww, job, info)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(:,:)
    real(wp), intent(out) :: s(:)
    real(wp),      intent(out), optional, target :: u(:,:), vt(:,:)
    real(wp), intent(out), optional, target :: ww(:)
    character, intent(in), optional :: job
    character :: local_job
    integer, intent(out), optional :: info
    integer :: local_info
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    real(wp),      target  :: s_work(1), l_a2(1,1)
    real(wp),      pointer :: local_u(:,:), local_vt(:,:), work(:)
    if (present(job)) then
        local_job = job
    else
        local_job = 'N'
    end if
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    if (present(u)) then
        ldu = max(1,size(u,1))
    else
        ldu = 1
    end if
    if (present(vt)) then
        ldvt = max(1,size(vt,1))
    else
        ldvt = 1
    end if
    if (present(u)) then
        if (size(u,2) == m) then
            jobu = 'A'
        else
            jobu = 'S'
        end if
        local_u => u
    else
        if (local_job == 'u' .or. local_job == 'U') then
            jobu = 'O'
        else
            jobu = 'N'
        end if
        local_u => l_a2
    end if
    if (present(vt)) then
        if (size(vt,1) == n) then
            jobvt = 'A'
        else
            jobvt = 'S'
        end if
        local_vt => vt
    else
        if (local_job == 'v' .or. local_job == 'V') then
            jobvt = 'O'
        else
            jobvt = 'N'
        end if
        local_vt => l_a2
    end if
    allocation_status = 0
    lwork = -1
    call dgesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,local_info)
    if (local_info /= 0) then
        goto 404
    end if
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call dgesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,local_info)
    else
        local_info = -1000
    end if

    if (present(ww)) then
        ww = real(work(2:min(m,n)-1))
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dgesvd', -local_info)
    end if
end subroutine
pure subroutine mfi_cgesvd(a, s, u, vt, ww, job, info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    real(wp), intent(out) :: s(:)
    complex(wp),      intent(out), optional, target :: u(:,:), vt(:,:)
    real(wp), intent(out), optional, target :: ww(:)
    character, intent(in), optional :: job
    character :: local_job
    integer, intent(out), optional :: info
    integer :: local_info
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    complex(wp),      target  :: s_work(1), l_a2(1,1)
    complex(wp),      pointer :: local_u(:,:), local_vt(:,:), work(:)
    real(wp), pointer :: rwork(:)
    if (present(job)) then
        local_job = job
    else
        local_job = 'N'
    end if
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    if (present(u)) then
        ldu = max(1,size(u,1))
    else
        ldu = 1
    end if
    if (present(vt)) then
        ldvt = max(1,size(vt,1))
    else
        ldvt = 1
    end if
    if (present(u)) then
        if (size(u,2) == m) then
            jobu = 'A'
        else
            jobu = 'S'
        end if
        local_u => u
    else
        if (local_job == 'u' .or. local_job == 'U') then
            jobu = 'O'
        else
            jobu = 'N'
        end if
        local_u => l_a2
    end if
    if (present(vt)) then
        if (size(vt,1) == n) then
            jobvt = 'A'
        else
            jobvt = 'S'
        end if
        local_vt => vt
    else
        if (local_job == 'v' .or. local_job == 'V') then
            jobvt = 'O'
        else
            jobvt = 'N'
        end if
        local_vt => l_a2
    end if
    allocation_status = 0
    lwork = -1
    allocate(rwork(5*min(m,n)), stat=allocation_status)
    call cgesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,rwork,local_info)
    if (local_info /= 0) then
        goto 404
    end if
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call cgesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,rwork,local_info)
    else
        local_info = -1000
    end if

    if (present(ww)) then
        ww = real(work(2:min(m,n)-1))
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    deallocate(rwork, stat=deallocation_status)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cgesvd', -local_info)
    end if
end subroutine
pure subroutine mfi_zgesvd(a, s, u, vt, ww, job, info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    real(wp), intent(out) :: s(:)
    complex(wp),      intent(out), optional, target :: u(:,:), vt(:,:)
    real(wp), intent(out), optional, target :: ww(:)
    character, intent(in), optional :: job
    character :: local_job
    integer, intent(out), optional :: info
    integer :: local_info
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    complex(wp),      target  :: s_work(1), l_a2(1,1)
    complex(wp),      pointer :: local_u(:,:), local_vt(:,:), work(:)
    real(wp), pointer :: rwork(:)
    if (present(job)) then
        local_job = job
    else
        local_job = 'N'
    end if
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    if (present(u)) then
        ldu = max(1,size(u,1))
    else
        ldu = 1
    end if
    if (present(vt)) then
        ldvt = max(1,size(vt,1))
    else
        ldvt = 1
    end if
    if (present(u)) then
        if (size(u,2) == m) then
            jobu = 'A'
        else
            jobu = 'S'
        end if
        local_u => u
    else
        if (local_job == 'u' .or. local_job == 'U') then
            jobu = 'O'
        else
            jobu = 'N'
        end if
        local_u => l_a2
    end if
    if (present(vt)) then
        if (size(vt,1) == n) then
            jobvt = 'A'
        else
            jobvt = 'S'
        end if
        local_vt => vt
    else
        if (local_job == 'v' .or. local_job == 'V') then
            jobvt = 'O'
        else
            jobvt = 'N'
        end if
        local_vt => l_a2
    end if
    allocation_status = 0
    lwork = -1
    allocate(rwork(5*min(m,n)), stat=allocation_status)
    call zgesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,rwork,local_info)
    if (local_info /= 0) then
        goto 404
    end if
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call zgesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,rwork,local_info)
    else
        local_info = -1000
    end if

    if (present(ww)) then
        ww = real(work(2:min(m,n)-1))
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    deallocate(rwork, stat=deallocation_status)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zgesvd', -local_info)
    end if
end subroutine
pure subroutine mfi_spotrf(a, info, uplo)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call spotrf(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('spotrf', local_info)
    end if
end subroutine
pure subroutine mfi_dpotrf(a, info, uplo)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call dpotrf(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('dpotrf', local_info)
    end if
end subroutine
pure subroutine mfi_cpotrf(a, info, uplo)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call cpotrf(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('cpotrf', local_info)
    end if
end subroutine
pure subroutine mfi_zpotrf(a, info, uplo)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call zpotrf(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('zpotrf', local_info)
    end if
end subroutine
pure subroutine mfi_spotri(a, info, uplo)
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call spotri(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('spotri', local_info)
    end if
end subroutine
pure subroutine mfi_dpotri(a, info, uplo)
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call dpotri(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('dpotri', local_info)
    end if
end subroutine
pure subroutine mfi_cpotri(a, info, uplo)
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call cpotri(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('cpotri', local_info)
    end if
end subroutine
pure subroutine mfi_zpotri(a, info, uplo)
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n = size(a,2)
    call zpotri(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('zpotri', local_info)
    end if
end subroutine
pure subroutine mfi_spotrs(a, b, uplo, info)
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, nrhs, lda, ldb
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call spotrs(local_uplo,n,nrhs,a,lda,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('spotrs',-local_info)
    end if
end subroutine
pure subroutine mfi_dpotrs(a, b, uplo, info)
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(:,:)
    real(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, nrhs, lda, ldb
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call dpotrs(local_uplo,n,nrhs,a,lda,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dpotrs',-local_info)
    end if
end subroutine
pure subroutine mfi_cpotrs(a, b, uplo, info)
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, nrhs, lda, ldb
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call cpotrs(local_uplo,n,nrhs,a,lda,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cpotrs',-local_info)
    end if
end subroutine
pure subroutine mfi_zpotrs(a, b, uplo, info)
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(inout) :: b(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, nrhs, lda, ldb
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call zpotrs(local_uplo,n,nrhs,a,lda,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zpotrs',-local_info)
    end if
end subroutine

    pure subroutine mfi_error(name, info)
        character(*), intent(in) :: name
        integer, intent(in) :: info
        call f77_xerbla(name, info)
    end subroutine

end module
