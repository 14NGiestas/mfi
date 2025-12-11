!> Modern fortran interfaces for LAPACK
module mfi_lapack
use iso_fortran_env
use f77_lapack
use f77_lapack, only: mfi_lartg => f77_lartg
implicit none

!> Generic modern interface for GEQRF.
!> Supports s, d, c, z.
!> See also:
!> [[f77_geqrf:sgeqrf]], [[f77_geqrf:dgeqrf]], [[f77_geqrf:cgeqrf]], [[f77_geqrf:zgeqrf]].
interface mfi_geqrf
    module procedure :: mfi_sgeqrf
    module procedure :: mfi_dgeqrf
    module procedure :: mfi_cgeqrf
    module procedure :: mfi_zgeqrf
end interface
!> Generic modern interface for GERQF.
!> Supports s, d, c, z.
!> See also:
!> [[f77_gerqf:sgerqf]], [[f77_gerqf:dgerqf]], [[f77_gerqf:cgerqf]], [[f77_gerqf:zgerqf]].
interface mfi_gerqf
    module procedure :: mfi_sgerqf
    module procedure :: mfi_dgerqf
    module procedure :: mfi_cgerqf
    module procedure :: mfi_zgerqf
end interface
!> Generic modern interface for GETRF.
!> Supports s, d, c, z.
!> See also:
!> [[f77_getrf:sgetrf]], [[f77_getrf:dgetrf]], [[f77_getrf:cgetrf]], [[f77_getrf:zgetrf]].
interface mfi_getrf
    module procedure :: mfi_sgetrf
    module procedure :: mfi_dgetrf
    module procedure :: mfi_cgetrf
    module procedure :: mfi_zgetrf
end interface
!> Generic modern interface for GETRI.
!> Supports s, d, c, z.
!> See also:
!> [[f77_getri:sgetri]], [[f77_getri:dgetri]], [[f77_getri:cgetri]], [[f77_getri:zgetri]].
interface mfi_getri
    module procedure :: mfi_sgetri
    module procedure :: mfi_dgetri
    module procedure :: mfi_cgetri
    module procedure :: mfi_zgetri
end interface
!> Generic modern interface for GETRS.
!> Supports s, d, c, z.
!> See also:
!> [[f77_getrs:sgetrs]], [[f77_getrs:dgetrs]], [[f77_getrs:cgetrs]], [[f77_getrs:zgetrs]].
interface mfi_getrs
    module procedure :: mfi_sgetrs
    module procedure :: mfi_dgetrs
    module procedure :: mfi_cgetrs
    module procedure :: mfi_zgetrs
end interface
!> Generic modern interface for HETRF.
!> Supports c, z.
!> See also:
!> [[f77_hetrf:chetrf]], [[f77_hetrf:zhetrf]].
interface mfi_hetrf
    module procedure :: mfi_chetrf
    module procedure :: mfi_zhetrf
end interface
!> Generic modern interface for HEGV.
!> Supports c, z.
!> See also:
!> [[f77_hegv:chegv]], [[f77_hegv:zhegv]].
interface mfi_hegv
    module procedure :: mfi_chegv
    module procedure :: mfi_zhegv
end interface
!> Generic modern interface for HEEVD.
!> Supports c, z.
!> See also:
!> [[f77_heevd:cheevd]], [[f77_heevd:zheevd]].
interface mfi_heevd
    module procedure :: mfi_cheevd
    module procedure :: mfi_zheevd
end interface
!> Generic modern interface for GESVD.
!> Supports s, d, c, z.
!> See also:
!> [[f77_gesvd:sgesvd]], [[f77_gesvd:dgesvd]], [[f77_gesvd:cgesvd]], [[f77_gesvd:zgesvd]].
interface mfi_gesvd
    module procedure :: mfi_sgesvd
    module procedure :: mfi_dgesvd
    module procedure :: mfi_cgesvd
    module procedure :: mfi_zgesvd
end interface
!> Generic modern interface for ORGQR.
!> Supports s, d.
!> See also:
!> [[f77_orgqr:sorgqr]], [[f77_orgqr:dorgqr]].
interface mfi_orgqr
    module procedure :: mfi_sorgqr
    module procedure :: mfi_dorgqr
end interface
!> Generic modern interface for ORGRQ.
!> Supports s, d.
!> See also:
!> [[f77_orgrq:sorgrq]], [[f77_orgrq:dorgrq]].
interface mfi_orgrq
    module procedure :: mfi_sorgrq
    module procedure :: mfi_dorgrq
end interface
!> Generic modern interface for UNGQR.
!> Supports c, z.
!> See also:
!> [[f77_ungqr:cungqr]], [[f77_ungqr:zungqr]].
interface mfi_ungqr
    module procedure :: mfi_cungqr
    module procedure :: mfi_zungqr
end interface
!> Generic modern interface for UNGRQ.
!> Supports c, z.
!> See also:
!> [[f77_ungrq:cungrq]], [[f77_ungrq:zungrq]].
interface mfi_ungrq
    module procedure :: mfi_cungrq
    module procedure :: mfi_zungrq
end interface
!> Generic modern interface for ORMQR.
!> Supports s, d.
!> See also:
!> [[f77_ormqr:sormqr]], [[f77_ormqr:dormqr]].
interface mfi_ormqr
    module procedure :: mfi_sormqr
    module procedure :: mfi_dormqr
end interface
!> Generic modern interface for ORMRQ.
!> Supports s, d.
!> See also:
!> [[f77_ormrq:sormrq]], [[f77_ormrq:dormrq]].
interface mfi_ormrq
    module procedure :: mfi_sormrq
    module procedure :: mfi_dormrq
end interface
!> Generic modern interface for UNMQR.
!> Supports c, z.
!> See also:
!> [[f77_unmqr:cunmqr]], [[f77_unmqr:zunmqr]].
interface mfi_unmqr
    module procedure :: mfi_cunmqr
    module procedure :: mfi_zunmqr
end interface
!> Generic modern interface for UNMRQ.
!> Supports c, z.
!> See also:
!> [[f77_unmrq:cunmrq]], [[f77_unmrq:zunmrq]].
interface mfi_unmrq
    module procedure :: mfi_cunmrq
    module procedure :: mfi_zunmrq
end interface
!> Generic modern interface for POTRF.
!> Supports s, d, c, z.
!> See also:
!> [[f77_potrf:spotrf]], [[f77_potrf:dpotrf]], [[f77_potrf:cpotrf]], [[f77_potrf:zpotrf]].
interface mfi_potrf
    module procedure :: mfi_spotrf
    module procedure :: mfi_dpotrf
    module procedure :: mfi_cpotrf
    module procedure :: mfi_zpotrf
end interface
!> Generic modern interface for POTRI.
!> Supports s, d, c, z.
!> See also:
!> [[f77_potri:spotri]], [[f77_potri:dpotri]], [[f77_potri:cpotri]], [[f77_potri:zpotri]].
interface mfi_potri
    module procedure :: mfi_spotri
    module procedure :: mfi_dpotri
    module procedure :: mfi_cpotri
    module procedure :: mfi_zpotri
end interface
!> Generic modern interface for POTRS.
!> Supports s, d, c, z.
!> See also:
!> [[f77_potrs:spotrs]], [[f77_potrs:dpotrs]], [[f77_potrs:cpotrs]], [[f77_potrs:zpotrs]].
interface mfi_potrs
    module procedure :: mfi_spotrs
    module procedure :: mfi_dpotrs
    module procedure :: mfi_cpotrs
    module procedure :: mfi_zpotrs
end interface
!> Generic modern interface for POCON.
!> Supports s, d, c, z.
!> See also:
!> [[f77_pocon:spocon]], [[f77_pocon:dpocon]], [[f77_pocon:cpocon]], [[f77_pocon:zpocon]].
interface mfi_pocon
    module procedure :: mfi_spocon
    module procedure :: mfi_dpocon
    module procedure :: mfi_cpocon
    module procedure :: mfi_zpocon
end interface
!> Generic modern interface for TRTRS.
!> Supports s, d, c, z.
!> See also:
!> [[f77_trtrs:strtrs]], [[f77_trtrs:dtrtrs]], [[f77_trtrs:ctrtrs]], [[f77_trtrs:ztrtrs]].
interface mfi_trtrs
    module procedure :: mfi_strtrs
    module procedure :: mfi_dtrtrs
    module procedure :: mfi_ctrtrs
    module procedure :: mfi_ztrtrs
end interface
!> Generic modern interface for SYTRF.
!> Supports s, d.
!> See also:
!> [[f77_sytrf:ssytrf]], [[f77_sytrf:dsytrf]].
interface mfi_sytrf
    module procedure :: mfi_ssytrf
    module procedure :: mfi_dsytrf
end interface

contains

!> Modern interface for [[f77_geqrf:sgeqrf]].
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
pure subroutine mfi_sgeqrf(a, tau, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(REAL32), pointer :: local_tau(:), work(:)
    real(REAL32), target  :: s_work(1)
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
!> Modern interface for [[f77_geqrf:dgeqrf]].
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
pure subroutine mfi_dgeqrf(a, tau, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(REAL64), pointer :: local_tau(:), work(:)
    real(REAL64), target  :: s_work(1)
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
!> Modern interface for [[f77_geqrf:cgeqrf]].
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
pure subroutine mfi_cgeqrf(a, tau, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(REAL32), pointer :: local_tau(:), work(:)
    complex(REAL32), target  :: s_work(1)
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
!> Modern interface for [[f77_geqrf:zgeqrf]].
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
pure subroutine mfi_zgeqrf(a, tau, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(REAL64), pointer :: local_tau(:), work(:)
    complex(REAL64), target  :: s_work(1)
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
!> Modern interface for [[f77_gerqf:sgerqf]].
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
pure subroutine mfi_sgerqf(a, tau, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(REAL32), pointer :: local_tau(:), work(:)
    real(REAL32), target  :: s_work(1)
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
!> Modern interface for [[f77_gerqf:dgerqf]].
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
pure subroutine mfi_dgerqf(a, tau, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(REAL64), pointer :: local_tau(:), work(:)
    real(REAL64), target  :: s_work(1)
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
!> Modern interface for [[f77_gerqf:cgerqf]].
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
pure subroutine mfi_cgerqf(a, tau, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(REAL32), pointer :: local_tau(:), work(:)
    complex(REAL32), target  :: s_work(1)
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
!> Modern interface for [[f77_gerqf:zgerqf]].
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
pure subroutine mfi_zgerqf(a, tau, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(out), optional, target :: tau(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(REAL64), pointer :: local_tau(:), work(:)
    complex(REAL64), target  :: s_work(1)
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
!> Modern interface for [[f77_getrf:sgetrf]].
!> See also: [[mfi_getrf]], [[f77_getrf]].
pure subroutine mfi_sgetrf(a, ipiv, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_getrf:dgetrf]].
!> See also: [[mfi_getrf]], [[f77_getrf]].
pure subroutine mfi_dgetrf(a, ipiv, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_getrf:cgetrf]].
!> See also: [[mfi_getrf]], [[f77_getrf]].
pure subroutine mfi_cgetrf(a, ipiv, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_getrf:zgetrf]].
!> See also: [[mfi_getrf]], [[f77_getrf]].
pure subroutine mfi_zgetrf(a, ipiv, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_getri:sgetri]].
!> See also: [[mfi_getri]], [[f77_getri]].
pure subroutine mfi_sgetri(a, ipiv, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    integer, intent(in) :: ipiv(:)
    real(REAL32), pointer :: work(:)
    real(REAL32) :: s_work(1)
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
!> Modern interface for [[f77_getri:dgetri]].
!> See also: [[mfi_getri]], [[f77_getri]].
pure subroutine mfi_dgetri(a, ipiv, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    integer, intent(in) :: ipiv(:)
    real(REAL64), pointer :: work(:)
    real(REAL64) :: s_work(1)
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
!> Modern interface for [[f77_getri:cgetri]].
!> See also: [[mfi_getri]], [[f77_getri]].
pure subroutine mfi_cgetri(a, ipiv, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    integer, intent(in) :: ipiv(:)
    complex(REAL32), pointer :: work(:)
    complex(REAL32) :: s_work(1)
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
!> Modern interface for [[f77_getri:zgetri]].
!> See also: [[mfi_getri]], [[f77_getri]].
pure subroutine mfi_zgetri(a, ipiv, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    integer, intent(in) :: ipiv(:)
    complex(REAL64), pointer :: work(:)
    complex(REAL64) :: s_work(1)
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
!> Modern interface for [[f77_getrs:sgetrs]].
!> See also: [[mfi_getrs]], [[f77_getrs]].
pure subroutine mfi_sgetrs(a,ipiv,b,trans,info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(inout) :: b(:,:)
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
!> Modern interface for [[f77_getrs:dgetrs]].
!> See also: [[mfi_getrs]], [[f77_getrs]].
pure subroutine mfi_dgetrs(a,ipiv,b,trans,info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(inout) :: b(:,:)
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
!> Modern interface for [[f77_getrs:cgetrs]].
!> See also: [[mfi_getrs]], [[f77_getrs]].
pure subroutine mfi_cgetrs(a,ipiv,b,trans,info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(inout) :: b(:,:)
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
!> Modern interface for [[f77_getrs:zgetrs]].
!> See also: [[mfi_getrs]], [[f77_getrs]].
pure subroutine mfi_zgetrs(a,ipiv,b,trans,info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(inout) :: b(:,:)
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
!> Modern interface for [[f77_hetrf:chetrf]].
!> See also: [[mfi_hetrf]], [[f77_hetrf]].
!> Computes the factorization of a Hermitian matrix using the Bunch-Kaufman diagonal pivoting method
!>
!> The factorization has the form:
!> - A = U*D*U**H (if uplo='U') or
!> - A = L*D*L**H (if uplo='L')
!>
!> where U (or L) is a product of permutation and unit upper (lower) triangular matrices,
!> and D is block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> Parameters:
!> - `a` (inout): On entry, the Hermitian matrix A. On exit, the block diagonal matrix D
!>                and the multipliers used to obtain the factor U or L.
!> - `uplo` (in, optional): Specifies whether the upper ('U') or lower ('L') triangular part
!>                of the Hermitian matrix A is stored. Default: 'U'
!> - `ipiv` (out, optional): The pivot indices that define the permutation matrix P.
!>                If ipiv is not provided, it will be allocated internally.
!> - `info` (out, optional): Output status: 0 for success, < 0 for illegal argument,
!>                > 0 if D(k,k) is exactly zero.
pure subroutine mfi_chetrf(a, uplo, ipiv, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional, target :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, lwork, allocation_status, deallocation_status
    integer, pointer :: local_ipiv(:)
    complex(REAL32), pointer :: work(:)
    complex(REAL32) :: s_work(1)  ! Work array for workspace query
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

    ! Retrieve work array size
    lwork = -1
    call chetrf(local_uplo, n, a, lda, local_ipiv, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call chetrf(local_uplo, n, a, lda, local_ipiv, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('chetrf', -local_info)
    end if

    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        if (local_info <= -1000) then
            call mfi_error('chetrf', -local_info)
        else
            call mfi_error('chetrf', local_info)
        end if
    end if
end subroutine
!> Modern interface for [[f77_hetrf:zhetrf]].
!> See also: [[mfi_hetrf]], [[f77_hetrf]].
!> Computes the factorization of a Hermitian matrix using the Bunch-Kaufman diagonal pivoting method
!>
!> The factorization has the form:
!> - A = U*D*U**H (if uplo='U') or
!> - A = L*D*L**H (if uplo='L')
!>
!> where U (or L) is a product of permutation and unit upper (lower) triangular matrices,
!> and D is block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> Parameters:
!> - `a` (inout): On entry, the Hermitian matrix A. On exit, the block diagonal matrix D
!>                and the multipliers used to obtain the factor U or L.
!> - `uplo` (in, optional): Specifies whether the upper ('U') or lower ('L') triangular part
!>                of the Hermitian matrix A is stored. Default: 'U'
!> - `ipiv` (out, optional): The pivot indices that define the permutation matrix P.
!>                If ipiv is not provided, it will be allocated internally.
!> - `info` (out, optional): Output status: 0 for success, < 0 for illegal argument,
!>                > 0 if D(k,k) is exactly zero.
pure subroutine mfi_zhetrf(a, uplo, ipiv, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional, target :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, lwork, allocation_status, deallocation_status
    integer, pointer :: local_ipiv(:)
    complex(REAL64), pointer :: work(:)
    complex(REAL64) :: s_work(1)  ! Work array for workspace query
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

    ! Retrieve work array size
    lwork = -1
    call zhetrf(local_uplo, n, a, lda, local_ipiv, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call zhetrf(local_uplo, n, a, lda, local_ipiv, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zhetrf', -local_info)
    end if

    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        if (local_info <= -1000) then
            call mfi_error('zhetrf', -local_info)
        else
            call mfi_error('zhetrf', local_info)
        end if
    end if
end subroutine
!> Modern interface for [[f77_hegv:chegv]].
!> See also: [[mfi_hegv]], [[f77_hegv]].
pure subroutine mfi_chegv(a, b, w, itype, jobz, uplo, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(inout) :: b(:,:)
    real(REAL32), intent(out) :: w(:)
    integer, intent(in), optional :: itype
    integer :: local_itype
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    complex(REAL32),      pointer :: work(:)
    real(REAL32), pointer :: rwork(:)
    complex(REAL32) :: s_work(1)
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
!> Modern interface for [[f77_hegv:zhegv]].
!> See also: [[mfi_hegv]], [[f77_hegv]].
pure subroutine mfi_zhegv(a, b, w, itype, jobz, uplo, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(inout) :: b(:,:)
    real(REAL64), intent(out) :: w(:)
    integer, intent(in), optional :: itype
    integer :: local_itype
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    complex(REAL64),      pointer :: work(:)
    real(REAL64), pointer :: rwork(:)
    complex(REAL64) :: s_work(1)
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
!> Modern interface for [[f77_heevd:cheevd]].
!> See also: [[mfi_heevd]], [[f77_heevd]].
pure subroutine mfi_cheevd(a, w, jobz, uplo, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(out) :: w(:)
    integer, intent(out), optional :: info
    integer :: local_info
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL32),      pointer :: work(:)
    real(REAL32), pointer :: rwork(:)
    integer,       pointer :: iwork(:)
    complex(REAL32)      :: s_work(1)
    real(REAL32) :: s_rwork(1)
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
!> Modern interface for [[f77_heevd:zheevd]].
!> See also: [[mfi_heevd]], [[f77_heevd]].
pure subroutine mfi_zheevd(a, w, jobz, uplo, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(out) :: w(:)
    integer, intent(out), optional :: info
    integer :: local_info
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    complex(REAL64),      pointer :: work(:)
    real(REAL64), pointer :: rwork(:)
    integer,       pointer :: iwork(:)
    complex(REAL64)      :: s_work(1)
    real(REAL64) :: s_rwork(1)
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
!> Modern interface for [[f77_gesvd:sgesvd]].
!> See also: [[mfi_gesvd]], [[f77_gesvd]].
pure subroutine mfi_sgesvd(a, s, u, vt, ww, job, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(out) :: s(:)
    real(REAL32),      intent(out), optional, target :: u(:,:), vt(:,:)
    real(REAL32), intent(out), optional, target :: ww(:)
    character, intent(in), optional :: job
    character :: local_job
    integer, intent(out), optional :: info
    integer :: local_info
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    real(REAL32),      target  :: s_work(1), l_a2(1,1)
    real(REAL32),      pointer :: local_u(:,:), local_vt(:,:), work(:)
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
!> Modern interface for [[f77_gesvd:dgesvd]].
!> See also: [[mfi_gesvd]], [[f77_gesvd]].
pure subroutine mfi_dgesvd(a, s, u, vt, ww, job, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(out) :: s(:)
    real(REAL64),      intent(out), optional, target :: u(:,:), vt(:,:)
    real(REAL64), intent(out), optional, target :: ww(:)
    character, intent(in), optional :: job
    character :: local_job
    integer, intent(out), optional :: info
    integer :: local_info
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    real(REAL64),      target  :: s_work(1), l_a2(1,1)
    real(REAL64),      pointer :: local_u(:,:), local_vt(:,:), work(:)
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
!> Modern interface for [[f77_gesvd:cgesvd]].
!> See also: [[mfi_gesvd]], [[f77_gesvd]].
pure subroutine mfi_cgesvd(a, s, u, vt, ww, job, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(out) :: s(:)
    complex(REAL32),      intent(out), optional, target :: u(:,:), vt(:,:)
    real(REAL32), intent(out), optional, target :: ww(:)
    character, intent(in), optional :: job
    character :: local_job
    integer, intent(out), optional :: info
    integer :: local_info
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    complex(REAL32),      target  :: s_work(1), l_a2(1,1)
    complex(REAL32),      pointer :: local_u(:,:), local_vt(:,:), work(:)
    real(REAL32), pointer :: rwork(:)
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
!> Modern interface for [[f77_gesvd:zgesvd]].
!> See also: [[mfi_gesvd]], [[f77_gesvd]].
pure subroutine mfi_zgesvd(a, s, u, vt, ww, job, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(out) :: s(:)
    complex(REAL64),      intent(out), optional, target :: u(:,:), vt(:,:)
    real(REAL64), intent(out), optional, target :: ww(:)
    character, intent(in), optional :: job
    character :: local_job
    integer, intent(out), optional :: info
    integer :: local_info
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    complex(REAL64),      target  :: s_work(1), l_a2(1,1)
    complex(REAL64),      pointer :: local_u(:,:), local_vt(:,:), work(:)
    real(REAL64), pointer :: rwork(:)
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
!> Modern interface for [[f77_orgqr:sorgqr]].
!> See also: [[mfi_orgqr]], [[f77_orgqr]].
!> Generates the real orthogonal matrix Q of the QR factorization
pure subroutine mfi_sorgqr(a, tau, k, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(REAL32), pointer :: work(:)
    real(REAL32), target  :: s_work(1)
    if (present(k)) then
        local_k = k
    else
        local_k = min(size(a,1), size(a,2))
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call sorgqr(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call sorgqr(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sorgqr', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_orgqr:dorgqr]].
!> See also: [[mfi_orgqr]], [[f77_orgqr]].
!> Generates the real orthogonal matrix Q of the QR factorization
pure subroutine mfi_dorgqr(a, tau, k, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(REAL64), pointer :: work(:)
    real(REAL64), target  :: s_work(1)
    if (present(k)) then
        local_k = k
    else
        local_k = min(size(a,1), size(a,2))
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call dorgqr(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call dorgqr(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dorgqr', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_orgrq:sorgrq]].
!> See also: [[mfi_orgrq]], [[f77_orgrq]].
!> Generates the real orthogonal matrix Q of the QR factorization
pure subroutine mfi_sorgrq(a, tau, k, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(REAL32), pointer :: work(:)
    real(REAL32), target  :: s_work(1)
    if (present(k)) then
        local_k = k
    else
        local_k = min(size(a,1), size(a,2))
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call sorgrq(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call sorgrq(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sorgrq', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_orgrq:dorgrq]].
!> See also: [[mfi_orgrq]], [[f77_orgrq]].
!> Generates the real orthogonal matrix Q of the QR factorization
pure subroutine mfi_dorgrq(a, tau, k, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    real(REAL64), pointer :: work(:)
    real(REAL64), target  :: s_work(1)
    if (present(k)) then
        local_k = k
    else
        local_k = min(size(a,1), size(a,2))
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call dorgrq(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call dorgrq(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dorgrq', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ungqr:cungqr]].
!> See also: [[mfi_ungqr]], [[f77_ungqr]].
!> Generates the real orthogonal matrix Q of the QR factorization
pure subroutine mfi_cungqr(a, tau, k, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(REAL32), pointer :: work(:)
    complex(REAL32), target  :: s_work(1)
    if (present(k)) then
        local_k = k
    else
        local_k = min(size(a,1), size(a,2))
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call cungqr(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call cungqr(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cungqr', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ungqr:zungqr]].
!> See also: [[mfi_ungqr]], [[f77_ungqr]].
!> Generates the real orthogonal matrix Q of the QR factorization
pure subroutine mfi_zungqr(a, tau, k, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(REAL64), pointer :: work(:)
    complex(REAL64), target  :: s_work(1)
    if (present(k)) then
        local_k = k
    else
        local_k = min(size(a,1), size(a,2))
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call zungqr(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call zungqr(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zungqr', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ungrq:cungrq]].
!> See also: [[mfi_ungrq]], [[f77_ungrq]].
!> Generates the real orthogonal matrix Q of the QR factorization
pure subroutine mfi_cungrq(a, tau, k, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(REAL32), pointer :: work(:)
    complex(REAL32), target  :: s_work(1)
    if (present(k)) then
        local_k = k
    else
        local_k = min(size(a,1), size(a,2))
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call cungrq(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call cungrq(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cungrq', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ungrq:zungrq]].
!> See also: [[mfi_ungrq]], [[f77_ungrq]].
!> Generates the real orthogonal matrix Q of the QR factorization
pure subroutine mfi_zungrq(a, tau, k, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda, lwork, allocation_status, deallocation_status
    complex(REAL64), pointer :: work(:)
    complex(REAL64), target  :: s_work(1)
    if (present(k)) then
        local_k = k
    else
        local_k = min(size(a,1), size(a,2))
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call zungrq(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call zungrq(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zungrq', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ormqr:sormqr]].
!> See also: [[mfi_ormqr]], [[f77_ormqr]].
!> Multiplies a matrix by the orthogonal matrix Q from geqrf
pure subroutine mfi_sormqr(a, tau, c, side, trans, info)
    integer, parameter :: wp = REAL32
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: trans
    character :: local_trans
    integer, intent(out), optional :: info
    integer :: local_info
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(in) :: tau(:)
    real(REAL32), intent(inout) :: c(:,:)
    integer :: m, n, k, lda, ldc, lwork, allocation_status, deallocation_status
    real(REAL32), pointer :: work(:)
    real(REAL32), target  :: s_work(1)
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1, size(a, 1))
    ldc = max(1, size(c, 1))
    m = size(c, 1)
    n = size(c, 2)

    if (local_side == 'L' .or. local_side == 'l') then
        k = size(a, 2)
    else
        k = size(a, 2)
    end if

    ! Query workspace size
    lwork = -1
    call sormqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call sormqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sormqr', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ormqr:dormqr]].
!> See also: [[mfi_ormqr]], [[f77_ormqr]].
!> Multiplies a matrix by the orthogonal matrix Q from geqrf
pure subroutine mfi_dormqr(a, tau, c, side, trans, info)
    integer, parameter :: wp = REAL64
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: trans
    character :: local_trans
    integer, intent(out), optional :: info
    integer :: local_info
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(in) :: tau(:)
    real(REAL64), intent(inout) :: c(:,:)
    integer :: m, n, k, lda, ldc, lwork, allocation_status, deallocation_status
    real(REAL64), pointer :: work(:)
    real(REAL64), target  :: s_work(1)
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1, size(a, 1))
    ldc = max(1, size(c, 1))
    m = size(c, 1)
    n = size(c, 2)

    if (local_side == 'L' .or. local_side == 'l') then
        k = size(a, 2)
    else
        k = size(a, 2)
    end if

    ! Query workspace size
    lwork = -1
    call dormqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call dormqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dormqr', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ormrq:sormrq]].
!> See also: [[mfi_ormrq]], [[f77_ormrq]].
!> Multiplies a matrix by the orthogonal matrix Q from geqrf
pure subroutine mfi_sormrq(a, tau, c, side, trans, info)
    integer, parameter :: wp = REAL32
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: trans
    character :: local_trans
    integer, intent(out), optional :: info
    integer :: local_info
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(in) :: tau(:)
    real(REAL32), intent(inout) :: c(:,:)
    integer :: m, n, k, lda, ldc, lwork, allocation_status, deallocation_status
    real(REAL32), pointer :: work(:)
    real(REAL32), target  :: s_work(1)
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1, size(a, 1))
    ldc = max(1, size(c, 1))
    m = size(c, 1)
    n = size(c, 2)

    if (local_side == 'L' .or. local_side == 'l') then
        k = size(a, 2)
    else
        k = size(a, 2)
    end if

    ! Query workspace size
    lwork = -1
    call sormrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call sormrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('sormrq', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ormrq:dormrq]].
!> See also: [[mfi_ormrq]], [[f77_ormrq]].
!> Multiplies a matrix by the orthogonal matrix Q from geqrf
pure subroutine mfi_dormrq(a, tau, c, side, trans, info)
    integer, parameter :: wp = REAL64
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: trans
    character :: local_trans
    integer, intent(out), optional :: info
    integer :: local_info
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(in) :: tau(:)
    real(REAL64), intent(inout) :: c(:,:)
    integer :: m, n, k, lda, ldc, lwork, allocation_status, deallocation_status
    real(REAL64), pointer :: work(:)
    real(REAL64), target  :: s_work(1)
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1, size(a, 1))
    ldc = max(1, size(c, 1))
    m = size(c, 1)
    n = size(c, 2)

    if (local_side == 'L' .or. local_side == 'l') then
        k = size(a, 2)
    else
        k = size(a, 2)
    end if

    ! Query workspace size
    lwork = -1
    call dormrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call dormrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dormrq', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_unmqr:cunmqr]].
!> See also: [[mfi_unmqr]], [[f77_unmqr]].
!> Multiplies a matrix by the orthogonal matrix Q from geqrf
pure subroutine mfi_cunmqr(a, tau, c, side, trans, info)
    integer, parameter :: wp = REAL32
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: trans
    character :: local_trans
    integer, intent(out), optional :: info
    integer :: local_info
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(in) :: tau(:)
    complex(REAL32), intent(inout) :: c(:,:)
    integer :: m, n, k, lda, ldc, lwork, allocation_status, deallocation_status
    complex(REAL32), pointer :: work(:)
    complex(REAL32), target  :: s_work(1)
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1, size(a, 1))
    ldc = max(1, size(c, 1))
    m = size(c, 1)
    n = size(c, 2)

    if (local_side == 'L' .or. local_side == 'l') then
        k = size(a, 2)
    else
        k = size(a, 2)
    end if

    ! Query workspace size
    lwork = -1
    call cunmqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call cunmqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cunmqr', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_unmqr:zunmqr]].
!> See also: [[mfi_unmqr]], [[f77_unmqr]].
!> Multiplies a matrix by the orthogonal matrix Q from geqrf
pure subroutine mfi_zunmqr(a, tau, c, side, trans, info)
    integer, parameter :: wp = REAL64
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: trans
    character :: local_trans
    integer, intent(out), optional :: info
    integer :: local_info
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(in) :: tau(:)
    complex(REAL64), intent(inout) :: c(:,:)
    integer :: m, n, k, lda, ldc, lwork, allocation_status, deallocation_status
    complex(REAL64), pointer :: work(:)
    complex(REAL64), target  :: s_work(1)
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1, size(a, 1))
    ldc = max(1, size(c, 1))
    m = size(c, 1)
    n = size(c, 2)

    if (local_side == 'L' .or. local_side == 'l') then
        k = size(a, 2)
    else
        k = size(a, 2)
    end if

    ! Query workspace size
    lwork = -1
    call zunmqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call zunmqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zunmqr', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_unmrq:cunmrq]].
!> See also: [[mfi_unmrq]], [[f77_unmrq]].
!> Multiplies a matrix by the orthogonal matrix Q from geqrf
pure subroutine mfi_cunmrq(a, tau, c, side, trans, info)
    integer, parameter :: wp = REAL32
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: trans
    character :: local_trans
    integer, intent(out), optional :: info
    integer :: local_info
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(in) :: tau(:)
    complex(REAL32), intent(inout) :: c(:,:)
    integer :: m, n, k, lda, ldc, lwork, allocation_status, deallocation_status
    complex(REAL32), pointer :: work(:)
    complex(REAL32), target  :: s_work(1)
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1, size(a, 1))
    ldc = max(1, size(c, 1))
    m = size(c, 1)
    n = size(c, 2)

    if (local_side == 'L' .or. local_side == 'l') then
        k = size(a, 2)
    else
        k = size(a, 2)
    end if

    ! Query workspace size
    lwork = -1
    call cunmrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call cunmrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cunmrq', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_unmrq:zunmrq]].
!> See also: [[mfi_unmrq]], [[f77_unmrq]].
!> Multiplies a matrix by the orthogonal matrix Q from geqrf
pure subroutine mfi_zunmrq(a, tau, c, side, trans, info)
    integer, parameter :: wp = REAL64
    character, intent(in), optional :: side
    character :: local_side
    character, intent(in), optional :: trans
    character :: local_trans
    integer, intent(out), optional :: info
    integer :: local_info
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(in) :: tau(:)
    complex(REAL64), intent(inout) :: c(:,:)
    integer :: m, n, k, lda, ldc, lwork, allocation_status, deallocation_status
    complex(REAL64), pointer :: work(:)
    complex(REAL64), target  :: s_work(1)
    if (present(side)) then
        local_side = side
    else
        local_side = 'L'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    lda = max(1, size(a, 1))
    ldc = max(1, size(c, 1))
    m = size(c, 1)
    n = size(c, 2)

    if (local_side == 'L' .or. local_side == 'l') then
        k = size(a, 2)
    else
        k = size(a, 2)
    end if

    ! Query workspace size
    lwork = -1
    call zunmrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call zunmrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zunmrq', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_potrf:spotrf]].
!> See also: [[mfi_potrf]], [[f77_potrf]].
pure subroutine mfi_spotrf(a, info, uplo)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_potrf:dpotrf]].
!> See also: [[mfi_potrf]], [[f77_potrf]].
pure subroutine mfi_dpotrf(a, info, uplo)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_potrf:cpotrf]].
!> See also: [[mfi_potrf]], [[f77_potrf]].
pure subroutine mfi_cpotrf(a, info, uplo)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_potrf:zpotrf]].
!> See also: [[mfi_potrf]], [[f77_potrf]].
pure subroutine mfi_zpotrf(a, info, uplo)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_potri:spotri]].
!> See also: [[mfi_potri]], [[f77_potri]].
pure subroutine mfi_spotri(a, info, uplo)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_potri:dpotri]].
!> See also: [[mfi_potri]], [[f77_potri]].
pure subroutine mfi_dpotri(a, info, uplo)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_potri:cpotri]].
!> See also: [[mfi_potri]], [[f77_potri]].
pure subroutine mfi_cpotri(a, info, uplo)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_potri:zpotri]].
!> See also: [[mfi_potri]], [[f77_potri]].
pure subroutine mfi_zpotri(a, info, uplo)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
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
!> Modern interface for [[f77_potrs:spotrs]].
!> See also: [[mfi_potrs]], [[f77_potrs]].
pure subroutine mfi_spotrs(a, b, uplo, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: b(:,:)
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
!> Modern interface for [[f77_potrs:dpotrs]].
!> See also: [[mfi_potrs]], [[f77_potrs]].
pure subroutine mfi_dpotrs(a, b, uplo, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: b(:,:)
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
!> Modern interface for [[f77_potrs:cpotrs]].
!> See also: [[mfi_potrs]], [[f77_potrs]].
pure subroutine mfi_cpotrs(a, b, uplo, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: b(:,:)
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
!> Modern interface for [[f77_potrs:zpotrs]].
!> See also: [[mfi_potrs]], [[f77_potrs]].
pure subroutine mfi_zpotrs(a, b, uplo, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: b(:,:)
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
!> Modern interface for [[f77_pocon:spocon]].
!> See also: [[mfi_pocon]], [[f77_pocon]].
!> Estimates the reciprocal of the condition number of a real symmetric / complex Hermitian positive definite matrix using the Cholesky factorization computed by ?POTRF
pure subroutine mfi_spocon(a, anorm, rcond, uplo, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(in) :: anorm
    real(REAL32), intent(out) :: rcond
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, allocation_status, deallocation_status
    real(REAL32), pointer :: work(:)
    integer, pointer :: xwork(:)
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    allocate(xwork(n), stat=allocation_status)
    if (allocation_status == 0) allocate(work(3*n), stat=allocation_status)

    if (allocation_status == 0) then
        call spocon(local_uplo, n, a, lda, anorm, rcond, work, xwork, local_info)
    else
        local_info = -1000
    end if

    deallocate(xwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)

    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('spocon',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_pocon:dpocon]].
!> See also: [[mfi_pocon]], [[f77_pocon]].
!> Estimates the reciprocal of the condition number of a real symmetric / complex Hermitian positive definite matrix using the Cholesky factorization computed by ?POTRF
pure subroutine mfi_dpocon(a, anorm, rcond, uplo, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(in) :: anorm
    real(REAL64), intent(out) :: rcond
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, allocation_status, deallocation_status
    real(REAL64), pointer :: work(:)
    integer, pointer :: xwork(:)
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    allocate(xwork(n), stat=allocation_status)
    if (allocation_status == 0) allocate(work(3*n), stat=allocation_status)

    if (allocation_status == 0) then
        call dpocon(local_uplo, n, a, lda, anorm, rcond, work, xwork, local_info)
    else
        local_info = -1000
    end if

    deallocate(xwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)

    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dpocon',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_pocon:cpocon]].
!> See also: [[mfi_pocon]], [[f77_pocon]].
!> Estimates the reciprocal of the condition number of a real symmetric / complex Hermitian positive definite matrix using the Cholesky factorization computed by ?POTRF
pure subroutine mfi_cpocon(a, anorm, rcond, uplo, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(in) :: anorm
    real(REAL32), intent(out) :: rcond
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, allocation_status, deallocation_status
    complex(REAL32), pointer :: work(:)
    real(REAL32), pointer :: xwork(:)
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    allocate(xwork(n), stat=allocation_status)
    if (allocation_status == 0) allocate(work(3*n), stat=allocation_status)

    if (allocation_status == 0) then
        call cpocon(local_uplo, n, a, lda, anorm, rcond, work, xwork, local_info)
    else
        local_info = -1000
    end if

    deallocate(xwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)

    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('cpocon',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_pocon:zpocon]].
!> See also: [[mfi_pocon]], [[f77_pocon]].
!> Estimates the reciprocal of the condition number of a real symmetric / complex Hermitian positive definite matrix using the Cholesky factorization computed by ?POTRF
pure subroutine mfi_zpocon(a, anorm, rcond, uplo, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(in) :: anorm
    real(REAL64), intent(out) :: rcond
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, allocation_status, deallocation_status
    complex(REAL64), pointer :: work(:)
    real(REAL64), pointer :: xwork(:)
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    allocate(xwork(n), stat=allocation_status)
    if (allocation_status == 0) allocate(work(3*n), stat=allocation_status)

    if (allocation_status == 0) then
        call zpocon(local_uplo, n, a, lda, anorm, rcond, work, xwork, local_info)
    else
        local_info = -1000
    end if

    deallocate(xwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)

    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('zpocon',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_trtrs:strtrs]].
!> See also: [[mfi_trtrs]], [[f77_trtrs]].
!> Solves a triangular linear system with multiple right-hand sides:
!>
!>     A * X = B  or  A**T * X = B  or  A**H * X = B,
!>
!> where A is a triangular matrix (stored in `a`), and B is overwritten by the solution X.
!>
!> Optional arguments:
!> - `uplo`: 'U' (upper triangular, default) or 'L' (lower triangular)
!> - `trans`: 'N' (no transpose), 'T' (transpose), or 'C' (conjugate transpose, default 'N')
!> - `diag`: 'N' (non-unit diagonal, default) or 'U' (unit diagonal)
!> - `info`: if not present and `info /= 0`, calls [[mfi_error]].
!>
!> The shapes are inferred from `a` (N-by-N) and `b` (N-by-NRHS).
pure subroutine mfi_strtrs(a, b, uplo, trans, diag, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(:,:)
    real(REAL32), intent(inout) :: b(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, nrhs, lda, ldb
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call strtrs(local_uplo, local_trans, local_diag, n, nrhs, a, lda, b, ldb, local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('strtrs', local_info)
    end if
end subroutine
!> Modern interface for [[f77_trtrs:dtrtrs]].
!> See also: [[mfi_trtrs]], [[f77_trtrs]].
!> Solves a triangular linear system with multiple right-hand sides:
!>
!>     A * X = B  or  A**T * X = B  or  A**H * X = B,
!>
!> where A is a triangular matrix (stored in `a`), and B is overwritten by the solution X.
!>
!> Optional arguments:
!> - `uplo`: 'U' (upper triangular, default) or 'L' (lower triangular)
!> - `trans`: 'N' (no transpose), 'T' (transpose), or 'C' (conjugate transpose, default 'N')
!> - `diag`: 'N' (non-unit diagonal, default) or 'U' (unit diagonal)
!> - `info`: if not present and `info /= 0`, calls [[mfi_error]].
!>
!> The shapes are inferred from `a` (N-by-N) and `b` (N-by-NRHS).
pure subroutine mfi_dtrtrs(a, b, uplo, trans, diag, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(:,:)
    real(REAL64), intent(inout) :: b(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, nrhs, lda, ldb
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call dtrtrs(local_uplo, local_trans, local_diag, n, nrhs, a, lda, b, ldb, local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('dtrtrs', local_info)
    end if
end subroutine
!> Modern interface for [[f77_trtrs:ctrtrs]].
!> See also: [[mfi_trtrs]], [[f77_trtrs]].
!> Solves a triangular linear system with multiple right-hand sides:
!>
!>     A * X = B  or  A**T * X = B  or  A**H * X = B,
!>
!> where A is a triangular matrix (stored in `a`), and B is overwritten by the solution X.
!>
!> Optional arguments:
!> - `uplo`: 'U' (upper triangular, default) or 'L' (lower triangular)
!> - `trans`: 'N' (no transpose), 'T' (transpose), or 'C' (conjugate transpose, default 'N')
!> - `diag`: 'N' (non-unit diagonal, default) or 'U' (unit diagonal)
!> - `info`: if not present and `info /= 0`, calls [[mfi_error]].
!>
!> The shapes are inferred from `a` (N-by-N) and `b` (N-by-NRHS).
pure subroutine mfi_ctrtrs(a, b, uplo, trans, diag, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(:,:)
    complex(REAL32), intent(inout) :: b(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, nrhs, lda, ldb
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call ctrtrs(local_uplo, local_trans, local_diag, n, nrhs, a, lda, b, ldb, local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('ctrtrs', local_info)
    end if
end subroutine
!> Modern interface for [[f77_trtrs:ztrtrs]].
!> See also: [[mfi_trtrs]], [[f77_trtrs]].
!> Solves a triangular linear system with multiple right-hand sides:
!>
!>     A * X = B  or  A**T * X = B  or  A**H * X = B,
!>
!> where A is a triangular matrix (stored in `a`), and B is overwritten by the solution X.
!>
!> Optional arguments:
!> - `uplo`: 'U' (upper triangular, default) or 'L' (lower triangular)
!> - `trans`: 'N' (no transpose), 'T' (transpose), or 'C' (conjugate transpose, default 'N')
!> - `diag`: 'N' (non-unit diagonal, default) or 'U' (unit diagonal)
!> - `info`: if not present and `info /= 0`, calls [[mfi_error]].
!>
!> The shapes are inferred from `a` (N-by-N) and `b` (N-by-NRHS).
pure subroutine mfi_ztrtrs(a, b, uplo, trans, diag, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(:,:)
    complex(REAL64), intent(inout) :: b(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: trans
    character :: local_trans
    character, intent(in), optional :: diag
    character :: local_diag
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, nrhs, lda, ldb
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    if (present(trans)) then
        local_trans = trans
    else
        local_trans = 'N'
    end if
    if (present(diag)) then
        local_diag = diag
    else
        local_diag = 'N'
    end if
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    nrhs = size(b,2)
    call ztrtrs(local_uplo, local_trans, local_diag, n, nrhs, a, lda, b, ldb, local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('ztrtrs', local_info)
    end if
end subroutine
!> Modern interface for [[f77_sytrf:ssytrf]].
!> See also: [[mfi_sytrf]], [[f77_sytrf]].
!> Computes the factorization of a symmetric matrix using the Bunch-Kaufman diagonal pivoting method
!> 
!> The factorization has the form:
!> - A = U*D*U**T (if uplo='U') or 
!> - A = L*D*L**T (if uplo='L')
!> 
!> where U (or L) is a product of permutation and unit upper (lower) triangular matrices,
!> and D is block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> Parameters:
!> - `a` (inout): On entry, the symmetric matrix A. On exit, the block diagonal matrix D 
!>                and the multipliers used to obtain the factor U or L.
!> - `uplo` (in, optional): Specifies whether the upper ('U') or lower ('L') triangular part 
!>                of the symmetric matrix A is stored. Default: 'U'
!> - `ipiv` (out, optional): The pivot indices that define the permutation matrix P. 
!>                If ipiv is not provided, it will be allocated internally.
!> - `info` (out, optional): Output status: 0 for success, < 0 for illegal argument, 
!>                > 0 if D(k,k) is exactly zero.
pure subroutine mfi_ssytrf(a, uplo, ipiv, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional, target :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, lwork, allocation_status, deallocation_status
    integer, pointer :: local_ipiv(:)
    real(REAL32), pointer :: work(:)
    real(REAL32) :: s_work(1)  ! Work array for workspace query
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

    ! Retrieve work array size
    lwork = -1
    call ssytrf(local_uplo, n, a, lda, local_ipiv, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call ssytrf(local_uplo, n, a, lda, local_ipiv, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('ssytrf', -local_info)
    end if
    
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        if (local_info <= -1000) then
            call mfi_error('ssytrf', -local_info)
        else
            call mfi_error('ssytrf', local_info)
        end if
    end if
end subroutine
!> Modern interface for [[f77_sytrf:dsytrf]].
!> See also: [[mfi_sytrf]], [[f77_sytrf]].
!> Computes the factorization of a symmetric matrix using the Bunch-Kaufman diagonal pivoting method
!> 
!> The factorization has the form:
!> - A = U*D*U**T (if uplo='U') or 
!> - A = L*D*L**T (if uplo='L')
!> 
!> where U (or L) is a product of permutation and unit upper (lower) triangular matrices,
!> and D is block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!>
!> Parameters:
!> - `a` (inout): On entry, the symmetric matrix A. On exit, the block diagonal matrix D 
!>                and the multipliers used to obtain the factor U or L.
!> - `uplo` (in, optional): Specifies whether the upper ('U') or lower ('L') triangular part 
!>                of the symmetric matrix A is stored. Default: 'U'
!> - `ipiv` (out, optional): The pivot indices that define the permutation matrix P. 
!>                If ipiv is not provided, it will be allocated internally.
!> - `info` (out, optional): Output status: 0 for success, < 0 for illegal argument, 
!>                > 0 if D(k,k) is exactly zero.
pure subroutine mfi_dsytrf(a, uplo, ipiv, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    character, intent(in), optional :: uplo
    character :: local_uplo
    integer, intent(out), optional, target :: ipiv(:)
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: n, lda, lwork, allocation_status, deallocation_status
    integer, pointer :: local_ipiv(:)
    real(REAL64), pointer :: work(:)
    real(REAL64) :: s_work(1)  ! Work array for workspace query
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

    ! Retrieve work array size
    lwork = -1
    call dsytrf(local_uplo, n, a, lda, local_ipiv, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call dsytrf(local_uplo, n, a, lda, local_ipiv, work, lwork, local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('dsytrf', -local_info)
    end if
    
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        if (local_info <= -1000) then
            call mfi_error('dsytrf', -local_info)
        else
            call mfi_error('dsytrf', local_info)
        end if
    end if
end subroutine

    pure subroutine mfi_error(name, info)
        character(*), intent(in) :: name
        integer, intent(in) :: info
        call f77_xerbla(name, info)
    end subroutine

end module
