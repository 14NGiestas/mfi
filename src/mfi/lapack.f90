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
!> Generic modern interface for HEEVR.
!> Supports c, z.
!> See also:
!> [[f77_heevr:cheevr]], [[f77_heevr:zheevr]].
interface mfi_heevr
    module procedure :: mfi_cheevr
    module procedure :: mfi_zheevr
end interface
!> Generic modern interface for HEEVX.
!> Supports c, z.
!> See also:
!> [[f77_heevx:cheevx]], [[f77_heevx:zheevx]].
interface mfi_heevx
    module procedure :: mfi_cheevx
    module procedure :: mfi_zheevx
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
!> Generic modern interface for ORG2R.
!> Supports s, d.
!> See also:
!> [[f77_org2r:sorg2r]], [[f77_org2r:dorg2r]].
interface mfi_org2r
    module procedure :: mfi_sorg2r
    module procedure :: mfi_dorg2r
end interface
!> Generic modern interface for UNG2R.
!> Supports c, z.
!> See also:
!> [[f77_ung2r:cung2r]], [[f77_ung2r:zung2r]].
interface mfi_ung2r
    module procedure :: mfi_cung2r
    module procedure :: mfi_zung2r
end interface
!> Generic modern interface for ORM2R.
!> Supports s, d.
!> See also:
!> [[f77_orm2r:sorm2r]], [[f77_orm2r:dorm2r]].
interface mfi_orm2r
    module procedure :: mfi_sorm2r
    module procedure :: mfi_dorm2r
end interface
!> Generic modern interface for UNM2R.
!> Supports c, z.
!> See also:
!> [[f77_unm2r:cunm2r]], [[f77_unm2r:zunm2r]].
interface mfi_unm2r
    module procedure :: mfi_cunm2r
    module procedure :: mfi_zunm2r
end interface
!> Generic modern interface for ORGR2.
!> Supports s, d.
!> See also:
!> [[f77_orgr2:sorgr2]], [[f77_orgr2:dorgr2]].
interface mfi_orgr2
    module procedure :: mfi_sorgr2
    module procedure :: mfi_dorgr2
end interface
!> Generic modern interface for UNGR2.
!> Supports c, z.
!> See also:
!> [[f77_ungr2:cungr2]], [[f77_ungr2:zungr2]].
interface mfi_ungr2
    module procedure :: mfi_cungr2
    module procedure :: mfi_zungr2
end interface
!> Generic modern interface for ORMR2.
!> Supports s, d.
!> See also:
!> [[f77_ormr2:sormr2]], [[f77_ormr2:dormr2]].
interface mfi_ormr2
    module procedure :: mfi_sormr2
    module procedure :: mfi_dormr2
end interface
!> Generic modern interface for UNMR2.
!> Supports c, z.
!> See also:
!> [[f77_unmr2:cunmr2]], [[f77_unmr2:zunmr2]].
interface mfi_unmr2
    module procedure :: mfi_cunmr2
    module procedure :: mfi_zunmr2
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

!> Modern interface for [[f77_geqrf:f77_geqrf]].
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
pure subroutine mfi_sgeqrf(a, tau, info)
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
    call f77_geqrf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_geqrf(m,n,a,lda,local_tau,work,lwork,local_info)
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
        call mfi_error('f77_geqrf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_geqrf:f77_geqrf]].
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
pure subroutine mfi_dgeqrf(a, tau, info)
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
    call f77_geqrf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_geqrf(m,n,a,lda,local_tau,work,lwork,local_info)
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
        call mfi_error('f77_geqrf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_geqrf:f77_geqrf]].
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
pure subroutine mfi_cgeqrf(a, tau, info)
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
    call f77_geqrf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_geqrf(m,n,a,lda,local_tau,work,lwork,local_info)
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
        call mfi_error('f77_geqrf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_geqrf:f77_geqrf]].
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
pure subroutine mfi_zgeqrf(a, tau, info)
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
    call f77_geqrf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_geqrf(m,n,a,lda,local_tau,work,lwork,local_info)
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
        call mfi_error('f77_geqrf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_gerqf:f77_gerqf]].
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
pure subroutine mfi_sgerqf(a, tau, info)
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
    call f77_gerqf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_gerqf(m,n,a,lda,local_tau,work,lwork,local_info)
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
        call mfi_error('f77_gerqf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_gerqf:f77_gerqf]].
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
pure subroutine mfi_dgerqf(a, tau, info)
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
    call f77_gerqf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_gerqf(m,n,a,lda,local_tau,work,lwork,local_info)
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
        call mfi_error('f77_gerqf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_gerqf:f77_gerqf]].
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
pure subroutine mfi_cgerqf(a, tau, info)
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
    call f77_gerqf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_gerqf(m,n,a,lda,local_tau,work,lwork,local_info)
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
        call mfi_error('f77_gerqf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_gerqf:f77_gerqf]].
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
pure subroutine mfi_zgerqf(a, tau, info)
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
    call f77_gerqf(m,n,a,lda,local_tau,s_work,lwork,local_info)
    if (local_info /= 0) goto 404

    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_gerqf(m,n,a,lda,local_tau,work,lwork,local_info)
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
        call mfi_error('f77_gerqf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_getrf:f77_getrf]].
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
        call f77_getrf(m,n,a,lda,local_ipiv,local_info)
    else
        local_info = -1000
    end if
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getrf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_getrf:f77_getrf]].
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
        call f77_getrf(m,n,a,lda,local_ipiv,local_info)
    else
        local_info = -1000
    end if
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getrf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_getrf:f77_getrf]].
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
        call f77_getrf(m,n,a,lda,local_ipiv,local_info)
    else
        local_info = -1000
    end if
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getrf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_getrf:f77_getrf]].
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
        call f77_getrf(m,n,a,lda,local_ipiv,local_info)
    else
        local_info = -1000
    end if
    if (.not. present(ipiv)) then
        deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getrf', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_getri:f77_getri]].
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
    call f77_getri(n,a,lda,ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_getri(n,a,lda,ipiv,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getri',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_getri:f77_getri]].
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
    call f77_getri(n,a,lda,ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_getri(n,a,lda,ipiv,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getri',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_getri:f77_getri]].
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
    call f77_getri(n,a,lda,ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_getri(n,a,lda,ipiv,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getri',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_getri:f77_getri]].
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
    call f77_getri(n,a,lda,ipiv,s_work,lwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_getri(n,a,lda,ipiv,work,lwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getri',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_getrs:f77_getrs]].
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
    call f77_getrs(local_trans,n,nrhs,a,lda,ipiv,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getrs',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_getrs:f77_getrs]].
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
    call f77_getrs(local_trans,n,nrhs,a,lda,ipiv,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getrs',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_getrs:f77_getrs]].
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
    call f77_getrs(local_trans,n,nrhs,a,lda,ipiv,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getrs',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_getrs:f77_getrs]].
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
    call f77_getrs(local_trans,n,nrhs,a,lda,ipiv,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_getrs',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_hetrf:f77_hetrf]].
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
    complex(REAL32) :: s_work(1)  ! Work array for workspace query
    integer, pointer :: local_ipiv(:)
    complex(REAL32), pointer :: work(:)
    local_ipiv => null()
    work => null()
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
    call f77_hetrf(local_uplo, n, a, lda, local_ipiv, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_hetrf(local_uplo, n, a, lda, local_ipiv, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(ipiv)) then
        if (associated(local_ipiv)) deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_hetrf', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_hetrf:f77_hetrf]].
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
    complex(REAL64) :: s_work(1)  ! Work array for workspace query
    integer, pointer :: local_ipiv(:)
    complex(REAL64), pointer :: work(:)
    local_ipiv => null()
    work => null()
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
    call f77_hetrf(local_uplo, n, a, lda, local_ipiv, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_hetrf(local_uplo, n, a, lda, local_ipiv, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(ipiv)) then
        if (associated(local_ipiv)) deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_hetrf', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_hegv:f77_hegv]].
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
    call f77_hegv(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,s_work,lwork,rwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_hegv(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,work,lwork,rwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    deallocate(rwork, stat=deallocation_status)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_hegv', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_hegv:f77_hegv]].
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
    call f77_hegv(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,s_work,lwork,rwork,local_info)
    if (local_info /= 0) goto 404
    lwork = int(s_work(1))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_hegv(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,work,lwork,rwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    deallocate(rwork, stat=deallocation_status)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_hegv', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_heevd:f77_heevd]].
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

    call f77_heevd(local_jobz,local_uplo,n,a,lda,w, &
                      s_work,lwork,s_rwork,lrwork,s_iwork,liwork,local_info)
    if (local_info /= 0) goto 404
    lwork  = int(s_work(1))
    lrwork = int(s_rwork(1))
    liwork = int(s_iwork(1))

    allocate(iwork(liwork), stat=allocation_status)

    if (allocation_status == 0) then
        allocate(rwork(lrwork), stat=allocation_status)
        allocate(work(lwork),   stat=allocation_status)
        call f77_heevd(local_jobz,local_uplo,n,a,lda,w, &
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
        call mfi_error('f77_heevd', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_heevd:f77_heevd]].
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

    call f77_heevd(local_jobz,local_uplo,n,a,lda,w, &
                      s_work,lwork,s_rwork,lrwork,s_iwork,liwork,local_info)
    if (local_info /= 0) goto 404
    lwork  = int(s_work(1))
    lrwork = int(s_rwork(1))
    liwork = int(s_iwork(1))

    allocate(iwork(liwork), stat=allocation_status)

    if (allocation_status == 0) then
        allocate(rwork(lrwork), stat=allocation_status)
        allocate(work(lwork),   stat=allocation_status)
        call f77_heevd(local_jobz,local_uplo,n,a,lda,w, &
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
        call mfi_error('f77_heevd', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_heevr:f77_heevr]].
!> See also: [[mfi_heevr]], [[f77_heevr]].
!> Computes selected eigenvalues and, optionally, eigenvectors of a Hermitian matrix using RRR
pure subroutine mfi_cheevr(a, w, jobz, uplo, range, vl, vu, il, iu, abstol, m, z, isuppz, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(out) :: w(:)
    complex(REAL32),      intent(out), optional, target :: z(:,:)
    integer,           intent(out), optional, target :: isuppz(:)
    integer,           intent(out), optional :: m
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: range
    character :: local_range
    real(REAL32), intent(in), optional :: vl
    real(REAL32) :: local_vl
    real(REAL32), intent(in), optional :: vu
    real(REAL32) :: local_vu
    real(REAL32), intent(in), optional :: abstol
    real(REAL32) :: local_abstol
    integer, intent(in), optional :: il
    integer :: local_il
    integer, intent(in), optional :: iu
    integer :: local_iu
    integer, intent(out), optional :: info
    integer :: local_info
    complex(REAL32),      pointer :: work(:), local_z(:,:)
    real(REAL32),      pointer :: rwork(:)
    integer,           pointer :: iwork(:), local_isuppz(:)
    complex(REAL32)      :: s_work(1)
    real(REAL32) :: s_rwork(1)
    integer       :: s_iwork(1)
    integer :: n, lda, ldz, lwork, lrwork, liwork, allocation_status, deallocation_status
    integer :: local_m
    complex(REAL32), target :: l_z(1,1)
    integer, target :: l_isuppz(2*size(a,1))  ! Default size for isuppz when not present
    work => null()
    rwork => null()
    iwork => null()
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
    if (present(range)) then
        local_range = range
    else
        local_range = 'A'
    end if
    if (present(vl)) then
        local_vl = vl
    else
        local_vl = 0.0_wp
    end if
    if (present(vu)) then
        local_vu = vu
    else
        local_vu = 0.0_wp
    end if
    if (present(il)) then
        local_il = il
    else
        local_il = 1
    end if
    if (present(iu)) then
        local_iu = iu
    else
        local_iu = size(a,1)
    end if
    if (present(abstol)) then
        local_abstol = abstol
    else
        local_abstol = 0.0_wp
    end if
    lda = max(1,size(a,1))
    n   = size(a,2)
    ldz = max(1,size(a,1))  ! Size for eigenvector matrix
    if (present(z)) ldz = max(1,size(z,1))
    allocation_status = 0
    lwork  = -1
    lrwork = -1
    liwork = -1
    
    ! Set up pointers
    if (present(z)) then
        local_z => z
    else
        local_z => l_z
    end if
    
    if (present(isuppz)) then
        local_isuppz => isuppz
    else
        local_isuppz => l_isuppz
    end if

    ! Query workspace size
    call f77_heevr(local_jobz, local_range, local_uplo, n, a, lda, &
                      local_vl, local_vu, local_il, local_iu, local_abstol, &
                      local_m, w, local_z, ldz, local_isuppz, s_work, lwork, s_rwork, lrwork, &
                      s_iwork, liwork, local_info)
    if (local_info /= 0) goto 404
    
    lwork  = max(1, int(real(s_work(1), wp)))
    lrwork = max(1, int(s_rwork(1)))
    liwork = max(1, int(s_iwork(1)))

    allocate(iwork(liwork), stat=allocation_status)
    if (allocation_status == 0) then
        allocate(rwork(lrwork), stat=allocation_status)
        if (allocation_status == 0) then
            allocate(work(lwork), stat=allocation_status)
            if (allocation_status == 0) then
                ! Call the actual routine
                call f77_heevr(local_jobz, local_range, local_uplo, n, a, lda, &
                              local_vl, local_vu, local_il, local_iu, local_abstol, &
                              local_m, w, local_z, ldz, local_isuppz, work, lwork, rwork, lrwork, &
                              iwork, liwork, local_info)

                ! Copy results if needed
                if (present(m)) m = local_m
            else
                local_info = -1000
            end if
        else
            local_info = -1000
        end if
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)
    if (associated(rwork)) deallocate(rwork, stat=deallocation_status)
    if (associated(iwork)) deallocate(iwork, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_heevr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_heevr:f77_heevr]].
!> See also: [[mfi_heevr]], [[f77_heevr]].
!> Computes selected eigenvalues and, optionally, eigenvectors of a Hermitian matrix using RRR
pure subroutine mfi_zheevr(a, w, jobz, uplo, range, vl, vu, il, iu, abstol, m, z, isuppz, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(out) :: w(:)
    complex(REAL64),      intent(out), optional, target :: z(:,:)
    integer,           intent(out), optional, target :: isuppz(:)
    integer,           intent(out), optional :: m
    character, intent(in), optional :: jobz
    character :: local_jobz
    character, intent(in), optional :: uplo
    character :: local_uplo
    character, intent(in), optional :: range
    character :: local_range
    real(REAL64), intent(in), optional :: vl
    real(REAL64) :: local_vl
    real(REAL64), intent(in), optional :: vu
    real(REAL64) :: local_vu
    real(REAL64), intent(in), optional :: abstol
    real(REAL64) :: local_abstol
    integer, intent(in), optional :: il
    integer :: local_il
    integer, intent(in), optional :: iu
    integer :: local_iu
    integer, intent(out), optional :: info
    integer :: local_info
    complex(REAL64),      pointer :: work(:), local_z(:,:)
    real(REAL64),      pointer :: rwork(:)
    integer,           pointer :: iwork(:), local_isuppz(:)
    complex(REAL64)      :: s_work(1)
    real(REAL64) :: s_rwork(1)
    integer       :: s_iwork(1)
    integer :: n, lda, ldz, lwork, lrwork, liwork, allocation_status, deallocation_status
    integer :: local_m
    complex(REAL64), target :: l_z(1,1)
    integer, target :: l_isuppz(2*size(a,1))  ! Default size for isuppz when not present
    work => null()
    rwork => null()
    iwork => null()
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
    if (present(range)) then
        local_range = range
    else
        local_range = 'A'
    end if
    if (present(vl)) then
        local_vl = vl
    else
        local_vl = 0.0_wp
    end if
    if (present(vu)) then
        local_vu = vu
    else
        local_vu = 0.0_wp
    end if
    if (present(il)) then
        local_il = il
    else
        local_il = 1
    end if
    if (present(iu)) then
        local_iu = iu
    else
        local_iu = size(a,1)
    end if
    if (present(abstol)) then
        local_abstol = abstol
    else
        local_abstol = 0.0_wp
    end if
    lda = max(1,size(a,1))
    n   = size(a,2)
    ldz = max(1,size(a,1))  ! Size for eigenvector matrix
    if (present(z)) ldz = max(1,size(z,1))
    allocation_status = 0
    lwork  = -1
    lrwork = -1
    liwork = -1
    
    ! Set up pointers
    if (present(z)) then
        local_z => z
    else
        local_z => l_z
    end if
    
    if (present(isuppz)) then
        local_isuppz => isuppz
    else
        local_isuppz => l_isuppz
    end if

    ! Query workspace size
    call f77_heevr(local_jobz, local_range, local_uplo, n, a, lda, &
                      local_vl, local_vu, local_il, local_iu, local_abstol, &
                      local_m, w, local_z, ldz, local_isuppz, s_work, lwork, s_rwork, lrwork, &
                      s_iwork, liwork, local_info)
    if (local_info /= 0) goto 404
    
    lwork  = max(1, int(real(s_work(1), wp)))
    lrwork = max(1, int(s_rwork(1)))
    liwork = max(1, int(s_iwork(1)))

    allocate(iwork(liwork), stat=allocation_status)
    if (allocation_status == 0) then
        allocate(rwork(lrwork), stat=allocation_status)
        if (allocation_status == 0) then
            allocate(work(lwork), stat=allocation_status)
            if (allocation_status == 0) then
                ! Call the actual routine
                call f77_heevr(local_jobz, local_range, local_uplo, n, a, lda, &
                              local_vl, local_vu, local_il, local_iu, local_abstol, &
                              local_m, w, local_z, ldz, local_isuppz, work, lwork, rwork, lrwork, &
                              iwork, liwork, local_info)

                ! Copy results if needed
                if (present(m)) m = local_m
            else
                local_info = -1000
            end if
        else
            local_info = -1000
        end if
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)
    if (associated(rwork)) deallocate(rwork, stat=deallocation_status)
    if (associated(iwork)) deallocate(iwork, stat=deallocation_status)
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_heevr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_heevx:f77_heevx]].
!> See also: [[mfi_heevx]], [[f77_heevx]].
!> Computes selected eigenvalues and, optionally, eigenvectors
!> of a complex Hermitian matrix.

pure subroutine mfi_cheevx(a, w, uplo, z, vl, vu, il, iu, m, ifail, abstol, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(out) :: w(:)
    character(len=1), intent(in), optional :: uplo
    complex(REAL32),      intent(out), optional :: z(:,:)
    real(REAL32),      intent(in), optional :: vl, vu, abstol
    integer,          intent(in), optional :: il, iu
    integer,          intent(out), optional :: m, info
    integer, intent(out), optional, target :: ifail(:)
    complex(REAL32),      pointer :: work(:)
    real(REAL32), pointer :: rwork(:)
    integer,       pointer :: iwork(:), local_ifail(:)
    complex(REAL32)      :: s_work(1)
    real(REAL32) :: s_rwork(1)
    integer       :: s_iwork(1)
    ! Define a dummy array to use when needed (declarations first)
    integer, target :: dummy_ifail(1)
    integer :: n, lda, ldz, lwork, allocation_status, deallocation_status
    logical :: owned_ifail
    character(1) :: jobz, range
    character(len=1) :: local_uplo
    real(REAL32) :: local_vl, local_vu, local_abstol
    integer :: local_il, local_iu
    integer :: local_m, local_info

    ! Initialize pointers to null to prevent deallocating unallocated pointers
    work => null()
    rwork => null()
    iwork => null()
    local_ifail => null()
    owned_ifail = .false.

    ! Set defaults
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    n = size(a, 1)
    lda = max(1, size(a, 1))

    ! Determine JOBZ based on presence of Z
    if (present(z)) then
        jobz = 'V'
        ldz = max(1, size(z, 1))
    else
        jobz = 'N'
        ldz = 1
    end if

    ! Determine RANGE based on which optional parameters are present
    ! Following reference implementation logic exactly:
    IF((PRESENT(vl).OR.PRESENT(vu)).AND.(PRESENT(il).OR.PRESENT(iu))) THEN
        local_info=-1001; GOTO 404
    ELSEIF((PRESENT(vl).OR.PRESENT(vu))) THEN
        range = 'V'
    ELSEIF((PRESENT(il).OR.PRESENT(iu))) THEN
        range = 'I'
    ELSE
        range = 'A'
    ENDIF

    ! Handle ifail allocation - following reference implementation logic exactly
    IF(.NOT.PRESENT(z)) THEN
        IF(PRESENT(ifail)) THEN
            local_info=-1001; GOTO 404  ! Error: IFAIL provided but Z not present
        ELSE
            ! Point to dummy array when Z not present and IFAIL not present
            local_ifail => dummy_ifail
        ENDIF
    ELSE
        ! Z is present
        IF(PRESENT(ifail)) THEN
            local_ifail => ifail
        ELSE
            ! Allocate IFAIL when Z is present but IFAIL not provided
            allocate(local_ifail(n), stat=allocation_status)
            if (allocation_status /= 0) then
                local_info = -1000
                goto 404
            end if
            owned_ifail = .true.
        ENDIF
    ENDIF

    ! Set defaults following reference implementation pattern exactly
    ! Reference: IF(PRESENT(IL)) THEN O_IL = IL ELSE O_IL = 1 ENDIF
    if (present(il)) then
        local_il = il
    else
        local_il = 1
    end if
    if (present(iu)) then
        local_iu = iu
    else
        local_iu = n
    end if
    if (present(vl)) then
        local_vl = vl
    else
        local_vl = -huge(0.0_wp)
    end if
    if (present(vu)) then
        local_vu = vu
    else
        local_vu = huge(0.0_wp)
    end if
    if (present(abstol)) then
        local_abstol = abstol
    else
        local_abstol = 0.0_wp
    end if

    allocation_status = 0

    ! Query workspace sizes
    lwork = -1

    ! For workspace query, use appropriate arrays based on presence
    if (present(z)) then
        call f77_heevx(jobz, range, local_uplo, n, a, lda, local_vl, local_vu, &
                          local_il, local_iu, local_abstol, local_m, w, z, ldz, &
                          s_work, lwork, s_rwork, s_iwork, local_ifail, local_info)
    else
        ! When z is not present, use a dummy array with correct size
        ! Using a 1x1 array is safe since jobz='N' means z is not accessed
        call f77_heevx(jobz, range, local_uplo, n, a, lda, local_vl, local_vu, &
                          local_il, local_iu, local_abstol, local_m, w, a, ldz, &
                          s_work, lwork, s_rwork, s_iwork, local_ifail, local_info)
    end if

    if (local_info /= 0) goto 404

    ! Get the optimal workspace sizes from the query
    lwork = max(1, int(real(s_work(1), wp)))

    ! Allocate workspace arrays
    allocate(work(lwork),   stat=allocation_status)
    if (allocation_status == 0) then
        allocate(rwork(max(1, int(s_rwork(1)))), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        allocate(iwork(max(1, int(s_iwork(1)))), stat=allocation_status)
    end if

    if (allocation_status == 0) then
        ! Call main routine with actual workspace
        if (present(z)) then
            call f77_heevx(jobz, range, local_uplo, n, a, lda, local_vl, local_vu, &
                              local_il, local_iu, local_abstol, local_m, w, z, ldz, &
                              work, lwork, rwork, iwork, local_ifail, local_info)
        else
            ! When z not present, pass a as dummy array (safe when jobz='N')
            call f77_heevx(jobz, range, local_uplo, n, a, lda, local_vl, local_vu, &
                              local_il, local_iu, local_abstol, local_m, w, a, ldz, &
                              work, lwork, rwork, iwork, local_ifail, local_info)
        end if
    else
        local_info = -1000
    end if

    ! Set optional output parameters
    if (present(m)) m = local_m

    ! Deallocate if local arrays were allocated
    ! Only deallocate if we own the allocation (pointer is associated)
    if (associated(work)) deallocate(work, stat=deallocation_status)
    if (associated(rwork)) deallocate(rwork, stat=deallocation_status)
    if (associated(iwork)) deallocate(iwork, stat=deallocation_status)
    if (owned_ifail) deallocate(local_ifail, stat=deallocation_status)

404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_heevx', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_heevx:f77_heevx]].
!> See also: [[mfi_heevx]], [[f77_heevx]].
!> Computes selected eigenvalues and, optionally, eigenvectors
!> of a complex Hermitian matrix.

pure subroutine mfi_zheevx(a, w, uplo, z, vl, vu, il, iu, m, ifail, abstol, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(out) :: w(:)
    character(len=1), intent(in), optional :: uplo
    complex(REAL64),      intent(out), optional :: z(:,:)
    real(REAL64),      intent(in), optional :: vl, vu, abstol
    integer,          intent(in), optional :: il, iu
    integer,          intent(out), optional :: m, info
    integer, intent(out), optional, target :: ifail(:)
    complex(REAL64),      pointer :: work(:)
    real(REAL64), pointer :: rwork(:)
    integer,       pointer :: iwork(:), local_ifail(:)
    complex(REAL64)      :: s_work(1)
    real(REAL64) :: s_rwork(1)
    integer       :: s_iwork(1)
    ! Define a dummy array to use when needed (declarations first)
    integer, target :: dummy_ifail(1)
    integer :: n, lda, ldz, lwork, allocation_status, deallocation_status
    logical :: owned_ifail
    character(1) :: jobz, range
    character(len=1) :: local_uplo
    real(REAL64) :: local_vl, local_vu, local_abstol
    integer :: local_il, local_iu
    integer :: local_m, local_info

    ! Initialize pointers to null to prevent deallocating unallocated pointers
    work => null()
    rwork => null()
    iwork => null()
    local_ifail => null()
    owned_ifail = .false.

    ! Set defaults
    if (present(uplo)) then
        local_uplo = uplo
    else
        local_uplo = 'U'
    end if
    n = size(a, 1)
    lda = max(1, size(a, 1))

    ! Determine JOBZ based on presence of Z
    if (present(z)) then
        jobz = 'V'
        ldz = max(1, size(z, 1))
    else
        jobz = 'N'
        ldz = 1
    end if

    ! Determine RANGE based on which optional parameters are present
    ! Following reference implementation logic exactly:
    IF((PRESENT(vl).OR.PRESENT(vu)).AND.(PRESENT(il).OR.PRESENT(iu))) THEN
        local_info=-1001; GOTO 404
    ELSEIF((PRESENT(vl).OR.PRESENT(vu))) THEN
        range = 'V'
    ELSEIF((PRESENT(il).OR.PRESENT(iu))) THEN
        range = 'I'
    ELSE
        range = 'A'
    ENDIF

    ! Handle ifail allocation - following reference implementation logic exactly
    IF(.NOT.PRESENT(z)) THEN
        IF(PRESENT(ifail)) THEN
            local_info=-1001; GOTO 404  ! Error: IFAIL provided but Z not present
        ELSE
            ! Point to dummy array when Z not present and IFAIL not present
            local_ifail => dummy_ifail
        ENDIF
    ELSE
        ! Z is present
        IF(PRESENT(ifail)) THEN
            local_ifail => ifail
        ELSE
            ! Allocate IFAIL when Z is present but IFAIL not provided
            allocate(local_ifail(n), stat=allocation_status)
            if (allocation_status /= 0) then
                local_info = -1000
                goto 404
            end if
            owned_ifail = .true.
        ENDIF
    ENDIF

    ! Set defaults following reference implementation pattern exactly
    ! Reference: IF(PRESENT(IL)) THEN O_IL = IL ELSE O_IL = 1 ENDIF
    if (present(il)) then
        local_il = il
    else
        local_il = 1
    end if
    if (present(iu)) then
        local_iu = iu
    else
        local_iu = n
    end if
    if (present(vl)) then
        local_vl = vl
    else
        local_vl = -huge(0.0_wp)
    end if
    if (present(vu)) then
        local_vu = vu
    else
        local_vu = huge(0.0_wp)
    end if
    if (present(abstol)) then
        local_abstol = abstol
    else
        local_abstol = 0.0_wp
    end if

    allocation_status = 0

    ! Query workspace sizes
    lwork = -1

    ! For workspace query, use appropriate arrays based on presence
    if (present(z)) then
        call f77_heevx(jobz, range, local_uplo, n, a, lda, local_vl, local_vu, &
                          local_il, local_iu, local_abstol, local_m, w, z, ldz, &
                          s_work, lwork, s_rwork, s_iwork, local_ifail, local_info)
    else
        ! When z is not present, use a dummy array with correct size
        ! Using a 1x1 array is safe since jobz='N' means z is not accessed
        call f77_heevx(jobz, range, local_uplo, n, a, lda, local_vl, local_vu, &
                          local_il, local_iu, local_abstol, local_m, w, a, ldz, &
                          s_work, lwork, s_rwork, s_iwork, local_ifail, local_info)
    end if

    if (local_info /= 0) goto 404

    ! Get the optimal workspace sizes from the query
    lwork = max(1, int(real(s_work(1), wp)))

    ! Allocate workspace arrays
    allocate(work(lwork),   stat=allocation_status)
    if (allocation_status == 0) then
        allocate(rwork(max(1, int(s_rwork(1)))), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        allocate(iwork(max(1, int(s_iwork(1)))), stat=allocation_status)
    end if

    if (allocation_status == 0) then
        ! Call main routine with actual workspace
        if (present(z)) then
            call f77_heevx(jobz, range, local_uplo, n, a, lda, local_vl, local_vu, &
                              local_il, local_iu, local_abstol, local_m, w, z, ldz, &
                              work, lwork, rwork, iwork, local_ifail, local_info)
        else
            ! When z not present, pass a as dummy array (safe when jobz='N')
            call f77_heevx(jobz, range, local_uplo, n, a, lda, local_vl, local_vu, &
                              local_il, local_iu, local_abstol, local_m, w, a, ldz, &
                              work, lwork, rwork, iwork, local_ifail, local_info)
        end if
    else
        local_info = -1000
    end if

    ! Set optional output parameters
    if (present(m)) m = local_m

    ! Deallocate if local arrays were allocated
    ! Only deallocate if we own the allocation (pointer is associated)
    if (associated(work)) deallocate(work, stat=deallocation_status)
    if (associated(rwork)) deallocate(rwork, stat=deallocation_status)
    if (associated(iwork)) deallocate(iwork, stat=deallocation_status)
    if (owned_ifail) deallocate(local_ifail, stat=deallocation_status)

404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_heevx', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_gesvd:f77_gesvd]].
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
    call f77_gesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,local_info)
    if (local_info /= 0) then
        goto 404
    end if
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_gesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,local_info)
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
        call mfi_error('f77_gesvd', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_gesvd:f77_gesvd]].
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
    call f77_gesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,local_info)
    if (local_info /= 0) then
        goto 404
    end if
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_gesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,local_info)
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
        call mfi_error('f77_gesvd', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_gesvd:f77_gesvd]].
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
    call f77_gesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,rwork,local_info)
    if (local_info /= 0) then
        goto 404
    end if
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_gesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,rwork,local_info)
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
        call mfi_error('f77_gesvd', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_gesvd:f77_gesvd]].
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
    call f77_gesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,rwork,local_info)
    if (local_info /= 0) then
        goto 404
    end if
    lwork = int(s_work(1))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_gesvd(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,rwork,local_info)
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
        call mfi_error('f77_gesvd', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_orgqr:f77_orgqr]].
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
    real(REAL32), target  :: s_work(1)
    real(REAL32), pointer :: work(:)
    work => null()
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call f77_orgqr(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_orgqr(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_orgqr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_orgqr:f77_orgqr]].
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
    real(REAL64), target  :: s_work(1)
    real(REAL64), pointer :: work(:)
    work => null()
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call f77_orgqr(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_orgqr(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_orgqr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_orgrq:f77_orgrq]].
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
    real(REAL32), target  :: s_work(1)
    real(REAL32), pointer :: work(:)
    work => null()
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call f77_orgrq(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_orgrq(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_orgrq', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_orgrq:f77_orgrq]].
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
    real(REAL64), target  :: s_work(1)
    real(REAL64), pointer :: work(:)
    work => null()
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call f77_orgrq(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_orgrq(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_orgrq', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_ungqr:f77_ungqr]].
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
    complex(REAL32), target  :: s_work(1)
    complex(REAL32), pointer :: work(:)
    work => null()
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call f77_ungqr(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_ungqr(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ungqr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_ungqr:f77_ungqr]].
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
    complex(REAL64), target  :: s_work(1)
    complex(REAL64), pointer :: work(:)
    work => null()
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call f77_ungqr(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_ungqr(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ungqr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_ungrq:f77_ungrq]].
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
    complex(REAL32), target  :: s_work(1)
    complex(REAL32), pointer :: work(:)
    work => null()
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call f77_ungrq(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_ungrq(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ungrq', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_ungrq:f77_ungrq]].
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
    complex(REAL64), target  :: s_work(1)
    complex(REAL64), pointer :: work(:)
    work => null()
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)

    ! Query workspace size
    lwork = -1
    call f77_ungrq(m, n, local_k, a, lda, tau, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_ungrq(m, n, local_k, a, lda, tau, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ungrq', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_ormqr:f77_ormqr]].
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
    real(REAL32), target  :: s_work(1)
    real(REAL32), pointer :: work(:)
    work => null()
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
    k = size(tau)

    ! Query workspace size
    lwork = -1
    call f77_ormqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_ormqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ormqr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_ormqr:f77_ormqr]].
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
    real(REAL64), target  :: s_work(1)
    real(REAL64), pointer :: work(:)
    work => null()
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
    k = size(tau)

    ! Query workspace size
    lwork = -1
    call f77_ormqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_ormqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ormqr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_ormrq:f77_ormrq]].
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
    real(REAL32), target  :: s_work(1)
    real(REAL32), pointer :: work(:)
    work => null()
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
    k = size(tau)

    ! Query workspace size
    lwork = -1
    call f77_ormrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_ormrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ormrq', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_ormrq:f77_ormrq]].
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
    real(REAL64), target  :: s_work(1)
    real(REAL64), pointer :: work(:)
    work => null()
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
    k = size(tau)

    ! Query workspace size
    lwork = -1
    call f77_ormrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_ormrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ormrq', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_unmqr:f77_unmqr]].
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
    complex(REAL32), target  :: s_work(1)
    complex(REAL32), pointer :: work(:)
    work => null()
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
    k = size(tau)

    ! Query workspace size
    lwork = -1
    call f77_unmqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_unmqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_unmqr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_unmqr:f77_unmqr]].
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
    complex(REAL64), target  :: s_work(1)
    complex(REAL64), pointer :: work(:)
    work => null()
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
    k = size(tau)

    ! Query workspace size
    lwork = -1
    call f77_unmqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_unmqr(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_unmqr', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_unmrq:f77_unmrq]].
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
    complex(REAL32), target  :: s_work(1)
    complex(REAL32), pointer :: work(:)
    work => null()
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
    k = size(tau)

    ! Query workspace size
    lwork = -1
    call f77_unmrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_unmrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_unmrq', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_unmrq:f77_unmrq]].
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
    complex(REAL64), target  :: s_work(1)
    complex(REAL64), pointer :: work(:)
    work => null()
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
    k = size(tau)

    ! Query workspace size
    lwork = -1
    call f77_unmrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(real(s_work(1), wp)))
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
        call f77_unmrq(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_unmrq', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_org2r:f77_org2r]].
!> See also: [[mfi_org2r]], [[f77_org2r]].
!> Generates the real orthogonal matrix Q of the QR factorization formed by geqr2
pure subroutine mfi_sorg2r(a, tau, k, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda
    real(REAL32), allocatable :: work(:)
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau,1)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)
    local_k = min(local_k, min(m, n))

    allocate(work(m)) ! Allocate workspace for generation
    call f77_org2r(m, n, local_k, a, lda, tau, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_org2r', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_org2r:f77_org2r]].
!> See also: [[mfi_org2r]], [[f77_org2r]].
!> Generates the real orthogonal matrix Q of the QR factorization formed by geqr2
pure subroutine mfi_dorg2r(a, tau, k, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda
    real(REAL64), allocatable :: work(:)
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau,1)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)
    local_k = min(local_k, min(m, n))

    allocate(work(m)) ! Allocate workspace for generation
    call f77_org2r(m, n, local_k, a, lda, tau, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_org2r', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ung2r:f77_ung2r]].
!> See also: [[mfi_ung2r]], [[f77_ung2r]].
!> Generates the real orthogonal matrix Q of the QR factorization formed by geqr2
pure subroutine mfi_cung2r(a, tau, k, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda
    complex(REAL32), allocatable :: work(:)
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau,1)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)
    local_k = min(local_k, min(m, n))

    allocate(work(m)) ! Allocate workspace for generation
    call f77_ung2r(m, n, local_k, a, lda, tau, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ung2r', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ung2r:f77_ung2r]].
!> See also: [[mfi_ung2r]], [[f77_ung2r]].
!> Generates the real orthogonal matrix Q of the QR factorization formed by geqr2
pure subroutine mfi_zung2r(a, tau, k, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda
    complex(REAL64), allocatable :: work(:)
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau,1)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)
    local_k = min(local_k, min(m, n))

    allocate(work(m)) ! Allocate workspace for generation
    call f77_ung2r(m, n, local_k, a, lda, tau, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ung2r', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_orm2r:f77_orm2r]].
!> See also: [[mfi_orm2r]], [[f77_orm2r]].
!> Multiplies a real matrix by the orthogonal matrix Q formed by geqr2
pure subroutine mfi_sorm2r(a, tau, c, side, trans, info)
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
    integer :: m, n, k, lda, ldc
    real(REAL32), allocatable :: work(:)
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
        k = size(tau, 1)
    else
        k = size(tau, 1)
    end if

    allocate(work(m*n)) ! Allocate sufficient workspace
    call f77_orm2r(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_orm2r', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_orm2r:f77_orm2r]].
!> See also: [[mfi_orm2r]], [[f77_orm2r]].
!> Multiplies a real matrix by the orthogonal matrix Q formed by geqr2
pure subroutine mfi_dorm2r(a, tau, c, side, trans, info)
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
    integer :: m, n, k, lda, ldc
    real(REAL64), allocatable :: work(:)
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
        k = size(tau, 1)
    else
        k = size(tau, 1)
    end if

    allocate(work(m*n)) ! Allocate sufficient workspace
    call f77_orm2r(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_orm2r', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_unm2r:f77_unm2r]].
!> See also: [[mfi_unm2r]], [[f77_unm2r]].
!> Multiplies a real matrix by the orthogonal matrix Q formed by geqr2
pure subroutine mfi_cunm2r(a, tau, c, side, trans, info)
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
    integer :: m, n, k, lda, ldc
    complex(REAL32), allocatable :: work(:)
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
        k = size(tau, 1)
    else
        k = size(tau, 1)
    end if

    allocate(work(m*n)) ! Allocate sufficient workspace
    call f77_unm2r(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_unm2r', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_unm2r:f77_unm2r]].
!> See also: [[mfi_unm2r]], [[f77_unm2r]].
!> Multiplies a real matrix by the orthogonal matrix Q formed by geqr2
pure subroutine mfi_zunm2r(a, tau, c, side, trans, info)
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
    integer :: m, n, k, lda, ldc
    complex(REAL64), allocatable :: work(:)
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
        k = size(tau, 1)
    else
        k = size(tau, 1)
    end if

    allocate(work(m*n)) ! Allocate sufficient workspace
    call f77_unm2r(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_unm2r', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_orgr2:f77_orgr2]].
!> See also: [[mfi_orgr2]], [[f77_orgr2]].
!> Generates the real orthogonal matrix Q of the RQ factorization formed by gerq2
pure subroutine mfi_sorgr2(a, tau, k, info)
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(:,:)
    real(REAL32), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda
    real(REAL32), allocatable :: work(:)
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau,1)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)
    local_k = min(local_k, min(m, n))

    allocate(work(m)) ! Allocate workspace for generation
    call f77_orgr2(m, n, local_k, a, lda, tau, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_orgr2', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_orgr2:f77_orgr2]].
!> See also: [[mfi_orgr2]], [[f77_orgr2]].
!> Generates the real orthogonal matrix Q of the RQ factorization formed by gerq2
pure subroutine mfi_dorgr2(a, tau, k, info)
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(:,:)
    real(REAL64), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda
    real(REAL64), allocatable :: work(:)
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau,1)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)
    local_k = min(local_k, min(m, n))

    allocate(work(m)) ! Allocate workspace for generation
    call f77_orgr2(m, n, local_k, a, lda, tau, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_orgr2', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ungr2:f77_ungr2]].
!> See also: [[mfi_ungr2]], [[f77_ungr2]].
!> Generates the real orthogonal matrix Q of the RQ factorization formed by gerq2
pure subroutine mfi_cungr2(a, tau, k, info)
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(:,:)
    complex(REAL32), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda
    complex(REAL32), allocatable :: work(:)
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau,1)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)
    local_k = min(local_k, min(m, n))

    allocate(work(m)) ! Allocate workspace for generation
    call f77_ungr2(m, n, local_k, a, lda, tau, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ungr2', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ungr2:f77_ungr2]].
!> See also: [[mfi_ungr2]], [[f77_ungr2]].
!> Generates the real orthogonal matrix Q of the RQ factorization formed by gerq2
pure subroutine mfi_zungr2(a, tau, k, info)
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(:,:)
    complex(REAL64), intent(in) :: tau(:)
    integer, intent(in), optional :: k
    integer :: local_k
    integer, intent(out), optional :: info
    integer :: local_info
    integer :: m, n, lda
    complex(REAL64), allocatable :: work(:)
    if (present(k)) then
        local_k = k
    else
        local_k = size(tau,1)
    end if
    lda = max(1, size(a,1))
    m = size(a, 1)
    n = size(a, 2)
    local_k = min(local_k, min(m, n))

    allocate(work(m)) ! Allocate workspace for generation
    call f77_ungr2(m, n, local_k, a, lda, tau, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ungr2', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ormr2:f77_ormr2]].
!> See also: [[mfi_ormr2]], [[f77_ormr2]].
!> Multiplies a real matrix by the orthogonal matrix Q formed by gerq2
pure subroutine mfi_sormr2(a, tau, c, side, trans, info)
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
    integer :: m, n, k, lda, ldc
    real(REAL32), allocatable :: work(:)
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
        k = size(tau, 1)
    else
        k = size(tau, 1)
    end if

    allocate(work(m*n)) ! Allocate sufficient workspace
    call f77_ormr2(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ormr2', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_ormr2:f77_ormr2]].
!> See also: [[mfi_ormr2]], [[f77_ormr2]].
!> Multiplies a real matrix by the orthogonal matrix Q formed by gerq2
pure subroutine mfi_dormr2(a, tau, c, side, trans, info)
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
    integer :: m, n, k, lda, ldc
    real(REAL64), allocatable :: work(:)
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
        k = size(tau, 1)
    else
        k = size(tau, 1)
    end if

    allocate(work(m*n)) ! Allocate sufficient workspace
    call f77_ormr2(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_ormr2', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_unmr2:f77_unmr2]].
!> See also: [[mfi_unmr2]], [[f77_unmr2]].
!> Multiplies a real matrix by the orthogonal matrix Q formed by gerq2
pure subroutine mfi_cunmr2(a, tau, c, side, trans, info)
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
    integer :: m, n, k, lda, ldc
    complex(REAL32), allocatable :: work(:)
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
        k = size(tau, 1)
    else
        k = size(tau, 1)
    end if

    allocate(work(m*n)) ! Allocate sufficient workspace
    call f77_unmr2(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_unmr2', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_unmr2:f77_unmr2]].
!> See also: [[mfi_unmr2]], [[f77_unmr2]].
!> Multiplies a real matrix by the orthogonal matrix Q formed by gerq2
pure subroutine mfi_zunmr2(a, tau, c, side, trans, info)
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
    integer :: m, n, k, lda, ldc
    complex(REAL64), allocatable :: work(:)
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
        k = size(tau, 1)
    else
        k = size(tau, 1)
    end if

    allocate(work(m*n)) ! Allocate sufficient workspace
    call f77_unmr2(local_side, local_trans, m, n, k, a, lda, tau, c, ldc, work, local_info)
    deallocate(work)

    ! Error handling
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_unmr2', -local_info)
    end if
end subroutine
!> Modern interface for [[f77_potrf:f77_potrf]].
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
    call f77_potrf(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_potrf', local_info)
    end if
end subroutine
!> Modern interface for [[f77_potrf:f77_potrf]].
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
    call f77_potrf(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_potrf', local_info)
    end if
end subroutine
!> Modern interface for [[f77_potrf:f77_potrf]].
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
    call f77_potrf(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_potrf', local_info)
    end if
end subroutine
!> Modern interface for [[f77_potrf:f77_potrf]].
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
    call f77_potrf(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_potrf', local_info)
    end if
end subroutine
!> Modern interface for [[f77_potri:f77_potri]].
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
    call f77_potri(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_potri', local_info)
    end if
end subroutine
!> Modern interface for [[f77_potri:f77_potri]].
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
    call f77_potri(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_potri', local_info)
    end if
end subroutine
!> Modern interface for [[f77_potri:f77_potri]].
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
    call f77_potri(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_potri', local_info)
    end if
end subroutine
!> Modern interface for [[f77_potri:f77_potri]].
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
    call f77_potri(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_potri', local_info)
    end if
end subroutine
!> Modern interface for [[f77_potrs:f77_potrs]].
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
    call f77_potrs(local_uplo,n,nrhs,a,lda,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_potrs',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_potrs:f77_potrs]].
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
    call f77_potrs(local_uplo,n,nrhs,a,lda,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_potrs',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_potrs:f77_potrs]].
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
    call f77_potrs(local_uplo,n,nrhs,a,lda,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_potrs',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_potrs:f77_potrs]].
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
    call f77_potrs(local_uplo,n,nrhs,a,lda,b,ldb,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_potrs',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_pocon:f77_pocon]].
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
        call f77_pocon(local_uplo, n, a, lda, anorm, rcond, work, xwork, local_info)
    else
        local_info = -1000
    end if

    deallocate(xwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)

    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_pocon',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_pocon:f77_pocon]].
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
        call f77_pocon(local_uplo, n, a, lda, anorm, rcond, work, xwork, local_info)
    else
        local_info = -1000
    end if

    deallocate(xwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)

    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_pocon',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_pocon:f77_pocon]].
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
        call f77_pocon(local_uplo, n, a, lda, anorm, rcond, work, xwork, local_info)
    else
        local_info = -1000
    end if

    deallocate(xwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)

    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_pocon',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_pocon:f77_pocon]].
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
        call f77_pocon(local_uplo, n, a, lda, anorm, rcond, work, xwork, local_info)
    else
        local_info = -1000
    end if

    deallocate(xwork, stat=deallocation_status)
    deallocate(work,  stat=deallocation_status)

    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('f77_pocon',-local_info)
    end if
end subroutine
!> Modern interface for [[f77_trtrs:f77_trtrs]].
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
    call f77_trtrs(local_uplo, local_trans, local_diag, n, nrhs, a, lda, b, ldb, local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_trtrs', local_info)
    end if
end subroutine
!> Modern interface for [[f77_trtrs:f77_trtrs]].
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
    call f77_trtrs(local_uplo, local_trans, local_diag, n, nrhs, a, lda, b, ldb, local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_trtrs', local_info)
    end if
end subroutine
!> Modern interface for [[f77_trtrs:f77_trtrs]].
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
    call f77_trtrs(local_uplo, local_trans, local_diag, n, nrhs, a, lda, b, ldb, local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_trtrs', local_info)
    end if
end subroutine
!> Modern interface for [[f77_trtrs:f77_trtrs]].
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
    call f77_trtrs(local_uplo, local_trans, local_diag, n, nrhs, a, lda, b, ldb, local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_trtrs', local_info)
    end if
end subroutine
!> Modern interface for [[f77_sytrf:f77_sytrf]].
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
    real(REAL32) :: s_work(1)  ! Work array for workspace query
    integer, pointer :: local_ipiv(:)
    real(REAL32), pointer :: work(:)
    local_ipiv => null()
    work => null()
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
    call f77_sytrf(local_uplo, n, a, lda, local_ipiv, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(s_work(1)))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_sytrf(local_uplo, n, a, lda, local_ipiv, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(ipiv)) then
        if (associated(local_ipiv)) deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_sytrf', abs(local_info))
    end if
end subroutine
!> Modern interface for [[f77_sytrf:f77_sytrf]].
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
    real(REAL64) :: s_work(1)  ! Work array for workspace query
    integer, pointer :: local_ipiv(:)
    real(REAL64), pointer :: work(:)
    local_ipiv => null()
    work => null()
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
    call f77_sytrf(local_uplo, n, a, lda, local_ipiv, s_work, lwork, local_info)
    if (local_info /= 0) goto 404

    lwork = max(1, int(s_work(1)))
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call f77_sytrf(local_uplo, n, a, lda, local_ipiv, work, lwork, local_info)
    else
        local_info = -1000
    end if
    if (associated(work)) deallocate(work, stat=deallocation_status)

    ! Error handling
404 continue
    if (.not. present(ipiv)) then
        if (associated(local_ipiv)) deallocate(local_ipiv, stat=deallocation_status)
    end if
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('f77_sytrf', abs(local_info))
    end if
end subroutine

    pure subroutine mfi_error(name, info)
        character(*), intent(in) :: name
        integer, intent(in) :: info
        call f77_xerbla(name, info)
    end subroutine

end module
