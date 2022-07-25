#:mute
#:include "common.fpp"
#:def gesvd(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, s, u, vt, ww, job, info)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,      inout, a(:,:))
@:args(${REAL_TYPE}$, out, s(:))
    ${TYPE}$,      intent(out), optional, target :: u(:,:), vt(:,:)
    ${REAL_TYPE}$, intent(out), optional, target :: ww(:)
@:optional(character, in, job)
@:optional(integer, out, info)
    character :: jobu, jobvt
    integer   :: m, n, lda, ldu, ldvt, lwork, allocation_status, deallocation_status
    ${TYPE}$,      target  :: s_work(1), l_a2(1,1)
    ${TYPE}$,      pointer :: local_u(:,:), local_vt(:,:), work(:)
#:if TYPE == COMPLEX_TYPE
    ${REAL_TYPE}$, pointer :: rwork(:)
#:endif
@:defaults(job='N', info=0)
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
@:optval(present(u),  ldu,  max(1,size(u,1)),  1)
@:optval(present(vt), ldvt, max(1,size(vt,1)), 1)
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
#:if TYPE == COMPLEX_TYPE
    allocate(rwork(5*min(m,n)), stat=allocation_status)
    call ${F77_NAME}$(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,rwork,local_info)
#:else
    call ${F77_NAME}$(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,s_work,lwork,local_info)
#:endif
    if (local_info /= 0) then
        goto 404
    end if
    lwork = s_work(1)
    allocate(work(lwork), stat=allocation_status)
    if (allocation_status == 0) then
#:if TYPE == COMPLEX_TYPE
        call ${F77_NAME}$(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,rwork,local_info)
#:else
        call ${F77_NAME}$(jobu,jobvt,m,n,a,lda,s,local_u,ldu,local_vt,ldvt,work,lwork,local_info)
#:endif
    else
        local_info = -1000
    end if

    if (present(ww)) then
        ww = work(2:min(m,n)-1)
    end if
    deallocate(work, stat=deallocation_status)
404 continue
#:if TYPE == COMPLEX_TYPE
    deallocate(rwork, stat=deallocation_status)
#:endif
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$', -local_info)
    end if
end subroutine
#:enddef

#:def hegv(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, b, w, itype, jobz, uplo, info)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, a(:,:), b(:,:))
@:args(${REAL_TYPE}$, out, w(:))
@:optional(integer,   in, itype)
@:optional(character, in, jobz, uplo)
@:optional(integer, out, info)
    ${TYPE}$,      pointer :: work(:)
    ${REAL_TYPE}$, pointer :: rwork(:)
    ${TYPE}$ :: s_work(1)
    integer :: n, lda, ldb, lwork, allocation_status, deallocation_status
@:defaults(itype=1, jobz='N', uplo='U', info=0)
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    n = size(a,2)
    allocation_status = 0
    allocate(rwork(max(1,3*N-2)), stat=allocation_status)
    lwork = -1
    call ${F77_NAME}$(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,s_work,lwork,rwork,local_info)
    if (local_info /= 0) goto 404
    lwork = s_work(1)
    if (allocation_status == 0) then
        allocate(work(lwork), stat=allocation_status)
    end if
    if (allocation_status == 0) then
        call ${F77_NAME}$(local_itype,local_jobz,local_uplo,n,a,lda,b,ldb,w,work,lwork,rwork,local_info)
    else
        local_info = -1000
    end if
    deallocate(work, stat=deallocation_status)
404 continue
    deallocate(rwork, stat=deallocation_status)
    if (present(info)) then
        info = local_info
    else if (local_info <= -1000) then
        call mfi_error('${F77_NAME}$', -local_info)
    end if
end subroutine
#:enddef

#:def heevd(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, w, jobz, uplo, info)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, a(:,:))
@:args(${REAL_TYPE}$, out, w(:))
@:optional(integer, out, info)
@:optional(character, in, jobz, uplo)
    ${TYPE}$,      pointer :: work(:)
    ${REAL_TYPE}$, pointer :: rwork(:)
    integer,       pointer :: iwork(:)
    ${TYPE}$      :: s_work(1)
    ${REAL_TYPE}$ :: s_rwork(1)
    integer       :: s_iwork(1)
    integer :: n, lda, lwork, lrwork, liwork, allocation_status, deallocation_status
@:defaults(jobz='N', uplo='U', info=0)
    lda = max(1,size(a,1))
    n   = size(a,2)
    allocation_status = 0
    lwork  = -1
    lrwork = -1
    liwork = -1

    call ${F77_NAME}$(local_jobz,local_uplo,n,a,lda,w, &
                      s_work,lwork,s_rwork,lrwork,s_iwork,liwork,local_info)
    if (local_info /= 0) goto 404
    lwork  = s_work(1)
    lrwork = s_rwork(1)
    liwork = s_iwork(1)

    allocate(iwork(liwork), stat=allocation_status)

    if (allocation_status == 0) then
        allocate(rwork(lrwork), stat=allocation_status)
        allocate(work(lwork),   stat=allocation_status)
        call ${F77_NAME}$(local_jobz,local_uplo,n,a,lda,w, &
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
        call mfi_error('${F77_NAME}$', -local_info)
    end if
end subroutine
#:enddef

#:def potrf_potri(MFI_NAME,F77_NAME,TYPE,KIND)
pure subroutine ${MFI_NAME}$(a, info, uplo)
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, a(:,:))
@:optional(character, in, uplo)
@:optional(integer, out, info)
    integer :: n, lda
@:defaults(uplo='U', info=0)
    lda = max(1,size(a,1))
    n = size(a,2)
    call ${F77_NAME}$(local_uplo,n,a,lda,local_info)
    if (present(info)) then
        info = local_info
    else if (local_info /= 0) then
        call mfi_error('${F77_NAME}$', local_info)
    end if
end subroutine
#:enddef
#:endmute
module mfi_lapack
use iso_fortran_env
use f77_lapack
implicit none

$:mfi_interface('?hegv',   COMPLEX_TYPES)
$:mfi_interface('?heevd',  COMPLEX_TYPES)
$:mfi_interface('?gesvd',  DEFAULT_TYPES)
$:mfi_interface('?potrf',  DEFAULT_TYPES)
$:mfi_interface('?potri',  DEFAULT_TYPES)

contains

$:mfi_implement('?hegv',   COMPLEX_TYPES, hegv)
$:mfi_implement('?heevd',  COMPLEX_TYPES, heevd)
$:mfi_implement('?gesvd',  DEFAULT_TYPES, gesvd)
$:mfi_implement('?potrf',  DEFAULT_TYPES, potrf_potri)
$:mfi_implement('?potri',  DEFAULT_TYPES, potrf_potri)

    pure subroutine mfi_error(name, info)
        character(*), intent(in) :: name
        integer, intent(in) :: info
        call f77_xerbla(name, info)
    end subroutine

end module
