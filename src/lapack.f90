!> Improved and original F77 interfaces for LAPACK
module f77_lapack
use iso_fortran_env
implicit none


!> ?geqrf supports s, d, c, z.
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
interface
pure subroutine sgeqrf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    real(REAL32), intent(inout) :: work(*)
end subroutine
pure subroutine dgeqrf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    real(REAL64), intent(inout) :: work(*)
end subroutine
pure subroutine cgeqrf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    complex(REAL32), intent(inout) :: work(*)
end subroutine
pure subroutine zgeqrf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    complex(REAL64), intent(inout) :: work(*)
end subroutine
end interface
!> ?gerqf supports s, d, c, z.
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
interface
pure subroutine sgerqf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    real(REAL32), intent(inout) :: work(*)
end subroutine
pure subroutine dgerqf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    real(REAL64), intent(inout) :: work(*)
end subroutine
pure subroutine cgerqf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    complex(REAL32), intent(inout) :: work(*)
end subroutine
pure subroutine zgerqf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    complex(REAL64), intent(inout) :: work(*)
end subroutine
end interface
!> ?getrf supports s, d, c, z.
!> See also: [[mfi_getrf]], [[f77_getrf]].
interface
pure subroutine sgetrf(m,n,a,lda,ipiv,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
end subroutine
pure subroutine dgetrf(m,n,a,lda,ipiv,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
end subroutine
pure subroutine cgetrf(m,n,a,lda,ipiv,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
end subroutine
pure subroutine zgetrf(m,n,a,lda,ipiv,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
end subroutine
end interface
!> ?getri supports s, d, c, z.
!> See also: [[mfi_getri]], [[f77_getri]].
interface
pure subroutine sgetri(n,a,lda,ipiv,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: work(*)
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
pure subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: work(*)
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
pure subroutine cgetri(n,a,lda,ipiv,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: work(*)
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
pure subroutine zgetri(n,a,lda,ipiv,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: work(*)
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?getrs supports s, d, c, z.
!> See also: [[mfi_getrs]], [[f77_getrs]].
interface
pure subroutine sgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    character, intent(in) :: trans
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine dgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    character, intent(in) :: trans
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine cgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    character, intent(in) :: trans
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine zgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    character, intent(in) :: trans
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
end interface
!> ?hetrf supports c, z.
!> See also: [[mfi_hetrf]], [[f77_hetrf]].
interface
pure subroutine chetrf(uplo, n, a, lda, ipiv, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: ipiv(*)
    complex(REAL32), intent(inout) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
pure subroutine zhetrf(uplo, n, a, lda, ipiv, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: ipiv(*)
    complex(REAL64), intent(inout) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?hegv supports c, z.
!> See also: [[mfi_hegv]], [[f77_hegv]].
interface
pure subroutine chegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    real(REAL32), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: itype
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
    complex(REAL32), intent(inout) :: work(*)
    real(REAL32), intent(in) :: rwork(*)
end subroutine
pure subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    real(REAL64), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: itype
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
    complex(REAL64), intent(inout) :: work(*)
    real(REAL64), intent(in) :: rwork(*)
end subroutine
end interface
!> ?heevd supports c, z.
!> See also: [[mfi_heevd]], [[f77_heevd]].
interface
pure subroutine cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(in) :: lrwork
    integer, intent(in) :: liwork
    complex(REAL32), intent(inout) :: work(*)
    real(REAL32), intent(inout) :: rwork(*)
    integer, intent(inout) :: iwork(*)
end subroutine
pure subroutine zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(in) :: lrwork
    integer, intent(in) :: liwork
    complex(REAL64), intent(inout) :: work(*)
    real(REAL64), intent(inout) :: rwork(*)
    integer, intent(inout) :: iwork(*)
end subroutine
end interface
!> ?gesvd supports s, d, c, z.
!> See also: [[mfi_gesvd]], [[f77_gesvd]].
interface
pure subroutine sgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(out) :: s(*)
    real(REAL32), intent(out) :: u(ldu,*)
    real(REAL32), intent(out) :: vt(ldvt,*)
    integer, intent(out) :: info
    character, intent(in) :: jobu
    character, intent(in) :: jobvt
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldu
    integer, intent(in) :: ldvt
    integer, intent(in) :: lwork
    real(REAL32), intent(inout) :: work(*)
end subroutine
pure subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(out) :: s(*)
    real(REAL64), intent(out) :: u(ldu,*)
    real(REAL64), intent(out) :: vt(ldvt,*)
    integer, intent(out) :: info
    character, intent(in) :: jobu
    character, intent(in) :: jobvt
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldu
    integer, intent(in) :: ldvt
    integer, intent(in) :: lwork
    real(REAL64), intent(inout) :: work(*)
end subroutine
pure subroutine cgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(out) :: s(*)
    complex(REAL32), intent(out) :: u(ldu,*)
    complex(REAL32), intent(out) :: vt(ldvt,*)
    integer, intent(out) :: info
    character, intent(in) :: jobu
    character, intent(in) :: jobvt
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldu
    integer, intent(in) :: ldvt
    integer, intent(in) :: lwork
    complex(REAL32), intent(inout) :: work(*)
    real(REAL32), intent(in) :: rwork(*)
end subroutine
pure subroutine zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(out) :: s(*)
    complex(REAL64), intent(out) :: u(ldu,*)
    complex(REAL64), intent(out) :: vt(ldvt,*)
    integer, intent(out) :: info
    character, intent(in) :: jobu
    character, intent(in) :: jobvt
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldu
    integer, intent(in) :: ldvt
    integer, intent(in) :: lwork
    complex(REAL64), intent(inout) :: work(*)
    real(REAL64), intent(in) :: rwork(*)
end subroutine
end interface
!> ?potrf supports s, d, c, z.
!> See also: [[mfi_potrf]], [[f77_potrf]].
interface
pure subroutine spotrf(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine dpotrf(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine cpotrf(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine zpotrf(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
end interface
!> ?potri supports s, d, c, z.
!> See also: [[mfi_potri]], [[f77_potri]].
interface
pure subroutine spotri(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine dpotri(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine cpotri(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine zpotri(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
end interface
!> ?potrs supports s, d, c, z.
!> See also: [[mfi_potrs]], [[f77_potrs]].
interface
pure subroutine spotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: b(ldb,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(out) :: info
end subroutine
pure subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: b(ldb,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(out) :: info
end subroutine
pure subroutine cpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: b(ldb,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(out) :: info
end subroutine
pure subroutine zpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: b(ldb,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(out) :: info
end subroutine
end interface
!> ?pocon supports s, d, c, z.
!> See also: [[mfi_pocon]], [[f77_pocon]].
interface
!> spocon estimates the reciprocal of the condition number (in the
!> 1-norm) of a real(REAL32) Hermitian positive definite matrix using the
!> Cholesky factorization A = U**H*U or A = L*L**H computed by sPOTRF.
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
pure subroutine spocon(uplo, n, a, lda, anorm, rcond, work, iwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(in) :: anorm
    real(REAL32), intent(out) :: rcond
    real(REAL32), intent(inout) :: work(*)
    integer, intent(inout) :: iwork(*)
    integer, intent(out) :: info
end subroutine
!> dpocon estimates the reciprocal of the condition number (in the
!> 1-norm) of a real(REAL64) Hermitian positive definite matrix using the
!> Cholesky factorization A = U**H*U or A = L*L**H computed by dPOTRF.
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
pure subroutine dpocon(uplo, n, a, lda, anorm, rcond, work, iwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(in) :: anorm
    real(REAL64), intent(out) :: rcond
    real(REAL64), intent(inout) :: work(*)
    integer, intent(inout) :: iwork(*)
    integer, intent(out) :: info
end subroutine
!> cpocon estimates the reciprocal of the condition number (in the
!> 1-norm) of a complex(REAL32) Hermitian positive definite matrix using the
!> Cholesky factorization A = U**H*U or A = L*L**H computed by cPOTRF.
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
pure subroutine cpocon(uplo, n, a, lda, anorm, rcond, work, rwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    complex(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(in) :: anorm
    real(REAL32), intent(out) :: rcond
    complex(REAL32), intent(inout) :: work(*)
    real(REAL32), intent(inout) :: rwork(*)
    integer, intent(out) :: info
end subroutine
!> zpocon estimates the reciprocal of the condition number (in the
!> 1-norm) of a complex(REAL64) Hermitian positive definite matrix using the
!> Cholesky factorization A = U**H*U or A = L*L**H computed by zPOTRF.
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
pure subroutine zpocon(uplo, n, a, lda, anorm, rcond, work, rwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    complex(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(in) :: anorm
    real(REAL64), intent(out) :: rcond
    complex(REAL64), intent(inout) :: work(*)
    real(REAL64), intent(inout) :: rwork(*)
    integer, intent(out) :: info
end subroutine
end interface
!> ?heevx supports c, z.
!> See also: [[mfi_heevx]], [[f77_heevx]].
interface
pure subroutine cheevx(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,&
                         work,lwork,rwork,lrwork,iwork,liwork,ifail,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: z(ldz, *)
    real(REAL32), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    character, intent(in) :: range
    real(REAL32), intent(in) :: vl
    real(REAL32), intent(in) :: vu
    real(REAL32), intent(in) :: abstol
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: lda
    integer, intent(in) :: ldz
    integer, intent(in) :: il
    integer, intent(in) :: iu
    integer, intent(in) :: lwork
    integer, intent(in) :: lrwork
    integer, intent(in) :: liwork
    integer, intent(in) :: ifail
    complex(REAL32), intent(inout) :: work(*)
    real(REAL32), intent(inout) :: rwork(*)
    integer, intent(inout) :: iwork(*)
end subroutine
pure subroutine zheevx(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,&
                         work,lwork,rwork,lrwork,iwork,liwork,ifail,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: z(ldz, *)
    real(REAL64), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    character, intent(in) :: range
    real(REAL64), intent(in) :: vl
    real(REAL64), intent(in) :: vu
    real(REAL64), intent(in) :: abstol
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: lda
    integer, intent(in) :: ldz
    integer, intent(in) :: il
    integer, intent(in) :: iu
    integer, intent(in) :: lwork
    integer, intent(in) :: lrwork
    integer, intent(in) :: liwork
    integer, intent(in) :: ifail
    complex(REAL64), intent(inout) :: work(*)
    real(REAL64), intent(inout) :: rwork(*)
    integer, intent(inout) :: iwork(*)
end subroutine
end interface
!> ?heevr supports c, z.
!> See also: [[mfi_heevr]], [[f77_heevr]].
interface
pure subroutine cheevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,&
                         isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: z(ldz, *)
    real(REAL32), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    character, intent(in) :: range
    real(REAL32), intent(in) :: vl
    real(REAL32), intent(in) :: vu
    real(REAL32), intent(in) :: abstol
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: lda
    integer, intent(in) :: ldz
    integer, intent(in) :: il
    integer, intent(in) :: iu
    integer, intent(in) :: lwork
    integer, intent(in) :: lrwork
    integer, intent(in) :: liwork
    integer, intent(in) :: isuppz(*)
    complex(REAL32), intent(inout) :: work(*)
    real(REAL32), intent(inout) :: rwork(*)
    integer, intent(inout) :: iwork(*)
end subroutine
pure subroutine zheevr(jobz,range,uplo,n,a,lda,vl,vu,il,iu,abstol,m,w,z,ldz,&
                         isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: z(ldz, *)
    real(REAL64), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    character, intent(in) :: range
    real(REAL64), intent(in) :: vl
    real(REAL64), intent(in) :: vu
    real(REAL64), intent(in) :: abstol
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: lda
    integer, intent(in) :: ldz
    integer, intent(in) :: il
    integer, intent(in) :: iu
    integer, intent(in) :: lwork
    integer, intent(in) :: lrwork
    integer, intent(in) :: liwork
    integer, intent(in) :: isuppz(*)
    complex(REAL64), intent(inout) :: work(*)
    real(REAL64), intent(inout) :: rwork(*)
    integer, intent(inout) :: iwork(*)
end subroutine
end interface
!> ?gels supports s, d, c, z.
!> See also: [[mfi_gels]], [[f77_gels]].
interface
!> SGELS solves overdetermined or underdetermined systems for GE matrices
pure subroutine sgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    character, intent(in) :: trans
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    real(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
!> DGELS solves overdetermined or underdetermined systems for GE matrices
pure subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    character, intent(in) :: trans
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    real(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
!> CGELS solves overdetermined or underdetermined systems for GE matrices
pure subroutine cgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    character, intent(in) :: trans
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    complex(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
!> ZGELS solves overdetermined or underdetermined systems for GE matrices
pure subroutine zgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    character, intent(in) :: trans
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    complex(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?gelst supports s, d, c, z.
!> See also: [[mfi_gelst]], [[f77_gelst]].
interface
!> SGELST solves overdetermined or underdetermined systems for GE matrices
!> using QR or LQ factorization with compact WY representation of Q.
pure subroutine sgelst(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    character, intent(in) :: trans
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    real(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
!> DGELST solves overdetermined or underdetermined systems for GE matrices
!> using QR or LQ factorization with compact WY representation of Q.
pure subroutine dgelst(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    character, intent(in) :: trans
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    real(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
!> CGELST solves overdetermined or underdetermined systems for GE matrices
!> using QR or LQ factorization with compact WY representation of Q.
pure subroutine cgelst(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    character, intent(in) :: trans
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    complex(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
!> ZGELST solves overdetermined or underdetermined systems for GE matrices
!> using QR or LQ factorization with compact WY representation of Q.
pure subroutine zgelst(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    character, intent(in) :: trans
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    complex(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?getsls supports s, d, c, z.
!> See also: [[mfi_getsls]], [[f77_getsls]].
interface
!> SGETSLS solves overdetermined or underdetermined systems for GE matrices
pure subroutine sgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    character, intent(in) :: trans
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    real(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
!> DGETSLS solves overdetermined or underdetermined systems for GE matrices
pure subroutine dgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    character, intent(in) :: trans
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    real(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
!> CGETSLS solves overdetermined or underdetermined systems for GE matrices
pure subroutine cgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    character, intent(in) :: trans
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    complex(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
!> ZGETSLS solves overdetermined or underdetermined systems for GE matrices
pure subroutine zgetsls(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    character, intent(in) :: trans
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    complex(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: nrhs
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?gelsd supports s, d, c, z.
!> See also: [[mfi_gelsd]], [[f77_gelsd]].
interface
!> SGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices
pure subroutine sgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: rcond
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    real(REAL32), intent(out) :: s(*)
    real(REAL32), intent(out) :: work(*)
    integer, intent(out) :: iwork(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
!> DGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices
pure subroutine dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: rcond
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    real(REAL64), intent(out) :: s(*)
    real(REAL64), intent(out) :: work(*)
    integer, intent(out) :: iwork(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
!> CGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices
pure subroutine cgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: rcond
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    complex(REAL32), intent(out) :: s(*)
    complex(REAL32), intent(out) :: work(*)
    integer, intent(out) :: iwork(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
!> ZGELSD computes the minimum-norm solution to a linear least squares problem for GE matrices
pure subroutine zgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: rcond
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    complex(REAL64), intent(out) :: s(*)
    complex(REAL64), intent(out) :: work(*)
    integer, intent(out) :: iwork(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?gelss supports s, d, c, z.
!> See also: [[mfi_gelss]], [[f77_gelss]].
interface
!> SGELSS solves overdetermined or underdetermined systems for GE matrices
pure subroutine sgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: rcond
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    real(REAL32), intent(out) :: s(*)
    real(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
!> DGELSS solves overdetermined or underdetermined systems for GE matrices
pure subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: rcond
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    real(REAL64), intent(out) :: s(*)
    real(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
!> CGELSS solves overdetermined or underdetermined systems for GE matrices
pure subroutine cgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: rcond
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    complex(REAL32), intent(out) :: s(*)
    complex(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
!> ZGELSS solves overdetermined or underdetermined systems for GE matrices
pure subroutine zgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: rcond
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    complex(REAL64), intent(out) :: s(*)
    complex(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?gelsy supports s, d, c, z.
!> See also: [[mfi_gelsy]], [[f77_gelsy]].
interface
!> SGELSY solves overdetermined or underdetermined systems for GE matrices
pure subroutine sgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: rcond
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    integer, intent(inout) :: jpvt(*)
    real(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
!> DGELSY solves overdetermined or underdetermined systems for GE matrices
pure subroutine dgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: rcond
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    integer, intent(inout) :: jpvt(*)
    real(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
!> CGELSY solves overdetermined or underdetermined systems for GE matrices
pure subroutine cgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: rcond
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    integer, intent(inout) :: jpvt(*)
    complex(REAL32), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
!> ZGELSY solves overdetermined or underdetermined systems for GE matrices
pure subroutine zgelsy(m, n, nrhs, a, lda, b, ldb, jpvt, rcond, rank, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: rcond
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    integer, intent(inout) :: jpvt(*)
    complex(REAL64), intent(out) :: work(*)
    integer, intent(out) :: info
    integer, intent(out) :: rank
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?gglse supports s, d, c, z.
!> See also: [[mfi_gglse]], [[f77_gglse]].
interface
pure subroutine sgglse(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    real(REAL32), intent(inout) :: c(*)
    real(REAL32), intent(inout) :: d(*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(out) :: x(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: p
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
pure subroutine dgglse(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    real(REAL64), intent(inout) :: c(*)
    real(REAL64), intent(inout) :: d(*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(out) :: x(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: p
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
pure subroutine cgglse(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    complex(REAL32), intent(inout) :: c(*)
    complex(REAL32), intent(inout) :: d(*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(out) :: x(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: p
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
pure subroutine zgglse(m, n, p, a, lda, b, ldb, c, d, x, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    complex(REAL64), intent(inout) :: c(*)
    complex(REAL64), intent(inout) :: d(*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(out) :: x(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: p
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?gglsm supports s, d, c, z.
!> See also: [[mfi_gglsm]], [[f77_gglsm]].
interface
pure subroutine sgglsm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    real(REAL32), intent(inout) :: d(*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(out) :: x(*)
    real(REAL32), intent(out) :: y(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: p
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
pure subroutine dgglsm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    real(REAL64), intent(inout) :: d(*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(out) :: x(*)
    real(REAL64), intent(out) :: y(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: p
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
pure subroutine cgglsm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    complex(REAL32), intent(inout) :: d(*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(out) :: x(*)
    complex(REAL32), intent(out) :: y(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: p
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
pure subroutine zgglsm(n, m, p, a, lda, b, ldb, d, x, y, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    complex(REAL64), intent(inout) :: d(*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(out) :: x(*)
    complex(REAL64), intent(out) :: y(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: p
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
end subroutine
end interface
!> ?lartg supports s, d, c, z.
!> See also: [[mfi_lartg]], [[f77_lartg]].
interface
pure subroutine slartg(f, g, c, s, r)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: c
    real(REAL32), intent(inout) :: f
    real(REAL32), intent(inout) :: g
    real(REAL32), intent(inout) :: r
    real(REAL32), intent(inout) :: s
end subroutine
pure subroutine dlartg(f, g, c, s, r)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: c
    real(REAL64), intent(inout) :: f
    real(REAL64), intent(inout) :: g
    real(REAL64), intent(inout) :: r
    real(REAL64), intent(inout) :: s
end subroutine
pure subroutine clartg(f, g, c, s, r)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: c
    complex(REAL32), intent(inout) :: f
    complex(REAL32), intent(inout) :: g
    complex(REAL32), intent(inout) :: r
    complex(REAL32), intent(inout) :: s
end subroutine
pure subroutine zlartg(f, g, c, s, r)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: c
    complex(REAL64), intent(inout) :: f
    complex(REAL64), intent(inout) :: g
    complex(REAL64), intent(inout) :: r
    complex(REAL64), intent(inout) :: s
end subroutine
end interface

interface f77_geqrf
    procedure :: sgeqrf
    procedure :: dgeqrf
    procedure :: cgeqrf
    procedure :: zgeqrf
end interface
interface f77_gerqf
    procedure :: sgerqf
    procedure :: dgerqf
    procedure :: cgerqf
    procedure :: zgerqf
end interface
interface f77_getrf
    procedure :: sgetrf
    procedure :: dgetrf
    procedure :: cgetrf
    procedure :: zgetrf
end interface
interface f77_getri
    procedure :: sgetri
    procedure :: dgetri
    procedure :: cgetri
    procedure :: zgetri
end interface
interface f77_getrs
    procedure :: sgetrs
    procedure :: dgetrs
    procedure :: cgetrs
    procedure :: zgetrs
end interface
interface f77_hetrf
    procedure :: chetrf
    procedure :: zhetrf
end interface
interface f77_hegv
    procedure :: chegv
    procedure :: zhegv
end interface
interface f77_heevd
    procedure :: cheevd
    procedure :: zheevd
end interface
interface f77_gesvd
    procedure :: sgesvd
    procedure :: dgesvd
    procedure :: cgesvd
    procedure :: zgesvd
end interface
interface f77_potrf
    procedure :: spotrf
    procedure :: dpotrf
    procedure :: cpotrf
    procedure :: zpotrf
end interface
interface f77_potri
    procedure :: spotri
    procedure :: dpotri
    procedure :: cpotri
    procedure :: zpotri
end interface
interface f77_potrs
    procedure :: spotrs
    procedure :: dpotrs
    procedure :: cpotrs
    procedure :: zpotrs
end interface
interface f77_pocon
    procedure :: spocon
    procedure :: dpocon
    procedure :: cpocon
    procedure :: zpocon
end interface
interface f77_heevx
    procedure :: cheevx
    procedure :: zheevx
end interface
interface f77_heevr
    procedure :: cheevr
    procedure :: zheevr
end interface
interface f77_gels
    procedure :: sgels
    procedure :: dgels
    procedure :: cgels
    procedure :: zgels
end interface
interface f77_gelst
    procedure :: sgelst
    procedure :: dgelst
    procedure :: cgelst
    procedure :: zgelst
end interface
interface f77_getsls
    procedure :: sgetsls
    procedure :: dgetsls
    procedure :: cgetsls
    procedure :: zgetsls
end interface
interface f77_gelsd
    procedure :: sgelsd
    procedure :: dgelsd
    procedure :: cgelsd
    procedure :: zgelsd
end interface
interface f77_gelss
    procedure :: sgelss
    procedure :: dgelss
    procedure :: cgelss
    procedure :: zgelss
end interface
interface f77_gelsy
    procedure :: sgelsy
    procedure :: dgelsy
    procedure :: cgelsy
    procedure :: zgelsy
end interface
interface f77_gglse
    procedure :: sgglse
    procedure :: dgglse
    procedure :: cgglse
    procedure :: zgglse
end interface
interface f77_gglsm
    procedure :: sgglsm
    procedure :: dgglsm
    procedure :: cgglsm
    procedure :: zgglsm
end interface
interface f77_lartg
    procedure :: slartg
    procedure :: dlartg
    procedure :: clartg
    procedure :: zlartg
end interface

    interface f77_xerbla
        pure subroutine xerbla(name,info)
            character(*), intent(in) :: name
            integer, intent(in) :: info
        end subroutine
    end interface f77_xerbla

end module

