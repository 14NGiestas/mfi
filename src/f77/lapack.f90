!> Improved and original F77 interfaces for LAPACK
module f77_lapack
use iso_fortran_env
implicit none

!> Generic old style interface for GEQRF.
!> Supports s, d, c, z.
!> See also: [[mfi_geqrf]], [[sgeqrf]], [[dgeqrf]], [[cgeqrf]], [[zgeqrf]].
interface f77_geqrf
!> Original interface for SGEQRF
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
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
!> Original interface for DGEQRF
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
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
!> Original interface for CGEQRF
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
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
!> Original interface for ZGEQRF
!> See also: [[mfi_geqrf]], [[f77_geqrf]].
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
!> Generic old style interface for GERQF.
!> Supports s, d, c, z.
!> See also: [[mfi_gerqf]], [[sgerqf]], [[dgerqf]], [[cgerqf]], [[zgerqf]].
interface f77_gerqf
!> Original interface for SGERQF
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
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
!> Original interface for DGERQF
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
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
!> Original interface for CGERQF
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
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
!> Original interface for ZGERQF
!> See also: [[mfi_gerqf]], [[f77_gerqf]].
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
!> Generic old style interface for GETRF.
!> Supports s, d, c, z.
!> See also: [[mfi_getrf]], [[sgetrf]], [[dgetrf]], [[cgetrf]], [[zgetrf]].
interface f77_getrf
!> Original interface for SGETRF
!> See also: [[mfi_getrf]], [[f77_getrf]].
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
!> Original interface for DGETRF
!> See also: [[mfi_getrf]], [[f77_getrf]].
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
!> Original interface for CGETRF
!> See also: [[mfi_getrf]], [[f77_getrf]].
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
!> Original interface for ZGETRF
!> See also: [[mfi_getrf]], [[f77_getrf]].
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
!> Generic old style interface for GETRI.
!> Supports s, d, c, z.
!> See also: [[mfi_getri]], [[sgetri]], [[dgetri]], [[cgetri]], [[zgetri]].
interface f77_getri
!> Original interface for SGETRI
!> See also: [[mfi_getri]], [[f77_getri]].
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
!> Original interface for DGETRI
!> See also: [[mfi_getri]], [[f77_getri]].
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
!> Original interface for CGETRI
!> See also: [[mfi_getri]], [[f77_getri]].
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
!> Original interface for ZGETRI
!> See also: [[mfi_getri]], [[f77_getri]].
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
!> Generic old style interface for GETRS.
!> Supports s, d, c, z.
!> See also: [[mfi_getrs]], [[sgetrs]], [[dgetrs]], [[cgetrs]], [[zgetrs]].
interface f77_getrs
!> Original interface for SGETRS
!> See also: [[mfi_getrs]], [[f77_getrs]].
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
!> Original interface for DGETRS
!> See also: [[mfi_getrs]], [[f77_getrs]].
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
!> Original interface for CGETRS
!> See also: [[mfi_getrs]], [[f77_getrs]].
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
!> Original interface for ZGETRS
!> See also: [[mfi_getrs]], [[f77_getrs]].
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
!> Generic old style interface for HETRF.
!> Supports c, z.
!> See also: [[mfi_hetrf]], [[chetrf]], [[zhetrf]].
interface f77_hetrf
!> Original interface for CHETRF
!> See also: [[mfi_hetrf]], [[f77_hetrf]].
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
!> Original interface for ZHETRF
!> See also: [[mfi_hetrf]], [[f77_hetrf]].
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
!> Generic old style interface for HEGV.
!> Supports c, z.
!> See also: [[mfi_hegv]], [[chegv]], [[zhegv]].
interface f77_hegv
!> Original interface for CHEGV
!> See also: [[mfi_hegv]], [[f77_hegv]].
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
!> Original interface for ZHEGV
!> See also: [[mfi_hegv]], [[f77_hegv]].
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
!> Generic old style interface for HEEVD.
!> Supports c, z.
!> See also: [[mfi_heevd]], [[cheevd]], [[zheevd]].
interface f77_heevd
!> Original interface for CHEEVD
!> See also: [[mfi_heevd]], [[f77_heevd]].
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
!> Original interface for ZHEEVD
!> See also: [[mfi_heevd]], [[f77_heevd]].
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
!> Generic old style interface for GESVD.
!> Supports s, d, c, z.
!> See also: [[mfi_gesvd]], [[sgesvd]], [[dgesvd]], [[cgesvd]], [[zgesvd]].
interface f77_gesvd
!> Original interface for SGESVD
!> See also: [[mfi_gesvd]], [[f77_gesvd]].
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
!> Original interface for DGESVD
!> See also: [[mfi_gesvd]], [[f77_gesvd]].
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
!> Original interface for CGESVD
!> See also: [[mfi_gesvd]], [[f77_gesvd]].
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
!> Original interface for ZGESVD
!> See also: [[mfi_gesvd]], [[f77_gesvd]].
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
!> Generic old style interface for POTRF.
!> Supports s, d, c, z.
!> See also: [[mfi_potrf]], [[spotrf]], [[dpotrf]], [[cpotrf]], [[zpotrf]].
interface f77_potrf
!> Original interface for SPOTRF
!> See also: [[mfi_potrf]], [[f77_potrf]].
pure subroutine spotrf(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for DPOTRF
!> See also: [[mfi_potrf]], [[f77_potrf]].
pure subroutine dpotrf(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for CPOTRF
!> See also: [[mfi_potrf]], [[f77_potrf]].
pure subroutine cpotrf(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for ZPOTRF
!> See also: [[mfi_potrf]], [[f77_potrf]].
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
!> Generic old style interface for POTRI.
!> Supports s, d, c, z.
!> See also: [[mfi_potri]], [[spotri]], [[dpotri]], [[cpotri]], [[zpotri]].
interface f77_potri
!> Original interface for SPOTRI
!> See also: [[mfi_potri]], [[f77_potri]].
pure subroutine spotri(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for DPOTRI
!> See also: [[mfi_potri]], [[f77_potri]].
pure subroutine dpotri(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for CPOTRI
!> See also: [[mfi_potri]], [[f77_potri]].
pure subroutine cpotri(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for ZPOTRI
!> See also: [[mfi_potri]], [[f77_potri]].
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
!> Generic old style interface for POTRS.
!> Supports s, d, c, z.
!> See also: [[mfi_potrs]], [[spotrs]], [[dpotrs]], [[cpotrs]], [[zpotrs]].
interface f77_potrs
!> Original interface for SPOTRS
!> See also: [[mfi_potrs]], [[f77_potrs]].
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
!> Original interface for DPOTRS
!> See also: [[mfi_potrs]], [[f77_potrs]].
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
!> Original interface for CPOTRS
!> See also: [[mfi_potrs]], [[f77_potrs]].
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
!> Original interface for ZPOTRS
!> See also: [[mfi_potrs]], [[f77_potrs]].
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
!> Generic old style interface for POCON.
!> Supports s, d, c, z.
!> See also: [[mfi_pocon]], [[spocon]], [[dpocon]], [[cpocon]], [[zpocon]].
interface f77_pocon
!> Original interface for SPOCON
!> See also: [[mfi_pocon]], [[f77_pocon]].
!> spocon estimates the reciprocal of the condition number (in the
!> 1-norm) of a real(REAL32) Hermitian positive definite matrix using the
!> Cholesky factorization \( A = U^\dagger U \) or \( A = LL^\dagger |) computed by SPOTRF.
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
!> Original interface for DPOCON
!> See also: [[mfi_pocon]], [[f77_pocon]].
!> dpocon estimates the reciprocal of the condition number (in the
!> 1-norm) of a real(REAL64) Hermitian positive definite matrix using the
!> Cholesky factorization \( A = U^\dagger U \) or \( A = LL^\dagger |) computed by DPOTRF.
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
!> Original interface for CPOCON
!> See also: [[mfi_pocon]], [[f77_pocon]].
!> cpocon estimates the reciprocal of the condition number (in the
!> 1-norm) of a complex(REAL32) Hermitian positive definite matrix using the
!> Cholesky factorization \( A = U^\dagger U \) or \( A = LL^\dagger |) computed by CPOTRF.
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
!> Original interface for ZPOCON
!> See also: [[mfi_pocon]], [[f77_pocon]].
!> zpocon estimates the reciprocal of the condition number (in the
!> 1-norm) of a complex(REAL64) Hermitian positive definite matrix using the
!> Cholesky factorization \( A = U^\dagger U \) or \( A = LL^\dagger |) computed by ZPOTRF.
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
!> Generic old style interface for HEEVX.
!> Supports c, z.
!> See also: [[mfi_heevx]], [[cheevx]], [[zheevx]].
interface f77_heevx
!> Original interface for CHEEVX
!> See also: [[mfi_heevx]], [[f77_heevx]].
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
!> Original interface for ZHEEVX
!> See also: [[mfi_heevx]], [[f77_heevx]].
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
!> Generic old style interface for HEEVR.
!> Supports c, z.
!> See also: [[mfi_heevr]], [[cheevr]], [[zheevr]].
interface f77_heevr
!> Original interface for CHEEVR
!> See also: [[mfi_heevr]], [[f77_heevr]].
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
!> Original interface for ZHEEVR
!> See also: [[mfi_heevr]], [[f77_heevr]].
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
!> Generic old style interface for GELS.
!> Supports s, d, c, z.
!> See also: [[mfi_gels]], [[sgels]], [[dgels]], [[cgels]], [[zgels]].
interface f77_gels
!> Original interface for SGELS
!> See also: [[mfi_gels]], [[f77_gels]].
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
!> Original interface for DGELS
!> See also: [[mfi_gels]], [[f77_gels]].
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
!> Original interface for CGELS
!> See also: [[mfi_gels]], [[f77_gels]].
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
!> Original interface for ZGELS
!> See also: [[mfi_gels]], [[f77_gels]].
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
!> Generic old style interface for GELST.
!> Supports s, d, c, z.
!> See also: [[mfi_gelst]], [[sgelst]], [[dgelst]], [[cgelst]], [[zgelst]].
interface f77_gelst
!> Original interface for SGELST
!> See also: [[mfi_gelst]], [[f77_gelst]].
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
!> Original interface for DGELST
!> See also: [[mfi_gelst]], [[f77_gelst]].
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
!> Original interface for CGELST
!> See also: [[mfi_gelst]], [[f77_gelst]].
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
!> Original interface for ZGELST
!> See also: [[mfi_gelst]], [[f77_gelst]].
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
!> Generic old style interface for GETSLS.
!> Supports s, d, c, z.
!> See also: [[mfi_getsls]], [[sgetsls]], [[dgetsls]], [[cgetsls]], [[zgetsls]].
interface f77_getsls
!> Original interface for SGETSLS
!> See also: [[mfi_getsls]], [[f77_getsls]].
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
!> Original interface for DGETSLS
!> See also: [[mfi_getsls]], [[f77_getsls]].
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
!> Original interface for CGETSLS
!> See also: [[mfi_getsls]], [[f77_getsls]].
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
!> Original interface for ZGETSLS
!> See also: [[mfi_getsls]], [[f77_getsls]].
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
!> Generic old style interface for GELSD.
!> Supports s, d, c, z.
!> See also: [[mfi_gelsd]], [[sgelsd]], [[dgelsd]], [[cgelsd]], [[zgelsd]].
interface f77_gelsd
!> Original interface for SGELSD
!> See also: [[mfi_gelsd]], [[f77_gelsd]].
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
!> Original interface for DGELSD
!> See also: [[mfi_gelsd]], [[f77_gelsd]].
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
!> Original interface for CGELSD
!> See also: [[mfi_gelsd]], [[f77_gelsd]].
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
!> Original interface for ZGELSD
!> See also: [[mfi_gelsd]], [[f77_gelsd]].
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
!> Generic old style interface for GELSS.
!> Supports s, d, c, z.
!> See also: [[mfi_gelss]], [[sgelss]], [[dgelss]], [[cgelss]], [[zgelss]].
interface f77_gelss
!> Original interface for SGELSS
!> See also: [[mfi_gelss]], [[f77_gelss]].
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
!> Original interface for DGELSS
!> See also: [[mfi_gelss]], [[f77_gelss]].
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
!> Original interface for CGELSS
!> See also: [[mfi_gelss]], [[f77_gelss]].
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
!> Original interface for ZGELSS
!> See also: [[mfi_gelss]], [[f77_gelss]].
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
!> Generic old style interface for GELSY.
!> Supports s, d, c, z.
!> See also: [[mfi_gelsy]], [[sgelsy]], [[dgelsy]], [[cgelsy]], [[zgelsy]].
interface f77_gelsy
!> Original interface for SGELSY
!> See also: [[mfi_gelsy]], [[f77_gelsy]].
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
!> Original interface for DGELSY
!> See also: [[mfi_gelsy]], [[f77_gelsy]].
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
!> Original interface for CGELSY
!> See also: [[mfi_gelsy]], [[f77_gelsy]].
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
!> Original interface for ZGELSY
!> See also: [[mfi_gelsy]], [[f77_gelsy]].
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
!> Generic old style interface for GGLSE.
!> Supports s, d, c, z.
!> See also: [[mfi_gglse]], [[sgglse]], [[dgglse]], [[cgglse]], [[zgglse]].
interface f77_gglse
!> Original interface for SGGLSE
!> See also: [[mfi_gglse]], [[f77_gglse]].
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
!> Original interface for DGGLSE
!> See also: [[mfi_gglse]], [[f77_gglse]].
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
!> Original interface for CGGLSE
!> See also: [[mfi_gglse]], [[f77_gglse]].
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
!> Original interface for ZGGLSE
!> See also: [[mfi_gglse]], [[f77_gglse]].
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
!> Generic old style interface for GGLSM.
!> Supports s, d, c, z.
!> See also: [[mfi_gglsm]], [[sgglsm]], [[dgglsm]], [[cgglsm]], [[zgglsm]].
interface f77_gglsm
!> Original interface for SGGLSM
!> See also: [[mfi_gglsm]], [[f77_gglsm]].
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
!> Original interface for DGGLSM
!> See also: [[mfi_gglsm]], [[f77_gglsm]].
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
!> Original interface for CGGLSM
!> See also: [[mfi_gglsm]], [[f77_gglsm]].
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
!> Original interface for ZGGLSM
!> See also: [[mfi_gglsm]], [[f77_gglsm]].
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
!> Generic old style interface for ORG2R.
!> Supports s, d.
!> See also: [[mfi_org2r]], [[sorg2r]], [[dorg2r]].
interface f77_org2r
!> Original interface for SORG2R
!> See also: [[mfi_org2r]], [[f77_org2r]].
!> This routine generates an \(M \times N \) real(REAL32)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:sgeqrf]].
pure subroutine sorg2r(m, n, k, a, lda, tau, work, info)
    import :: REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for DORG2R
!> See also: [[mfi_org2r]], [[f77_org2r]].
!> This routine generates an \(M \times N \) real(REAL64)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:dgeqrf]].
pure subroutine dorg2r(m, n, k, a, lda, tau, work, info)
    import :: REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for ORGR2.
!> Supports s, d.
!> See also: [[mfi_orgr2]], [[sorgr2]], [[dorgr2]].
interface f77_orgr2
!> Original interface for SORGR2
!> See also: [[mfi_orgr2]], [[f77_orgr2]].
!> This routine generates an \(M \times N \) real(REAL32)
!> matrix \( Q \) with orthonormal rows,
!> which is defined as the last \( M \) rows of a product of \( K \) elementary
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:sgerqf]].
pure subroutine sorgr2(m, n, k, a, lda, tau, work, info)
    import :: REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for DORGR2
!> See also: [[mfi_orgr2]], [[f77_orgr2]].
!> This routine generates an \(M \times N \) real(REAL64)
!> matrix \( Q \) with orthonormal rows,
!> which is defined as the last \( M \) rows of a product of \( K \) elementary
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:dgerqf]].
pure subroutine dorgr2(m, n, k, a, lda, tau, work, info)
    import :: REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for ORM2R.
!> Supports s, d.
!> See also: [[mfi_orm2r]], [[sorm2r]], [[dorm2r]].
interface f77_orm2r
!> Original interface for SORM2R
!> See also: [[mfi_orm2r]], [[f77_orm2r]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:sgeqrf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine sorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    import :: REAL32
    character, intent(in) :: side
    character, intent(in) :: trans
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(out) :: info
end subroutine
!> Original interface for DORM2R
!> See also: [[mfi_orm2r]], [[f77_orm2r]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:dgeqrf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine dorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    import :: REAL64
    character, intent(in) :: side
    character, intent(in) :: trans
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for ORMR2.
!> Supports s, d.
!> See also: [[mfi_ormr2]], [[sormr2]], [[dormr2]].
interface f77_ormr2
!> Original interface for SORMR2
!> See also: [[mfi_ormr2]], [[f77_ormr2]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:sgerqf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine sormr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    import :: REAL32
    character, intent(in) :: side
    character, intent(in) :: trans
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(out) :: info
end subroutine
!> Original interface for DORMR2
!> See also: [[mfi_ormr2]], [[f77_ormr2]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:dgerqf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine dormr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    import :: REAL64
    character, intent(in) :: side
    character, intent(in) :: trans
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for ORMQR.
!> Supports s, d.
!> See also: [[mfi_ormqr]], [[sormqr]], [[dormqr]].
interface f77_ormqr
!> Original interface for SORMQR
!> See also: [[mfi_ormqr]], [[f77_ormqr]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:sgeqrf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine sormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    import :: REAL32
    character, intent(in) :: side
    character, intent(in) :: trans
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
!> Original interface for DORMQR
!> See also: [[mfi_ormqr]], [[f77_ormqr]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:dgeqrf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    import :: REAL64
    character, intent(in) :: side
    character, intent(in) :: trans
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for ORMRQ.
!> Supports s, d.
!> See also: [[mfi_ormrq]], [[sormrq]], [[dormrq]].
interface f77_ormrq
!> Original interface for SORMRQ
!> See also: [[mfi_ormrq]], [[f77_ormrq]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:sgerqf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine sormrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    import :: REAL32
    character, intent(in) :: side
    character, intent(in) :: trans
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
!> Original interface for DORMRQ
!> See also: [[mfi_ormrq]], [[f77_ormrq]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:dgerqf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine dormrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    import :: REAL64
    character, intent(in) :: side
    character, intent(in) :: trans
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for ORGQR.
!> Supports s, d.
!> See also: [[mfi_orgqr]], [[sorgqr]], [[dorgqr]].
interface f77_orgqr
!> Original interface for SORGQR
!> See also: [[mfi_orgqr]], [[f77_orgqr]].
!> This routine generates an \(M \times N \) real(REAL32)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:sgeqrf]].
pure subroutine sorgqr(m, n, k, a, lda, tau, work, lwork, info)
    import :: REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
!> Original interface for DORGQR
!> See also: [[mfi_orgqr]], [[f77_orgqr]].
!> This routine generates an \(M \times N \) real(REAL64)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:dgeqrf]].
pure subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
    import :: REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for ORGRQ.
!> Supports s, d.
!> See also: [[mfi_orgrq]], [[sorgrq]], [[dorgrq]].
interface f77_orgrq
!> Original interface for SORGRQ
!> See also: [[mfi_orgrq]], [[f77_orgrq]].
!> This routine generates an \(M \times N \) real(REAL32)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger . . . H(k)^\dagger \)
!> as returned by [[f77_gerqf:sgerqf]].
pure subroutine sorgrq(m, n, k, a, lda, tau, work, lwork, info)
    import :: REAL32
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(out) :: work(*)
    real(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
!> Original interface for DORGRQ
!> See also: [[mfi_orgrq]], [[f77_orgrq]].
!> This routine generates an \(M \times N \) real(REAL64)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger . . . H(k)^\dagger \)
!> as returned by [[f77_gerqf:dgerqf]].
pure subroutine dorgrq(m, n, k, a, lda, tau, work, lwork, info)
    import :: REAL64
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(out) :: work(*)
    real(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for UNG2R.
!> Supports c, z.
!> See also: [[mfi_ung2r]], [[cung2r]], [[zung2r]].
interface f77_ung2r
!> Original interface for CUNG2R
!> See also: [[mfi_ung2r]], [[f77_ung2r]].
!> This routine generates an \(M \times N \) complex(REAL32)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:cgeqrf]].
pure subroutine cung2r(m, n, k, a, lda, tau, work, info)
    import :: REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for ZUNG2R
!> See also: [[mfi_ung2r]], [[f77_ung2r]].
!> This routine generates an \(M \times N \) complex(REAL64)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:zgeqrf]].
pure subroutine zung2r(m, n, k, a, lda, tau, work, info)
    import :: REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for UNGR2.
!> Supports c, z.
!> See also: [[mfi_ungr2]], [[cungr2]], [[zungr2]].
interface f77_ungr2
!> Original interface for CUNGR2
!> See also: [[mfi_ungr2]], [[f77_ungr2]].
!> This routine generates an \(M \times N \) complex(REAL32)
!> matrix \( Q \) with orthonormal rows,
!> which is defined as the last \( M \) rows of a product of \( K \) elementary
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:cgerqf]].
pure subroutine cungr2(m, n, k, a, lda, tau, work, info)
    import :: REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
!> Original interface for ZUNGR2
!> See also: [[mfi_ungr2]], [[f77_ungr2]].
!> This routine generates an \(M \times N \) complex(REAL64)
!> matrix \( Q \) with orthonormal rows,
!> which is defined as the last \( M \) rows of a product of \( K \) elementary
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:zgerqf]].
pure subroutine zungr2(m, n, k, a, lda, tau, work, info)
    import :: REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for UNM2R.
!> Supports c, z.
!> See also: [[mfi_unm2r]], [[cunm2r]], [[zunm2r]].
interface f77_unm2r
!> Original interface for CUNM2R
!> See also: [[mfi_unm2r]], [[f77_unm2r]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:cgeqrf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine cunm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    import :: REAL32
    character, intent(in) :: side
    character, intent(in) :: trans
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(out) :: info
end subroutine
!> Original interface for ZUNM2R
!> See also: [[mfi_unm2r]], [[f77_unm2r]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:zgeqrf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine zunm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    import :: REAL64
    character, intent(in) :: side
    character, intent(in) :: trans
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for UNMR2.
!> Supports c, z.
!> See also: [[mfi_unmr2]], [[cunmr2]], [[zunmr2]].
interface f77_unmr2
!> Original interface for CUNMR2
!> See also: [[mfi_unmr2]], [[f77_unmr2]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:cgerqf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine cunmr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    import :: REAL32
    character, intent(in) :: side
    character, intent(in) :: trans
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(out) :: info
end subroutine
!> Original interface for ZUNMR2
!> See also: [[mfi_unmr2]], [[f77_unmr2]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:zgerqf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine zunmr2(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
    import :: REAL64
    character, intent(in) :: side
    character, intent(in) :: trans
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for UNMQR.
!> Supports c, z.
!> See also: [[mfi_unmqr]], [[cunmqr]], [[zunmqr]].
interface f77_unmqr
!> Original interface for CUNMQR
!> See also: [[mfi_unmqr]], [[f77_unmqr]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:cgeqrf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine cunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    import :: REAL32
    character, intent(in) :: side
    character, intent(in) :: trans
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
!> Original interface for ZUNMQR
!> See also: [[mfi_unmqr]], [[f77_unmqr]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1) H(2) \cdots H(k) \)
!> as returned by [[f77_geqrf:zgeqrf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine zunmqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    import :: REAL64
    character, intent(in) :: side
    character, intent(in) :: trans
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for UNMRQ.
!> Supports c, z.
!> See also: [[mfi_unmrq]], [[cunmrq]], [[zunmrq]].
interface f77_unmrq
!> Original interface for CUNMRQ
!> See also: [[mfi_unmrq]], [[f77_unmrq]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:cgerqf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine cunmrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    import :: REAL32
    character, intent(in) :: side
    character, intent(in) :: trans
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
!> Original interface for ZUNMRQ
!> See also: [[mfi_unmrq]], [[f77_unmrq]].
!> This routine overwrites the general complex \(M \times N\) matrix \( C \) with
!>```fortran
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>```
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!> \( Q = H(1)^\dagger H(2)^\dagger \cdots H(k)^\dagger \)
!> as returned by [[f77_gerqf:zgerqf]].
!> \( Q \) is of order \( M \) if `SIDE = 'L'`
!> and of order \( N \) if `SIDE = 'R'`.
pure subroutine zunmrq(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    import :: REAL64
    character, intent(in) :: side
    character, intent(in) :: trans
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for UNGQR.
!> Supports c, z.
!> See also: [[mfi_ungqr]], [[cungqr]], [[zungqr]].
interface f77_ungqr
!> Original interface for CUNGQR
!> See also: [[mfi_ungqr]], [[f77_ungqr]].
!> This routine generates an \(M \times N \) complex(REAL32)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:cgeqrf]].
pure subroutine cungqr(m, n, k, a, lda, tau, work, lwork, info)
    import :: REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
!> Original interface for ZUNGQR
!> See also: [[mfi_ungqr]], [[f77_ungqr]].
!> This routine generates an \(M \times N \) complex(REAL64)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( M \).
!> \( Q  =  H(1) H(2) . . . H(k) \)
!> as returned by [[f77_geqrf:zgeqrf]].
pure subroutine zungqr(m, n, k, a, lda, tau, work, lwork, info)
    import :: REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for UNGRQ.
!> Supports c, z.
!> See also: [[mfi_ungrq]], [[cungrq]], [[zungrq]].
interface f77_ungrq
!> Original interface for CUNGRQ
!> See also: [[mfi_ungrq]], [[f77_ungrq]].
!> This routine generates an \(M \times N \) complex(REAL32)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger . . . H(k)^\dagger \)
!> as returned by [[f77_gerqf:cgerqf]].
pure subroutine cungrq(m, n, k, a, lda, tau, work, lwork, info)
    import :: REAL32
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(out) :: work(*)
    complex(REAL32), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
!> Original interface for ZUNGRQ
!> See also: [[mfi_ungrq]], [[f77_ungrq]].
!> This routine generates an \(M \times N \) complex(REAL64)
!> matrix \( Q \) with orthonormal columns,
!> which is defined as the first \( N \) columns of a product of \( K \) elementary
!> reflectors of order \( N \).
!> \( Q  =  H(1)^\dagger H(2)^\dagger . . . H(k)^\dagger \)
!> as returned by [[f77_gerqf:zgerqf]].
pure subroutine zungrq(m, n, k, a, lda, tau, work, lwork, info)
    import :: REAL64
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(out) :: work(*)
    complex(REAL64), intent(in) :: tau(*)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(out) :: info
end subroutine
end interface
!> Generic old style interface for LARTG.
!> Supports s, d, c, z.
!> See also: [[mfi_lartg]], [[slartg]], [[dlartg]], [[clartg]], [[zlartg]].
interface f77_lartg
!> Original interface for SLARTG
!> See also: [[mfi_lartg]], [[f77_lartg]].
pure subroutine slartg(f, g, c, s, r)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: c
    real(REAL32), intent(inout) :: f
    real(REAL32), intent(inout) :: g
    real(REAL32), intent(inout) :: r
    real(REAL32), intent(inout) :: s
end subroutine
!> Original interface for DLARTG
!> See also: [[mfi_lartg]], [[f77_lartg]].
pure subroutine dlartg(f, g, c, s, r)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: c
    real(REAL64), intent(inout) :: f
    real(REAL64), intent(inout) :: g
    real(REAL64), intent(inout) :: r
    real(REAL64), intent(inout) :: s
end subroutine
!> Original interface for CLARTG
!> See also: [[mfi_lartg]], [[f77_lartg]].
pure subroutine clartg(f, g, c, s, r)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: c
    complex(REAL32), intent(inout) :: f
    complex(REAL32), intent(inout) :: g
    complex(REAL32), intent(inout) :: r
    complex(REAL32), intent(inout) :: s
end subroutine
!> Original interface for ZLARTG
!> See also: [[mfi_lartg]], [[f77_lartg]].
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

    interface f77_xerbla
        pure subroutine xerbla(name,info)
            character(*), intent(in) :: name
            integer, intent(in) :: info
        end subroutine
    end interface f77_xerbla

end module

