module f77_lapack
use iso_fortran_env
implicit none

interface f77_geqrf
pure subroutine sgeqrf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    real(wp), intent(inout) :: work(*)
end subroutine
pure subroutine dgeqrf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    real(wp), intent(inout) :: work(*)
end subroutine
pure subroutine cgeqrf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    complex(wp), intent(inout) :: work(*)
end subroutine
pure subroutine zgeqrf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    complex(wp), intent(inout) :: work(*)
end subroutine
end interface
interface f77_gerqf
pure subroutine sgerqf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    real(wp), intent(inout) :: work(*)
end subroutine
pure subroutine dgerqf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    real(wp), intent(inout) :: work(*)
end subroutine
pure subroutine cgerqf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    complex(wp), intent(inout) :: work(*)
end subroutine
pure subroutine zgerqf(m,n,a,lda,tau,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(out) :: tau(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    complex(wp), intent(inout) :: work(*)
end subroutine
end interface
interface f77_getrf
pure subroutine sgetrf(m,n,a,lda,ipiv,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
end subroutine
pure subroutine dgetrf(m,n,a,lda,ipiv,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
end subroutine
pure subroutine cgetrf(m,n,a,lda,ipiv,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
end subroutine
pure subroutine zgetrf(m,n,a,lda,ipiv,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(lda,*)
    integer, intent(out) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
end subroutine
end interface
interface f77_getri
pure subroutine sgetri(n,a,lda,ipiv,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(inout) :: work(*)
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
pure subroutine dgetri(n,a,lda,ipiv,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(inout) :: work(*)
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
pure subroutine cgetri(n,a,lda,ipiv,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(inout) :: work(*)
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
pure subroutine zgetri(n,a,lda,ipiv,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(inout) :: work(*)
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
end interface
interface f77_getrs
pure subroutine sgetrs(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(inout) :: b(ldb,*)
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
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(inout) :: b(ldb,*)
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
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(inout) :: b(ldb,*)
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
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(inout) :: b(ldb,*)
    character, intent(in) :: trans
    integer, intent(in) :: ipiv(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
end interface
interface f77_hetrf
pure subroutine chetrf(uplo, n, a, lda, ipiv, work, lwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: ipiv(*)
    complex(wp), intent(inout) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
pure subroutine zhetrf(uplo, n, a, lda, ipiv, work, lwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: ipiv(*)
    complex(wp), intent(inout) :: work(*)
    integer, intent(out) :: info
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
end subroutine
end interface
interface f77_hegv
pure subroutine chegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(inout) :: b(ldb,*)
    real(wp), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: itype
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
    complex(wp), intent(inout) :: work(*)
    real(wp), intent(in) :: rwork(*)
end subroutine
pure subroutine zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(inout) :: b(ldb,*)
    real(wp), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: itype
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: lwork
    complex(wp), intent(inout) :: work(*)
    real(wp), intent(in) :: rwork(*)
end subroutine
end interface
interface f77_heevd
pure subroutine cheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(in) :: lrwork
    integer, intent(in) :: liwork
    complex(wp), intent(inout) :: work(*)
    real(wp), intent(inout) :: rwork(*)
    integer, intent(inout) :: iwork(*)
end subroutine
pure subroutine zheevd(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: w(*)
    integer, intent(out) :: info
    character, intent(in) :: jobz
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: lwork
    integer, intent(in) :: lrwork
    integer, intent(in) :: liwork
    complex(wp), intent(inout) :: work(*)
    real(wp), intent(inout) :: rwork(*)
    integer, intent(inout) :: iwork(*)
end subroutine
end interface
interface f77_gesvd
pure subroutine sgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: s(*)
    real(wp), intent(out) :: u(ldu,*)
    real(wp), intent(out) :: vt(ldvt,*)
    integer, intent(out) :: info
    character, intent(in) :: jobu
    character, intent(in) :: jobvt
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldu
    integer, intent(in) :: ldvt
    integer, intent(in) :: lwork
    real(wp), intent(inout) :: work(*)
end subroutine
pure subroutine dgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: s(*)
    real(wp), intent(out) :: u(ldu,*)
    real(wp), intent(out) :: vt(ldvt,*)
    integer, intent(out) :: info
    character, intent(in) :: jobu
    character, intent(in) :: jobvt
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldu
    integer, intent(in) :: ldvt
    integer, intent(in) :: lwork
    real(wp), intent(inout) :: work(*)
end subroutine
pure subroutine cgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: s(*)
    complex(wp), intent(out) :: u(ldu,*)
    complex(wp), intent(out) :: vt(ldvt,*)
    integer, intent(out) :: info
    character, intent(in) :: jobu
    character, intent(in) :: jobvt
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldu
    integer, intent(in) :: ldvt
    integer, intent(in) :: lwork
    complex(wp), intent(inout) :: work(*)
    real(wp), intent(in) :: rwork(*)
end subroutine
pure subroutine zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(inout) :: a(lda,*)
    real(wp), intent(out) :: s(*)
    complex(wp), intent(out) :: u(ldu,*)
    complex(wp), intent(out) :: vt(ldvt,*)
    integer, intent(out) :: info
    character, intent(in) :: jobu
    character, intent(in) :: jobvt
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldu
    integer, intent(in) :: ldvt
    integer, intent(in) :: lwork
    complex(wp), intent(inout) :: work(*)
    real(wp), intent(in) :: rwork(*)
end subroutine
end interface
interface f77_potrf
pure subroutine spotrf(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine dpotrf(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine cpotrf(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine zpotrf(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
end interface
interface f77_potri
pure subroutine spotri(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine dpotri(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine cpotri(uplo, n, a, lda, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
pure subroutine zpotri(uplo, n, a, lda, info)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(out) :: info
end subroutine
end interface
interface f77_potrs
pure subroutine spotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: b(ldb,*)
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
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: b(ldb,*)
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
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: b(ldb,*)
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
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: b(ldb,*)
    character, intent(in) :: uplo
    integer, intent(in) :: n
    integer, intent(in) :: nrhs
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(out) :: info
end subroutine
end interface

! Other Auxiliary Routines
interface f77_lartg
pure subroutine slartg(f, g, c, s, r)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: c
    real(wp), intent(inout) :: f
    real(wp), intent(inout) :: g
    real(wp), intent(inout) :: r
    real(wp), intent(inout) :: s
end subroutine
pure subroutine dlartg(f, g, c, s, r)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: c
    real(wp), intent(inout) :: f
    real(wp), intent(inout) :: g
    real(wp), intent(inout) :: r
    real(wp), intent(inout) :: s
end subroutine
pure subroutine clartg(f, g, c, s, r)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: c
    complex(wp), intent(inout) :: f
    complex(wp), intent(inout) :: g
    complex(wp), intent(inout) :: r
    complex(wp), intent(inout) :: s
end subroutine
pure subroutine zlartg(f, g, c, s, r)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: c
    complex(wp), intent(inout) :: f
    complex(wp), intent(inout) :: g
    complex(wp), intent(inout) :: r
    complex(wp), intent(inout) :: s
end subroutine
end interface

    interface f77_xerbla
        pure subroutine xerbla(name,info)
            character(*), intent(in) :: name
            integer, intent(in) :: info
        end subroutine
    end interface f77_xerbla

end module

