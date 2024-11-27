module f77_blas
use iso_fortran_env
implicit none

!FIXME rot, dot, rotg, nrm2: problem with functions that have TYPE /= TYPE_result
!https://spec.oneapi.com/versions/latest/elements/oneMKL/source/domains/blas/asum.html#onemkl-blas-asum
!FIXME sdsdot: Weird specific interface: computes a vec-vect dot product but perform a sum
!FIXME scal: problem with functions that have TYPE /= TYPE_scalar
!https://spec.oneapi.com/versions/latest/elements/oneMKL/source/domains/blas/scal.html#onemkl-blas-scal

! BLAS level 1
!!$:f77_interface('?asum',  DEFAULT_TYPES, asum, result=REAL_TYPES)
interface
pure subroutine saxpy(n, a, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(*)
    real(wp), intent(in) :: a
    real(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine daxpy(n, a, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(*)
    real(wp), intent(in) :: a
    real(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine caxpy(n, a, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: a
    complex(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zaxpy(n, a, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: a
    complex(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_axpy
    procedure :: saxpy
    procedure :: daxpy
    procedure :: caxpy
    procedure :: zaxpy
end interface
interface
pure subroutine scopy(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dcopy(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine ccopy(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zcopy(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_copy
    procedure :: scopy
    procedure :: dcopy
    procedure :: ccopy
    procedure :: zcopy
end interface
!$:f77_interface('?dot',  REAL_TYPES, dot_product, result=REAL_TYPES)
!$:f77_interface('sdsdot')
interface
pure function cdotu(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp) :: cdotu
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
pure function zdotu(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp) :: zdotu
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
end interface
interface f77_dotu
    procedure :: cdotu
    procedure :: zdotu
end interface
interface
pure function cdotc(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp) :: cdotc
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
pure function zdotc(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp) :: zdotc
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
end interface
interface f77_dotc
    procedure :: cdotc
    procedure :: zdotc
end interface
!$:f77_interface('?nrm2', DEFAULT_TYPES, nrm2, result=REAL_TYPES)
!$:f77_interface('?rot',  DEFAULT_TYPES, rot,  result=REAL_TYPES)
!$:f77_interface('?rotg', DEFAULT_TYPES, rotg, result=REAL_TYPES)
interface
pure subroutine srotm(n, x, incx, y, incy, param)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(inout) :: x(*)
    real(wp), intent(inout) :: y(*)
    real(wp), intent(in) :: param(5)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine drotm(n, x, incx, y, incy, param)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(inout) :: x(*)
    real(wp), intent(inout) :: y(*)
    real(wp), intent(in) :: param(5)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_rotm
    procedure :: srotm
    procedure :: drotm
end interface
interface
pure subroutine srotmg(d1, d2, x1, y1, param)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: y1
    real(wp), intent(out) :: param(5)
    real(wp), intent(inout) :: d1
    real(wp), intent(inout) :: d2
    real(wp), intent(inout) :: x1
end subroutine
pure subroutine drotmg(d1, d2, x1, y1, param)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: y1
    real(wp), intent(out) :: param(5)
    real(wp), intent(inout) :: d1
    real(wp), intent(inout) :: d2
    real(wp), intent(inout) :: x1
end subroutine
end interface
interface f77_rotmg
    procedure :: srotmg
    procedure :: drotmg
end interface
!$:f77_interface('?scal')
interface
pure subroutine sswap(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dswap(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine cswap(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zswap(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_swap
    procedure :: sswap
    procedure :: dswap
    procedure :: cswap
    procedure :: zswap
end interface

! BLAS level 2
interface
pure subroutine sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_gbmv
    procedure :: sgbmv
    procedure :: dgbmv
    procedure :: cgbmv
    procedure :: zgbmv
end interface
interface
pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_gemv
    procedure :: sgemv
    procedure :: dgemv
    procedure :: cgemv
    procedure :: zgemv
end interface
interface
pure subroutine sger(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(*)
    real(wp), intent(in) :: y(*)
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(*)
    real(wp), intent(in) :: y(*)
    real(wp), intent(inout) :: a(lda,*)
    real(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_ger
    procedure :: sger
    procedure :: dger
end interface
interface
pure subroutine cgerc(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_gerc
    procedure :: cgerc
    procedure :: zgerc
end interface
interface
pure subroutine cgeru(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zgeru(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    complex(wp), intent(inout) :: a(lda,*)
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_geru
    procedure :: cgeru
    procedure :: zgeru
end interface
interface
pure subroutine chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_hbmv
    procedure :: chbmv
    procedure :: zhbmv
end interface
interface
pure subroutine chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_hemv
    procedure :: chemv
    procedure :: zhemv
end interface
interface
pure subroutine cher(uplo, n, alpha, x, incx, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine zher(uplo, n, alpha, x, incx, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_her
    procedure :: cher
    procedure :: zher
end interface
interface
pure subroutine cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    complex(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    complex(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_her2
    procedure :: cher2
    procedure :: zher2
end interface
interface
pure subroutine chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: ap(*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: ap(*)
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_hpmv
    procedure :: chpmv
    procedure :: zhpmv
end interface
interface
pure subroutine chpr(uplo, n, alpha, x, incx, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine zhpr(uplo, n, alpha, x, incx, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_hpr
    procedure :: chpr
    procedure :: zhpr
end interface
interface
pure subroutine chpr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    complex(wp), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: x(*)
    complex(wp), intent(in) :: y(*)
    complex(wp), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_hpr2
    procedure :: chpr2
    procedure :: zhpr2
end interface
interface
pure subroutine ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_sbmv
    procedure :: ssbmv
    procedure :: dsbmv
end interface
interface
pure subroutine sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: ap(*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: ap(*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_spmv
    procedure :: sspmv
    procedure :: dspmv
end interface
interface
pure subroutine sspr(uplo, n, alpha, x, incx, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine dspr(uplo, n, alpha, x, incx, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_spr
    procedure :: sspr
    procedure :: dspr
end interface
interface
pure subroutine sspr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(*)
    real(wp), intent(in) :: y(*)
    real(wp), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dspr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(*)
    real(wp), intent(in) :: y(*)
    real(wp), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_spr2
    procedure :: sspr2
    procedure :: dspr2
end interface
interface
pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_symv
    procedure :: ssymv
    procedure :: dsymv
end interface
interface
pure subroutine ssyr(uplo, n, alpha, x, incx, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine dsyr(uplo, n, alpha, x, incx, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(*)
    real(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_syr
    procedure :: ssyr
    procedure :: dsyr
end interface
interface
pure subroutine ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: x(*)
    real(wp), intent(in) :: y(*)
    real(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: x(*)
    real(wp), intent(in) :: y(*)
    real(wp), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
interface f77_syr2
    procedure :: ssyr2
    procedure :: dsyr2
end interface
interface
pure subroutine stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine ctbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine ztbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_tbmv
    procedure :: stbmv
    procedure :: dtbmv
    procedure :: ctbmv
    procedure :: ztbmv
end interface
interface
pure subroutine stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine ctbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine ztbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_tbsv
    procedure :: stbsv
    procedure :: dtbsv
    procedure :: ctbsv
    procedure :: ztbsv
end interface
interface
pure subroutine stpmv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: ap(*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: ap(*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine ctpmv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: ap(*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine ztpmv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: ap(*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_tpmv
    procedure :: stpmv
    procedure :: dtpmv
    procedure :: ctpmv
    procedure :: ztpmv
end interface
interface
pure subroutine stpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: ap(*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine dtpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: ap(*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine ctpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: ap(*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine ztpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: ap(*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_tpsv
    procedure :: stpsv
    procedure :: dtpsv
    procedure :: ctpsv
    procedure :: ztpsv
end interface
interface
pure subroutine strmv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine ctrmv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine ztrmv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_trmv
    procedure :: strmv
    procedure :: dtrmv
    procedure :: ctrmv
    procedure :: ztrmv
end interface
interface
pure subroutine strsv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine ctrsv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
interface f77_trsv
    procedure :: strsv
    procedure :: dtrsv
    procedure :: ctrsv
    procedure :: ztrsv
end interface

! BLAS level 3
interface
pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: b(ldb,*)
    real(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: b(ldb,*)
    real(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: b(ldb,*)
    complex(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: b(ldb,*)
    complex(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
interface f77_gemm
    procedure :: sgemm
    procedure :: dgemm
    procedure :: cgemm
    procedure :: zgemm
end interface
interface
pure subroutine chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: b(ldb,*)
    complex(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: b(ldb,*)
    complex(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    complex(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
interface f77_hemm
    procedure :: chemm
    procedure :: zhemm
end interface
interface
pure subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
pure subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
end interface
interface f77_herk
    procedure :: cherk
    procedure :: zherk
end interface
interface
pure subroutine cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: b(ldb,*)
    complex(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(in) :: b(ldb,*)
    complex(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    complex(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
interface f77_her2k
    procedure :: cher2k
    procedure :: zher2k
end interface
interface
pure subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: b(ldb,*)
    real(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: b(ldb,*)
    real(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
interface f77_symm
    procedure :: ssymm
    procedure :: dsymm
end interface
interface
pure subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
pure subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
end interface
interface f77_syrk
    procedure :: ssyrk
    procedure :: dsyrk
end interface
interface
pure subroutine ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: b(ldb,*)
    real(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(in) :: b(ldb,*)
    real(wp), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
interface f77_syr2k
    procedure :: ssyr2k
    procedure :: dsyr2k
end interface
interface
pure subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
end interface
interface f77_trmm
    procedure :: strmm
    procedure :: dtrmm
    procedure :: ctrmm
    procedure :: ztrmm
end interface
interface
pure subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(wp), intent(in) :: a(lda,*)
    real(wp), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(wp), intent(in) :: a(lda,*)
    complex(wp), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(wp), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
end interface
interface f77_trsm
    procedure :: strsm
    procedure :: dtrsm
    procedure :: ctrsm
    procedure :: ztrsm
end interface


! Specific interfaces for slamch and dlamch
! as fortran can't differentiate them by the return kind
interface
    pure real(REAL32) function slamch(cmach)
        import :: REAL32
        character, intent(in) :: cmach
    end function

    pure real(REAL64) function dlamch(cmach)
        import :: REAL64
        character, intent(in) :: cmach
    end function
end interface

! Extensions
! BLAS Level 1 - Utils / Extensions
interface
pure function isamax(n, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    integer :: isamax
    real(wp), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
pure function idamax(n, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    integer :: idamax
    real(wp), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
pure function icamax(n, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    integer :: icamax
    complex(wp), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
pure function izamax(n, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    integer :: izamax
    complex(wp), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
end interface
interface f77_iamax
    procedure :: isamax
    procedure :: idamax
    procedure :: icamax
    procedure :: izamax
end interface
! Implement the blas extensions in
interface f77_iamin
    procedure :: isamin
    procedure :: idamin
    procedure :: icamin
    procedure :: izamin
end interface
contains
pure function isamin(n, x, incx)
    integer, parameter :: wp = REAL32
    integer :: isamin
    real(wp), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        isamin = 0
        return
    end if
    isamin = minloc(x(1:n:incx),dim=1)
end function
pure function idamin(n, x, incx)
    integer, parameter :: wp = REAL64
    integer :: idamin
    real(wp), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        idamin = 0
        return
    end if
    idamin = minloc(x(1:n:incx),dim=1)
end function
pure function icamin(n, x, incx)
    integer, parameter :: wp = REAL32
    integer :: icamin
    complex(wp), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        icamin = 0
        return
    end if
    icamin = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
end function
pure function izamin(n, x, incx)
    integer, parameter :: wp = REAL64
    integer :: izamin
    complex(wp), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        izamin = 0
        return
    end if
    izamin = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
end function

end module

