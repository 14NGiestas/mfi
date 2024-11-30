!> Improved and original F77 interfaces for BLAS
module f77_blas
use iso_fortran_env
implicit none

!> ?copy supports s, d, c, z.
!> See also: [[mfi_copy]], [[f77_copy]].
interface
pure subroutine scopy(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dcopy(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine ccopy(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zcopy(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?swap supports s, d, c, z.
!> See also: [[mfi_swap]], [[f77_swap]].
interface
pure subroutine sswap(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dswap(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine cswap(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zswap(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?axpy supports s, d, c, z.
!> See also: [[mfi_axpy]], [[f77_axpy]].
interface
pure subroutine saxpy(n, a, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: a
    real(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine daxpy(n, a, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: a
    real(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine caxpy(n, a, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: a
    complex(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zaxpy(n, a, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: a
    complex(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?dot supports s, d.
!> See also: [[mfi_dot]], [[f77_dot]].
interface
pure function sdot(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32) :: sdot
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
pure function ddot(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64) :: ddot
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
end interface
!> ?dotc supports c, z.
!> See also: [[mfi_dotc]], [[f77_dotc]].
interface
pure function cdotc(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32) :: cdotc
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
pure function zdotc(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64) :: zdotc
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
end interface
!> ?dotu supports c, z.
!> See also: [[mfi_dotu]], [[f77_dotu]].
interface
pure function cdotu(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32) :: cdotu
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
pure function zdotu(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64) :: zdotu
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end function
end interface
!> ?asum supports s, d, sc, dz.
!> See also: [[mfi_asum]], [[f77_asum]].
interface
pure function sasum(n, x, incx)
    import :: REAL32
    real(REAL32) :: sasum
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
pure function dasum(n, x, incx)
    import :: REAL64
    real(REAL64) :: dasum
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
pure function scasum(n, x, incx)
    import :: REAL32
    real(REAL32) :: scasum
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
pure function dzasum(n, x, incx)
    import :: REAL64
    real(REAL64) :: dzasum
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
end interface
!> ?nrm2 supports s, d, sc, dz.
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
interface
pure function snrm2(n, x, incx)
    import :: REAL32
    real(REAL32) :: snrm2
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
pure function dnrm2(n, x, incx)
    import :: REAL64
    real(REAL64) :: dnrm2
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
pure function scnrm2(n, x, incx)
    import :: REAL32
    real(REAL32) :: scnrm2
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
pure function dznrm2(n, x, incx)
    import :: REAL64
    real(REAL64) :: dznrm2
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
end interface
!> ?rot supports s, d, c, z, cs, zd.
!> See also: [[mfi_rot]], [[f77_rot]].
interface
!> SROT applies a plane rotation.
!> ['s']
pure subroutine srot(n, x, incx, y, incy, c, s)
    import :: REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL32), intent(in) :: c
    real(REAL32), intent(in) :: s
end subroutine
!> DROT applies a plane rotation.
!> ['d']
pure subroutine drot(n, x, incx, y, incy, c, s)
    import :: REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL64), intent(in) :: c
    real(REAL64), intent(in) :: s
end subroutine
!> CROT applies a plane rotation.
!> ['c']
pure subroutine crot(n, x, incx, y, incy, c, s)
    import :: REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL32), intent(in) :: c
    complex(REAL32), intent(in) :: s
end subroutine
!> ZROT applies a plane rotation.
!> ['z']
pure subroutine zrot(n, x, incx, y, incy, c, s)
    import :: REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL64), intent(in) :: c
    complex(REAL64), intent(in) :: s
end subroutine
!> CSROT applies a plane rotation.
!> ['c', 's']
pure subroutine csrot(n, x, incx, y, incy, c, s)
    import :: REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL32), intent(in) :: c
    real(REAL32), intent(in) :: s
end subroutine
!> ZDROT applies a plane rotation.
!> ['z', 'd']
pure subroutine zdrot(n, x, incx, y, incy, c, s)
    import :: REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
    real(REAL64), intent(in) :: c
    real(REAL64), intent(in) :: s
end subroutine
end interface
!> ?rotg supports s, d, c, z.
!> See also: [[mfi_rotg]], [[f77_rotg]].
interface
!>srotg generates a Givens rotation with real cosine and complex sine:
!>```
!> [  c  s ] [ a ] = [ r ]
!> [ -s  c ] [ b ]   [ 0 ]
!>```
!> satisfying `c**2 + s**2 = 1`.
pure subroutine srotg(a, b, c, s)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: a
    real(REAL32), intent(inout) :: b
    real(REAL32), intent(out) :: c
    real(REAL32), intent(out) :: s
end subroutine

!>drotg generates a Givens rotation with real cosine and complex sine:
!>```
!> [  c  s ] [ a ] = [ r ]
!> [ -s  c ] [ b ]   [ 0 ]
!>```
!> satisfying `c**2 + s**2 = 1`.
pure subroutine drotg(a, b, c, s)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: a
    real(REAL64), intent(inout) :: b
    real(REAL64), intent(out) :: c
    real(REAL64), intent(out) :: s
end subroutine

!>crotg generates a Givens rotation with real cosine and complex sine:
!>```
!>  [  c         s ] [ a ] = [ r ]
!>  [ -conjg(s)  c ] [ b ]   [ 0 ]
!>```
!> where c is real, s is complex, and `c**2 + conjg(s)*s = 1`.
pure subroutine crotg(a, b, c, s)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(inout) :: a
    complex(REAL32), intent(inout) :: b
    real(REAL32), intent(out) :: c
    complex(REAL32), intent(out) :: s
end subroutine

!>zrotg generates a Givens rotation with real cosine and complex sine:
!>```
!>  [  c         s ] [ a ] = [ r ]
!>  [ -conjg(s)  c ] [ b ]   [ 0 ]
!>```
!> where c is real, s is complex, and `c**2 + conjg(s)*s = 1`.
pure subroutine zrotg(a, b, c, s)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(inout) :: a
    complex(REAL64), intent(inout) :: b
    real(REAL64), intent(out) :: c
    complex(REAL64), intent(out) :: s
end subroutine

end interface
!> ?rotm supports s, d.
!> See also: [[mfi_rotm]], [[f77_rotm]].
interface
pure subroutine srotm(n, x, incx, y, incy, param)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(inout) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    real(REAL32), intent(in) :: param(5)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine drotm(n, x, incx, y, incy, param)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(inout) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    real(REAL64), intent(in) :: param(5)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?rotmg supports s, d.
!> See also: [[mfi_rotmg]], [[f77_rotmg]].
interface
pure subroutine srotmg(d1, d2, x1, y1, param)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: y1
    real(REAL32), intent(out) :: param(5)
    real(REAL32), intent(inout) :: d1
    real(REAL32), intent(inout) :: d2
    real(REAL32), intent(inout) :: x1
end subroutine
pure subroutine drotmg(d1, d2, x1, y1, param)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: y1
    real(REAL64), intent(out) :: param(5)
    real(REAL64), intent(inout) :: d1
    real(REAL64), intent(inout) :: d2
    real(REAL64), intent(inout) :: x1
end subroutine
end interface
!> ?scal supports s, d, c, z, cs, zd.
!> See also: [[mfi_scal]], [[f77_scal]].
interface
!> SSCAL scales a vector by a constant.
pure subroutine sscal(n, a, x, incx)
    import :: REAL32
    real(REAL32), intent(inout) :: x(*)
    real(REAL32), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> DSCAL scales a vector by a constant.
pure subroutine dscal(n, a, x, incx)
    import :: REAL64
    real(REAL64), intent(inout) :: x(*)
    real(REAL64), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> CSCAL scales a vector by a constant.
pure subroutine cscal(n, a, x, incx)
    import :: REAL32
    complex(REAL32), intent(inout) :: x(*)
    complex(REAL32), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> ZSCAL scales a vector by a constant.
pure subroutine zscal(n, a, x, incx)
    import :: REAL64
    complex(REAL64), intent(inout) :: x(*)
    complex(REAL64), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> CSSCAL scales a vector by a constant.
pure subroutine csscal(n, a, x, incx)
    import :: REAL32
    complex(REAL32), intent(inout) :: x(*)
    real(REAL32), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> ZDSCAL scales a vector by a constant.
pure subroutine zdscal(n, a, x, incx)
    import :: REAL64
    complex(REAL64), intent(inout) :: x(*)
    real(REAL64), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
!> ?gbmv supports s, d, c, z.
!> See also: [[mfi_gbmv]], [[f77_gbmv]].
interface
pure subroutine sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
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
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
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
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
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
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: kl
    integer, intent(in) :: ku
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?gemv supports s, d, c, z.
!> See also: [[mfi_gemv]], [[f77_gemv]].
interface
pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    character, intent(in) :: trans
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    character, intent(in) :: trans
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?ger supports s, d.
!> See also: [[mfi_ger]], [[f77_ger]].
interface
pure subroutine sger(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    real(REAL32), intent(inout) :: a(lda,*)
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    real(REAL64), intent(inout) :: a(lda,*)
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?gerc supports c, z.
!> See also: [[mfi_gerc]], [[f77_gerc]].
interface
pure subroutine cgerc(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?geru supports c, z.
!> See also: [[mfi_geru]], [[f77_geru]].
interface
pure subroutine cgeru(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    complex(REAL32), intent(inout) :: a(lda,*)
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zgeru(m, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    complex(REAL64), intent(inout) :: a(lda,*)
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?hbmv supports c, z.
!> See also: [[mfi_hbmv]], [[f77_hbmv]].
interface
pure subroutine chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?hemv supports c, z.
!> See also: [[mfi_hemv]], [[f77_hemv]].
interface
pure subroutine chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?her supports c, z.
!> See also: [[mfi_her]], [[f77_her]].
interface
pure subroutine cher(uplo, n, alpha, x, incx, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine zher(uplo, n, alpha, x, incx, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
!> ?her2 supports c, z.
!> See also: [[mfi_her2]], [[f77_her2]].
interface
pure subroutine cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    complex(REAL32), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    complex(REAL64), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?hpmv supports c, z.
!> See also: [[mfi_hpmv]], [[f77_hpmv]].
interface
pure subroutine chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: ap(*)
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: ap(*)
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?hpr supports c, z.
!> See also: [[mfi_hpr]], [[f77_hpr]].
interface
pure subroutine chpr(uplo, n, alpha, x, incx, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine zhpr(uplo, n, alpha, x, incx, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(wp), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
!> ?hpr2 supports c, z.
!> See also: [[mfi_hpr2]], [[f77_hpr2]].
interface
pure subroutine chpr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(in) :: y(*)
    complex(REAL32), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: x(*)
    complex(REAL64), intent(in) :: y(*)
    complex(REAL64), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?sbmv supports s, d.
!> See also: [[mfi_sbmv]], [[f77_sbmv]].
interface
pure subroutine ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?spmv supports s, d.
!> See also: [[mfi_spmv]], [[f77_spmv]].
interface
pure subroutine sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: ap(*)
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: ap(*)
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?spr supports s, d.
!> See also: [[mfi_spr]], [[f77_spr]].
interface
pure subroutine sspr(uplo, n, alpha, x, incx, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine dspr(uplo, n, alpha, x, incx, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
!> ?spr2 supports s, d.
!> See also: [[mfi_spr2]], [[f77_spr2]].
interface
pure subroutine sspr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    real(REAL32), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dspr2(uplo, n, alpha, x, incx, y, incy, ap)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    real(REAL64), intent(inout) :: ap(*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?symv supports s, d.
!> See also: [[mfi_symv]], [[f77_symv]].
interface
pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?syr supports s, d.
!> See also: [[mfi_syr]], [[f77_syr]].
interface
pure subroutine ssyr(uplo, n, alpha, x, incx, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
pure subroutine dsyr(uplo, n, alpha, x, incx, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
!> ?syr2 supports s, d.
!> See also: [[mfi_syr2]], [[f77_syr2]].
interface
pure subroutine ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(in) :: y(*)
    real(REAL32), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
pure subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(in) :: y(*)
    real(REAL64), intent(inout) :: a(lda,*)
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
end interface
!> ?tbmv supports s, d, c, z.
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
interface
pure subroutine stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: x(*)
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
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: x(*)
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
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: x(*)
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
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
!> ?tbsv supports s, d, c, z.
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
interface
pure subroutine stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: x(*)
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
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: x(*)
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
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: x(*)
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
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
!> ?tpmv supports s, d, c, z.
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
interface
pure subroutine stpmv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: ap(*)
    real(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: ap(*)
    real(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine ctpmv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: ap(*)
    complex(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine ztpmv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: ap(*)
    complex(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
!> ?tpsv supports s, d, c, z.
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
interface
pure subroutine stpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: ap(*)
    real(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine dtpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: ap(*)
    real(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine ctpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: ap(*)
    complex(REAL32), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
pure subroutine ztpsv(uplo, trans, diag, n, ap, x, incx)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: ap(*)
    complex(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
!> ?trmv supports s, d, c, z.
!> See also: [[mfi_trmv]], [[f77_trmv]].
interface
pure subroutine strmv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: x(*)
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
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: x(*)
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
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: x(*)
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
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
!> ?trsv supports s, d, c, z.
!> See also: [[mfi_trsv]], [[f77_trsv]].
interface
pure subroutine strsv(uplo, trans, diag, n, a, lda, x, incx)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: x(*)
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
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: x(*)
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
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: x(*)
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
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: x(*)
    character, intent(in) :: uplo
    character, intent(in) :: trans
    character, intent(in) :: diag
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: incx
end subroutine
end interface
!> ?gemm supports s, d, c, z.
!> See also: [[mfi_gemm]], [[f77_gemm]].
interface
pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: b(ldb,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
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
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: b(ldb,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
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
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: b(ldb,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
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
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: b(ldb,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: transa
    character, intent(in) :: transb
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
!> ?hemm supports c, z.
!> See also: [[mfi_hemm]], [[f77_hemm]].
interface
pure subroutine chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: b(ldb,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
    complex(REAL32), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: b(ldb,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    complex(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
!> ?herk supports c, z.
!> See also: [[mfi_herk]], [[f77_herk]].
interface
pure subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
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
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
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
!> ?her2k supports c, z.
!> See also: [[mfi_her2k]], [[f77_her2k]].
interface
pure subroutine cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(in) :: b(ldb,*)
    complex(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    complex(REAL32), intent(in) :: alpha
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
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(in) :: b(ldb,*)
    complex(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    complex(REAL64), intent(in) :: alpha
    real(wp), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
!> ?symm supports s, d.
!> See also: [[mfi_symm]], [[f77_symm]].
interface
pure subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: b(ldb,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: b(ldb,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
!> ?syrk supports s, d.
!> See also: [[mfi_syrk]], [[f77_syrk]].
interface
pure subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
pure subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldc
end subroutine
end interface
!> ?syr2k supports s, d.
!> See also: [[mfi_syr2k]], [[f77_syr2k]].
interface
pure subroutine ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(in) :: b(ldb,*)
    real(REAL32), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(REAL32), intent(in) :: alpha
    real(REAL32), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
pure subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(in) :: b(ldb,*)
    real(REAL64), intent(inout) :: c(ldc,*)
    character, intent(in) :: trans
    character, intent(in) :: uplo
    real(REAL64), intent(in) :: alpha
    real(REAL64), intent(in) :: beta
    integer, intent(in) :: n
    integer, intent(in) :: k
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
    integer, intent(in) :: ldc
end subroutine
end interface
!> ?trmm supports s, d, c, z.
!> See also: [[mfi_trmm]], [[f77_trmm]].
interface
pure subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
end interface
!> ?trsm supports s, d, c, z.
!> See also: [[mfi_trsm]], [[f77_trsm]].
interface
pure subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: a(lda,*)
    real(REAL32), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: a(lda,*)
    real(REAL64), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    real(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: a(lda,*)
    complex(REAL32), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(REAL32), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
pure subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    import :: REAL64
    integer, parameter :: wp = REAL64
    complex(REAL64), intent(in) :: a(lda,*)
    complex(REAL64), intent(inout) :: b(ldb,*)
    character, intent(in) :: side
    character, intent(in) :: uplo
    character, intent(in) :: transa
    character, intent(in) :: diag
    complex(REAL64), intent(in) :: alpha
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: lda
    integer, intent(in) :: ldb
end subroutine
end interface

interface f77_copy
    procedure :: scopy
    procedure :: dcopy
    procedure :: ccopy
    procedure :: zcopy
end interface
interface f77_swap
    procedure :: sswap
    procedure :: dswap
    procedure :: cswap
    procedure :: zswap
end interface
interface f77_axpy
    procedure :: saxpy
    procedure :: daxpy
    procedure :: caxpy
    procedure :: zaxpy
end interface
interface f77_dot
    procedure :: sdot
    procedure :: ddot
end interface
interface f77_dotc
    procedure :: cdotc
    procedure :: zdotc
end interface
interface f77_dotu
    procedure :: cdotu
    procedure :: zdotu
end interface
interface f77_asum
    procedure :: sasum
    procedure :: dasum
    procedure :: scasum
    procedure :: dzasum
end interface
interface f77_nrm2
    procedure :: snrm2
    procedure :: dnrm2
    procedure :: scnrm2
    procedure :: dznrm2
end interface
interface f77_rot
    procedure :: srot
    procedure :: drot
    procedure :: crot
    procedure :: zrot
    procedure :: csrot
    procedure :: zdrot
end interface
interface f77_rotg
    procedure :: srotg
    procedure :: drotg
    procedure :: crotg
    procedure :: zrotg
end interface
interface f77_rotm
    procedure :: srotm
    procedure :: drotm
end interface
interface f77_rotmg
    procedure :: srotmg
    procedure :: drotmg
end interface
interface f77_scal
    procedure :: sscal
    procedure :: dscal
    procedure :: cscal
    procedure :: zscal
    procedure :: csscal
    procedure :: zdscal
end interface
interface f77_gbmv
    procedure :: sgbmv
    procedure :: dgbmv
    procedure :: cgbmv
    procedure :: zgbmv
end interface
interface f77_gemv
    procedure :: sgemv
    procedure :: dgemv
    procedure :: cgemv
    procedure :: zgemv
end interface
interface f77_ger
    procedure :: sger
    procedure :: dger
end interface
interface f77_gerc
    procedure :: cgerc
    procedure :: zgerc
end interface
interface f77_geru
    procedure :: cgeru
    procedure :: zgeru
end interface
interface f77_hbmv
    procedure :: chbmv
    procedure :: zhbmv
end interface
interface f77_hemv
    procedure :: chemv
    procedure :: zhemv
end interface
interface f77_her
    procedure :: cher
    procedure :: zher
end interface
interface f77_her2
    procedure :: cher2
    procedure :: zher2
end interface
interface f77_hpmv
    procedure :: chpmv
    procedure :: zhpmv
end interface
interface f77_hpr
    procedure :: chpr
    procedure :: zhpr
end interface
interface f77_hpr2
    procedure :: chpr2
    procedure :: zhpr2
end interface
interface f77_sbmv
    procedure :: ssbmv
    procedure :: dsbmv
end interface
interface f77_spmv
    procedure :: sspmv
    procedure :: dspmv
end interface
interface f77_spr
    procedure :: sspr
    procedure :: dspr
end interface
interface f77_spr2
    procedure :: sspr2
    procedure :: dspr2
end interface
interface f77_symv
    procedure :: ssymv
    procedure :: dsymv
end interface
interface f77_syr
    procedure :: ssyr
    procedure :: dsyr
end interface
interface f77_syr2
    procedure :: ssyr2
    procedure :: dsyr2
end interface
interface f77_tbmv
    procedure :: stbmv
    procedure :: dtbmv
    procedure :: ctbmv
    procedure :: ztbmv
end interface
interface f77_tbsv
    procedure :: stbsv
    procedure :: dtbsv
    procedure :: ctbsv
    procedure :: ztbsv
end interface
interface f77_tpmv
    procedure :: stpmv
    procedure :: dtpmv
    procedure :: ctpmv
    procedure :: ztpmv
end interface
interface f77_tpsv
    procedure :: stpsv
    procedure :: dtpsv
    procedure :: ctpsv
    procedure :: ztpsv
end interface
interface f77_trmv
    procedure :: strmv
    procedure :: dtrmv
    procedure :: ctrmv
    procedure :: ztrmv
end interface
interface f77_trsv
    procedure :: strsv
    procedure :: dtrsv
    procedure :: ctrsv
    procedure :: ztrsv
end interface
interface f77_gemm
    procedure :: sgemm
    procedure :: dgemm
    procedure :: cgemm
    procedure :: zgemm
end interface
interface f77_hemm
    procedure :: chemm
    procedure :: zhemm
end interface
interface f77_herk
    procedure :: cherk
    procedure :: zherk
end interface
interface f77_her2k
    procedure :: cher2k
    procedure :: zher2k
end interface
interface f77_symm
    procedure :: ssymm
    procedure :: dsymm
end interface
interface f77_syrk
    procedure :: ssyrk
    procedure :: dsyrk
end interface
interface f77_syr2k
    procedure :: ssyr2k
    procedure :: dsyr2k
end interface
interface f77_trmm
    procedure :: strmm
    procedure :: dtrmm
    procedure :: ctrmm
    procedure :: ztrmm
end interface
interface f77_trsm
    procedure :: strsm
    procedure :: dtrsm
    procedure :: ctrsm
    procedure :: ztrsm
end interface

!> ?lamch supports s, d. See [[mfi_lamch]] for the modern version.
interface
    !> SLAMCH determines single precision machine parameters.
    pure real(REAL32) function slamch(cmach)
        import :: REAL32
        character, intent(in) :: cmach
    end function

    !> DLAMCH determines double precision machine parameters.
    pure real(REAL64) function dlamch(cmach)
        import :: REAL64
        character, intent(in) :: cmach
    end function
end interface

interface
    !> Compute the inner product of two vectors with extended
    !> precision accumulation.
    !>
    !> Returns S.P. result with dot product accumulated in D.P.
    !> SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
    !> where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
    !> defined in a similar way using INCY.
    pure function sdsdot(n, sb, sx, incx, sy, incy)
        import :: REAL32
        integer, parameter :: wp = REAL32
        real(wp) :: sdsdot
        real(wp), intent(in) :: sx(*)
        real(wp), intent(in) :: sy(*)
        real(wp), intent(in) :: sb
        integer, intent(in) :: n
        integer, intent(in) :: incx
        integer, intent(in) :: incy
    end function

    !> Compute the inner product of two vectors with extended
    !> precision accumulation and result.
    !>
    !> Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
    !> DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
    !> where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
    !> defined in a similar way using INCY.
    pure function dsdot(n, sx, incx, sy, incy)
        import :: REAL32, REAL64
        integer, parameter :: sp = REAL32
        integer, parameter :: dp = REAL64
        real(dp) :: dsdot
        real(sp), intent(in) :: sx(*)
        real(sp), intent(in) :: sy(*)
        integer,  intent(in) :: n
        integer,  intent(in) :: incx
        integer,  intent(in) :: incy
    end function
end interface

! Extensions
! BLAS Level 1 - Utils / Extensions
! Implement the blas extensions in
interface f77_iamax
    procedure :: isamax
    procedure :: idamax
    procedure :: icamax
    procedure :: izamax
end interface
interface f77_iamin
    procedure :: isamin
    procedure :: idamin
    procedure :: icamin
    procedure :: izamin
end interface
contains
pure function isamax(n, x, incx)
    integer, parameter :: wp = REAL32
    integer :: isamax
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        isamax = 0
        return
    end if
    isamax = minloc(x(1:n:incx),dim=1)
end function
pure function idamax(n, x, incx)
    integer, parameter :: wp = REAL64
    integer :: idamax
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        idamax = 0
        return
    end if
    idamax = minloc(x(1:n:incx),dim=1)
end function
pure function icamax(n, x, incx)
    integer, parameter :: wp = REAL32
    integer :: icamax
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        icamax = 0
        return
    end if
    icamax = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
end function
pure function izamax(n, x, incx)
    integer, parameter :: wp = REAL64
    integer :: izamax
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    !If either n or incx are not positive, the routine returns 0.
    if (n <= 0 .or. incx <= 0) then
        izamax = 0
        return
    end if
    izamax = minloc(abs(real(x(1:n:incx))) + abs(aimag(x(1:n:incx))),dim=1)
end function
pure function isamin(n, x, incx)
    integer, parameter :: wp = REAL32
    integer :: isamin
    real(REAL32), intent(in) :: x(*)
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
    real(REAL64), intent(in) :: x(*)
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
    complex(REAL32), intent(in) :: x(*)
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
    complex(REAL64), intent(in) :: x(*)
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

