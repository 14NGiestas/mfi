!> Improved and original F77 interfaces for BLAS
module f77_blas
use iso_fortran_env
implicit none

!> Generic old style interface for COPY.
!> Supports s, d, c, z.
!> See also: [[mfi_copy]], [[scopy]],[[dcopy]],[[ccopy]],[[zcopy]].
interface f77_copy
!> Original interface for SCOPY
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine scopy(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DCOPY
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine dcopy(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for CCOPY
!> See also: [[mfi_copy]], [[f77_copy]].
pure subroutine ccopy(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZCOPY
!> See also: [[mfi_copy]], [[f77_copy]].
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
!> Generic old style interface for SWAP.
!> Supports s, d, c, z.
!> See also: [[mfi_swap]], [[sswap]],[[dswap]],[[cswap]],[[zswap]].
interface f77_swap
!> Original interface for SSWAP
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine sswap(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: x(*)
    real(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for DSWAP
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine dswap(n, x, incx, y, incy)
    import :: REAL64
    integer, parameter :: wp = REAL64
    real(REAL64), intent(in) :: x(*)
    real(REAL64), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for CSWAP
!> See also: [[mfi_swap]], [[f77_swap]].
pure subroutine cswap(n, x, incx, y, incy)
    import :: REAL32
    integer, parameter :: wp = REAL32
    complex(REAL32), intent(in) :: x(*)
    complex(REAL32), intent(inout) :: y(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
    integer, intent(in) :: incy
end subroutine
!> Original interface for ZSWAP
!> See also: [[mfi_swap]], [[f77_swap]].
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
!> Generic old style interface for AXPY.
!> Supports s, d, c, z.
!> See also: [[mfi_axpy]], [[saxpy]],[[daxpy]],[[caxpy]],[[zaxpy]].
interface f77_axpy
!> Original interface for SAXPY
!> See also: [[mfi_axpy]], [[f77_axpy]].
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
!> Original interface for DAXPY
!> See also: [[mfi_axpy]], [[f77_axpy]].
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
!> Original interface for CAXPY
!> See also: [[mfi_axpy]], [[f77_axpy]].
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
!> Original interface for ZAXPY
!> See also: [[mfi_axpy]], [[f77_axpy]].
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
!> Generic old style interface for DOT.
!> Supports s, d.
!> See also: [[mfi_dot]], [[sdot]],[[ddot]].
interface f77_dot
!> Original interface for SDOT
!> See also: [[mfi_dot]], [[f77_dot]].
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
!> Original interface for DDOT
!> See also: [[mfi_dot]], [[f77_dot]].
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
!> Generic old style interface for DOTC.
!> Supports c, z.
!> See also: [[mfi_dotc]], [[cdotc]],[[zdotc]].
interface f77_dotc
!> Original interface for CDOTC
!> See also: [[mfi_dotc]], [[f77_dotc]].
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
!> Original interface for ZDOTC
!> See also: [[mfi_dotc]], [[f77_dotc]].
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
!> Generic old style interface for DOTU.
!> Supports c, z.
!> See also: [[mfi_dotu]], [[cdotu]],[[zdotu]].
interface f77_dotu
!> Original interface for CDOTU
!> See also: [[mfi_dotu]], [[f77_dotu]].
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
!> Original interface for ZDOTU
!> See also: [[mfi_dotu]], [[f77_dotu]].
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
!> Generic old style interface for ASUM.
!> Supports s, d, sc, dz.
!> See also: [[mfi_asum]], [[sasum]],[[dasum]],[[scasum]],[[dzasum]].
interface f77_asum
!> Original interface for SASUM
!> See also: [[mfi_asum]], [[f77_asum]].
pure function sasum(n, x, incx)
    import :: REAL32
    real(REAL32) :: sasum
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for DASUM
!> See also: [[mfi_asum]], [[f77_asum]].
pure function dasum(n, x, incx)
    import :: REAL64
    real(REAL64) :: dasum
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for SCASUM
!> See also: [[mfi_asum]], [[f77_asum]].
pure function scasum(n, x, incx)
    import :: REAL32
    real(REAL32) :: scasum
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for DZASUM
!> See also: [[mfi_asum]], [[f77_asum]].
pure function dzasum(n, x, incx)
    import :: REAL64
    real(REAL64) :: dzasum
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
end interface
!> Generic old style interface for NRM2.
!> Supports s, d, sc, dz.
!> See also: [[mfi_nrm2]], [[snrm2]],[[dnrm2]],[[scnrm2]],[[dznrm2]].
interface f77_nrm2
!> Original interface for SNRM2
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function snrm2(n, x, incx)
    import :: REAL32
    real(REAL32) :: snrm2
    real(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for DNRM2
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function dnrm2(n, x, incx)
    import :: REAL64
    real(REAL64) :: dnrm2
    real(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for SCNRM2
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function scnrm2(n, x, incx)
    import :: REAL32
    real(REAL32) :: scnrm2
    complex(REAL32), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
!> Original interface for DZNRM2
!> See also: [[mfi_nrm2]], [[f77_nrm2]].
pure function dznrm2(n, x, incx)
    import :: REAL64
    real(REAL64) :: dznrm2
    complex(REAL64), intent(in) :: x(*)
    integer, intent(in) :: n
    integer, intent(in) :: incx
end function
end interface
!> Generic old style interface for ROT.
!> Supports s, d, c, z, cs, zd.
!> See also: [[mfi_rot]], [[srot]],[[drot]],[[crot]],[[zrot]],[[csrot]],[[zdrot]].
interface f77_rot
!> Original interface for SROT
!> See also: [[mfi_rot]], [[f77_rot]].
!> SROT applies a plane rotation.
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
!> Original interface for DROT
!> See also: [[mfi_rot]], [[f77_rot]].
!> DROT applies a plane rotation.
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
!> Original interface for CROT
!> See also: [[mfi_rot]], [[f77_rot]].
!> CROT applies a plane rotation.
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
!> Original interface for ZROT
!> See also: [[mfi_rot]], [[f77_rot]].
!> ZROT applies a plane rotation.
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
!> Original interface for CSROT
!> See also: [[mfi_rot]], [[f77_rot]].
!> CSROT applies a plane rotation.
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
!> Original interface for ZDROT
!> See also: [[mfi_rot]], [[f77_rot]].
!> ZDROT applies a plane rotation.
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
!> Generic old style interface for ROTG.
!> Supports s, d, c, z.
!> See also: [[mfi_rotg]], [[srotg]],[[drotg]],[[crotg]],[[zrotg]].
interface f77_rotg
!> Original interface for SROTG
!> See also: [[mfi_rotg]], [[f77_rotg]].
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

!> Original interface for DROTG
!> See also: [[mfi_rotg]], [[f77_rotg]].
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

!> Original interface for CROTG
!> See also: [[mfi_rotg]], [[f77_rotg]].
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

!> Original interface for ZROTG
!> See also: [[mfi_rotg]], [[f77_rotg]].
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
!> Generic old style interface for ROTM.
!> Supports s, d.
!> See also: [[mfi_rotm]], [[srotm]],[[drotm]].
interface f77_rotm
!> Original interface for SROTM
!> See also: [[mfi_rotm]], [[f77_rotm]].
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
!> Original interface for DROTM
!> See also: [[mfi_rotm]], [[f77_rotm]].
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
!> Generic old style interface for ROTMG.
!> Supports s, d.
!> See also: [[mfi_rotmg]], [[srotmg]],[[drotmg]].
interface f77_rotmg
!> Original interface for SROTMG
!> See also: [[mfi_rotmg]], [[f77_rotmg]].
pure subroutine srotmg(d1, d2, x1, y1, param)
    import :: REAL32
    integer, parameter :: wp = REAL32
    real(REAL32), intent(in) :: y1
    real(REAL32), intent(out) :: param(5)
    real(REAL32), intent(inout) :: d1
    real(REAL32), intent(inout) :: d2
    real(REAL32), intent(inout) :: x1
end subroutine
!> Original interface for DROTMG
!> See also: [[mfi_rotmg]], [[f77_rotmg]].
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
!> Generic old style interface for SCAL.
!> Supports s, d, c, z, cs, zd.
!> See also: [[mfi_scal]], [[sscal]],[[dscal]],[[cscal]],[[zscal]],[[csscal]],[[zdscal]].
interface f77_scal
!> Original interface for SSCAL
!> See also: [[mfi_scal]], [[f77_scal]].
!> SSCAL scales a vector by a constant.
pure subroutine sscal(n, a, x, incx)
    import :: REAL32
    real(REAL32), intent(inout) :: x(*)
    real(REAL32), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for DSCAL
!> See also: [[mfi_scal]], [[f77_scal]].
!> DSCAL scales a vector by a constant.
pure subroutine dscal(n, a, x, incx)
    import :: REAL64
    real(REAL64), intent(inout) :: x(*)
    real(REAL64), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for CSCAL
!> See also: [[mfi_scal]], [[f77_scal]].
!> CSCAL scales a vector by a constant.
pure subroutine cscal(n, a, x, incx)
    import :: REAL32
    complex(REAL32), intent(inout) :: x(*)
    complex(REAL32), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for ZSCAL
!> See also: [[mfi_scal]], [[f77_scal]].
!> ZSCAL scales a vector by a constant.
pure subroutine zscal(n, a, x, incx)
    import :: REAL64
    complex(REAL64), intent(inout) :: x(*)
    complex(REAL64), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for CSSCAL
!> See also: [[mfi_scal]], [[f77_scal]].
!> CSSCAL scales a vector by a constant.
pure subroutine csscal(n, a, x, incx)
    import :: REAL32
    complex(REAL32), intent(inout) :: x(*)
    real(REAL32), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
!> Original interface for ZDSCAL
!> See also: [[mfi_scal]], [[f77_scal]].
!> ZDSCAL scales a vector by a constant.
pure subroutine zdscal(n, a, x, incx)
    import :: REAL64
    complex(REAL64), intent(inout) :: x(*)
    real(REAL64), intent(in) :: a
    integer, intent(in) :: n
    integer, intent(in) :: incx
end subroutine
end interface
!> Generic old style interface for GBMV.
!> Supports s, d, c, z.
!> See also: [[mfi_gbmv]], [[sgbmv]],[[dgbmv]],[[cgbmv]],[[zgbmv]].
interface f77_gbmv
!> Original interface for SGBMV
!> See also: [[mfi_gbmv]], [[f77_gbmv]].
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
!> Original interface for DGBMV
!> See also: [[mfi_gbmv]], [[f77_gbmv]].
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
!> Original interface for CGBMV
!> See also: [[mfi_gbmv]], [[f77_gbmv]].
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
!> Original interface for ZGBMV
!> See also: [[mfi_gbmv]], [[f77_gbmv]].
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
!> Generic old style interface for GEMV.
!> Supports s, d, c, z.
!> See also: [[mfi_gemv]], [[sgemv]],[[dgemv]],[[cgemv]],[[zgemv]].
interface f77_gemv
!> Original interface for SGEMV
!> See also: [[mfi_gemv]], [[f77_gemv]].
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
!> Original interface for DGEMV
!> See also: [[mfi_gemv]], [[f77_gemv]].
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
!> Original interface for CGEMV
!> See also: [[mfi_gemv]], [[f77_gemv]].
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
!> Original interface for ZGEMV
!> See also: [[mfi_gemv]], [[f77_gemv]].
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
!> Generic old style interface for GER.
!> Supports s, d.
!> See also: [[mfi_ger]], [[sger]],[[dger]].
interface f77_ger
!> Original interface for SGER
!> See also: [[mfi_ger]], [[f77_ger]].
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
!> Original interface for DGER
!> See also: [[mfi_ger]], [[f77_ger]].
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
!> Generic old style interface for GERC.
!> Supports c, z.
!> See also: [[mfi_gerc]], [[cgerc]],[[zgerc]].
interface f77_gerc
!> Original interface for CGERC
!> See also: [[mfi_gerc]], [[f77_gerc]].
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
!> Original interface for ZGERC
!> See also: [[mfi_gerc]], [[f77_gerc]].
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
!> Generic old style interface for GERU.
!> Supports c, z.
!> See also: [[mfi_geru]], [[cgeru]],[[zgeru]].
interface f77_geru
!> Original interface for CGERU
!> See also: [[mfi_geru]], [[f77_geru]].
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
!> Original interface for ZGERU
!> See also: [[mfi_geru]], [[f77_geru]].
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
!> Generic old style interface for HBMV.
!> Supports c, z.
!> See also: [[mfi_hbmv]], [[chbmv]],[[zhbmv]].
interface f77_hbmv
!> Original interface for CHBMV
!> See also: [[mfi_hbmv]], [[f77_hbmv]].
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
!> Original interface for ZHBMV
!> See also: [[mfi_hbmv]], [[f77_hbmv]].
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
!> Generic old style interface for HEMV.
!> Supports c, z.
!> See also: [[mfi_hemv]], [[chemv]],[[zhemv]].
interface f77_hemv
!> Original interface for CHEMV
!> See also: [[mfi_hemv]], [[f77_hemv]].
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
!> Original interface for ZHEMV
!> See also: [[mfi_hemv]], [[f77_hemv]].
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
!> Generic old style interface for HER.
!> Supports c, z.
!> See also: [[mfi_her]], [[cher]],[[zher]].
interface f77_her
!> Original interface for CHER
!> See also: [[mfi_her]], [[f77_her]].
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
!> Original interface for ZHER
!> See also: [[mfi_her]], [[f77_her]].
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
!> Generic old style interface for HER2.
!> Supports c, z.
!> See also: [[mfi_her2]], [[cher2]],[[zher2]].
interface f77_her2
!> Original interface for CHER2
!> See also: [[mfi_her2]], [[f77_her2]].
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
!> Original interface for ZHER2
!> See also: [[mfi_her2]], [[f77_her2]].
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
!> Generic old style interface for HPMV.
!> Supports c, z.
!> See also: [[mfi_hpmv]], [[chpmv]],[[zhpmv]].
interface f77_hpmv
!> Original interface for CHPMV
!> See also: [[mfi_hpmv]], [[f77_hpmv]].
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
!> Original interface for ZHPMV
!> See also: [[mfi_hpmv]], [[f77_hpmv]].
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
!> Generic old style interface for HPR.
!> Supports c, z.
!> See also: [[mfi_hpr]], [[chpr]],[[zhpr]].
interface f77_hpr
!> Original interface for CHPR
!> See also: [[mfi_hpr]], [[f77_hpr]].
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
!> Original interface for ZHPR
!> See also: [[mfi_hpr]], [[f77_hpr]].
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
!> Generic old style interface for HPR2.
!> Supports c, z.
!> See also: [[mfi_hpr2]], [[chpr2]],[[zhpr2]].
interface f77_hpr2
!> Original interface for CHPR2
!> See also: [[mfi_hpr2]], [[f77_hpr2]].
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
!> Original interface for ZHPR2
!> See also: [[mfi_hpr2]], [[f77_hpr2]].
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
!> Generic old style interface for SBMV.
!> Supports s, d.
!> See also: [[mfi_sbmv]], [[ssbmv]],[[dsbmv]].
interface f77_sbmv
!> Original interface for SSBMV
!> See also: [[mfi_sbmv]], [[f77_sbmv]].
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
!> Original interface for DSBMV
!> See also: [[mfi_sbmv]], [[f77_sbmv]].
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
!> Generic old style interface for SPMV.
!> Supports s, d.
!> See also: [[mfi_spmv]], [[sspmv]],[[dspmv]].
interface f77_spmv
!> Original interface for SSPMV
!> See also: [[mfi_spmv]], [[f77_spmv]].
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
!> Original interface for DSPMV
!> See also: [[mfi_spmv]], [[f77_spmv]].
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
!> Generic old style interface for SPR.
!> Supports s, d.
!> See also: [[mfi_spr]], [[sspr]],[[dspr]].
interface f77_spr
!> Original interface for SSPR
!> See also: [[mfi_spr]], [[f77_spr]].
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
!> Original interface for DSPR
!> See also: [[mfi_spr]], [[f77_spr]].
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
!> Generic old style interface for SPR2.
!> Supports s, d.
!> See also: [[mfi_spr2]], [[sspr2]],[[dspr2]].
interface f77_spr2
!> Original interface for SSPR2
!> See also: [[mfi_spr2]], [[f77_spr2]].
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
!> Original interface for DSPR2
!> See also: [[mfi_spr2]], [[f77_spr2]].
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
!> Generic old style interface for SYMV.
!> Supports s, d.
!> See also: [[mfi_symv]], [[ssymv]],[[dsymv]].
interface f77_symv
!> Original interface for SSYMV
!> See also: [[mfi_symv]], [[f77_symv]].
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
!> Original interface for DSYMV
!> See also: [[mfi_symv]], [[f77_symv]].
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
!> Generic old style interface for SYR.
!> Supports s, d.
!> See also: [[mfi_syr]], [[ssyr]],[[dsyr]].
interface f77_syr
!> Original interface for SSYR
!> See also: [[mfi_syr]], [[f77_syr]].
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
!> Original interface for DSYR
!> See also: [[mfi_syr]], [[f77_syr]].
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
!> Generic old style interface for SYR2.
!> Supports s, d.
!> See also: [[mfi_syr2]], [[ssyr2]],[[dsyr2]].
interface f77_syr2
!> Original interface for SSYR2
!> See also: [[mfi_syr2]], [[f77_syr2]].
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
!> Original interface for DSYR2
!> See also: [[mfi_syr2]], [[f77_syr2]].
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
!> Generic old style interface for TBMV.
!> Supports s, d, c, z.
!> See also: [[mfi_tbmv]], [[stbmv]],[[dtbmv]],[[ctbmv]],[[ztbmv]].
interface f77_tbmv
!> Original interface for STBMV
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
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
!> Original interface for DTBMV
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
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
!> Original interface for CTBMV
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
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
!> Original interface for ZTBMV
!> See also: [[mfi_tbmv]], [[f77_tbmv]].
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
!> Generic old style interface for TBSV.
!> Supports s, d, c, z.
!> See also: [[mfi_tbsv]], [[stbsv]],[[dtbsv]],[[ctbsv]],[[ztbsv]].
interface f77_tbsv
!> Original interface for STBSV
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
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
!> Original interface for DTBSV
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
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
!> Original interface for CTBSV
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
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
!> Original interface for ZTBSV
!> See also: [[mfi_tbsv]], [[f77_tbsv]].
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
!> Generic old style interface for TPMV.
!> Supports s, d, c, z.
!> See also: [[mfi_tpmv]], [[stpmv]],[[dtpmv]],[[ctpmv]],[[ztpmv]].
interface f77_tpmv
!> Original interface for STPMV
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
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
!> Original interface for DTPMV
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
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
!> Original interface for CTPMV
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
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
!> Original interface for ZTPMV
!> See also: [[mfi_tpmv]], [[f77_tpmv]].
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
!> Generic old style interface for TPSV.
!> Supports s, d, c, z.
!> See also: [[mfi_tpsv]], [[stpsv]],[[dtpsv]],[[ctpsv]],[[ztpsv]].
interface f77_tpsv
!> Original interface for STPSV
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
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
!> Original interface for DTPSV
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
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
!> Original interface for CTPSV
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
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
!> Original interface for ZTPSV
!> See also: [[mfi_tpsv]], [[f77_tpsv]].
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
!> Generic old style interface for TRMV.
!> Supports s, d, c, z.
!> See also: [[mfi_trmv]], [[strmv]],[[dtrmv]],[[ctrmv]],[[ztrmv]].
interface f77_trmv
!> Original interface for STRMV
!> See also: [[mfi_trmv]], [[f77_trmv]].
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
!> Original interface for DTRMV
!> See also: [[mfi_trmv]], [[f77_trmv]].
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
!> Original interface for CTRMV
!> See also: [[mfi_trmv]], [[f77_trmv]].
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
!> Original interface for ZTRMV
!> See also: [[mfi_trmv]], [[f77_trmv]].
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
!> Generic old style interface for TRSV.
!> Supports s, d, c, z.
!> See also: [[mfi_trsv]], [[strsv]],[[dtrsv]],[[ctrsv]],[[ztrsv]].
interface f77_trsv
!> Original interface for STRSV
!> See also: [[mfi_trsv]], [[f77_trsv]].
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
!> Original interface for DTRSV
!> See also: [[mfi_trsv]], [[f77_trsv]].
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
!> Original interface for CTRSV
!> See also: [[mfi_trsv]], [[f77_trsv]].
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
!> Original interface for ZTRSV
!> See also: [[mfi_trsv]], [[f77_trsv]].
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
!> Generic old style interface for GEMM.
!> Supports s, d, c, z.
!> See also: [[mfi_gemm]], [[sgemm]],[[dgemm]],[[cgemm]],[[zgemm]].
interface f77_gemm
!> Original interface for SGEMM
!> See also: [[mfi_gemm]], [[f77_gemm]].
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
!> Original interface for DGEMM
!> See also: [[mfi_gemm]], [[f77_gemm]].
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
!> Original interface for CGEMM
!> See also: [[mfi_gemm]], [[f77_gemm]].
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
!> Original interface for ZGEMM
!> See also: [[mfi_gemm]], [[f77_gemm]].
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
!> Generic old style interface for HEMM.
!> Supports c, z.
!> See also: [[mfi_hemm]], [[chemm]],[[zhemm]].
interface f77_hemm
!> Original interface for CHEMM
!> See also: [[mfi_hemm]], [[f77_hemm]].
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
!> Original interface for ZHEMM
!> See also: [[mfi_hemm]], [[f77_hemm]].
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
!> Generic old style interface for HERK.
!> Supports c, z.
!> See also: [[mfi_herk]], [[cherk]],[[zherk]].
interface f77_herk
!> Original interface for CHERK
!> See also: [[mfi_herk]], [[f77_herk]].
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
!> Original interface for ZHERK
!> See also: [[mfi_herk]], [[f77_herk]].
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
!> Generic old style interface for HER2K.
!> Supports c, z.
!> See also: [[mfi_her2k]], [[cher2k]],[[zher2k]].
interface f77_her2k
!> Original interface for CHER2K
!> See also: [[mfi_her2k]], [[f77_her2k]].
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
!> Original interface for ZHER2K
!> See also: [[mfi_her2k]], [[f77_her2k]].
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
!> Generic old style interface for SYMM.
!> Supports s, d.
!> See also: [[mfi_symm]], [[ssymm]],[[dsymm]].
interface f77_symm
!> Original interface for SSYMM
!> See also: [[mfi_symm]], [[f77_symm]].
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
!> Original interface for DSYMM
!> See also: [[mfi_symm]], [[f77_symm]].
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
!> Generic old style interface for SYRK.
!> Supports s, d.
!> See also: [[mfi_syrk]], [[ssyrk]],[[dsyrk]].
interface f77_syrk
!> Original interface for SSYRK
!> See also: [[mfi_syrk]], [[f77_syrk]].
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
!> Original interface for DSYRK
!> See also: [[mfi_syrk]], [[f77_syrk]].
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
!> Generic old style interface for SYR2K.
!> Supports s, d.
!> See also: [[mfi_syr2k]], [[ssyr2k]],[[dsyr2k]].
interface f77_syr2k
!> Original interface for SSYR2K
!> See also: [[mfi_syr2k]], [[f77_syr2k]].
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
!> Original interface for DSYR2K
!> See also: [[mfi_syr2k]], [[f77_syr2k]].
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
!> Generic old style interface for TRMM.
!> Supports s, d, c, z.
!> See also: [[mfi_trmm]], [[strmm]],[[dtrmm]],[[ctrmm]],[[ztrmm]].
interface f77_trmm
!> Original interface for STRMM
!> See also: [[mfi_trmm]], [[f77_trmm]].
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
!> Original interface for DTRMM
!> See also: [[mfi_trmm]], [[f77_trmm]].
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
!> Original interface for CTRMM
!> See also: [[mfi_trmm]], [[f77_trmm]].
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
!> Original interface for ZTRMM
!> See also: [[mfi_trmm]], [[f77_trmm]].
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
!> Generic old style interface for TRSM.
!> Supports s, d, c, z.
!> See also: [[mfi_trsm]], [[strsm]],[[dtrsm]],[[ctrsm]],[[ztrsm]].
interface f77_trsm
!> Original interface for STRSM
!> See also: [[mfi_trsm]], [[f77_trsm]].
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
!> Original interface for DTRSM
!> See also: [[mfi_trsm]], [[f77_trsm]].
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
!> Original interface for CTRSM
!> See also: [[mfi_trsm]], [[f77_trsm]].
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
!> Original interface for ZTRSM
!> See also: [[mfi_trsm]], [[f77_trsm]].
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

