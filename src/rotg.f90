module f77_blas_rotg
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for ROTG.
!> Supports s, d, c, z.
!> See also: [[mfi_rotg]], [[srotg]], [[drotg]], [[crotg]], [[zrotg]].
interface f77_rotg
!> Original interface for SROTG
!> See also: [[mfi_rotg]], [[rotg]].
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
!> See also: [[mfi_rotg]], [[rotg]].
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
!> See also: [[mfi_rotg]], [[rotg]].
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
!> See also: [[mfi_rotg]], [[rotg]].
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
end module

