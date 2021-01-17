#:include "common.fpp"

#:def gesvd(NAME,TYPE,KIND)
#:if TYPE == COMPLEX_TYPE
pure subroutine ${NAME}$(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
#:else
pure subroutine ${NAME}$(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,info)
#:endif
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  inout, a(lda,*))
@:args(${REAL_TYPE}$, out, s(*))
@:args(${TYPE}$,  out,   u(ldu,*), vt(ldvt,*))
@:args(integer,   out,   info)
@:args(character, in,    jobu, jobvt)
@:args(integer,   in,    m, n, lda, ldu, ldvt, lwork)
@:args(${TYPE}$,  inout, work(*))
#:if TYPE == COMPLEX_TYPE
@:args(${REAL_TYPE}$, in, rwork(*))
#:endif
end subroutine
#:enddef

#:def potrf_potri(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, a, lda, info)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in, a(lda,*))
@:args(character, in,  uplo)
@:args(integer,   in,  n, lda)
@:args(integer,   out, info)
end subroutine
#:enddef

module f77_lapack
use iso_fortran_env
implicit none

$:f77_interface('?gesvd',  DEFAULT_TYPES, gesvd)
$:f77_interface('?potrf',  DEFAULT_TYPES, potrf_potri)
$:f77_interface('?potri',  DEFAULT_TYPES, potrf_potri)

end module

