#:mute
#:include "common.fpp"

#:def geqrf_gerqf(NAME,TYPE,KIND)
pure subroutine ${NAME}$(m,n,a,lda,tau,work,lwork,info)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, a(lda,*))
@:args(${TYPE}$,   out, tau(*))
@:args(integer,    out, info)
@:args(integer,     in, m, n, lda, lwork)
@:args(${TYPE}$, inout, work(*))
end subroutine
#:enddef

#:def getrf(NAME,TYPE,KIND)
pure subroutine ${NAME}$(m,n,a,lda,ipiv,info)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, a(lda,*))
@:args(integer,    out, ipiv(*))
@:args(integer,    out, info)
@:args(integer,     in, m, n, lda)
end subroutine
#:enddef

#:def getri(NAME,TYPE,KIND)
pure subroutine ${NAME}$(n,a,lda,ipiv,work,lwork,info)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, a(lda,*))
@:args(${TYPE}$, inout, work(*))
@:args(integer,     in, ipiv(*))
@:args(integer,    out, info)
@:args(integer,     in, n, lda, lwork)
end subroutine
#:enddef

#:def getrs(NAME,TYPE,KIND)
pure subroutine ${NAME}$(trans,n,nrhs,a,lda,ipiv,b,ldb,info)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, a(lda,*))
@:args(${TYPE}$, inout, b(ldb,*))
@:args(character,   in, trans)
@:args(integer,     in, ipiv(*))
@:args(integer,    out, info)
@:args(integer,     in, n, nrhs, lda, ldb)
end subroutine
#:enddef

#:def hetrf(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, a, lda, ipiv, work, lwork, info)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$, inout, a(lda,*))
@:args(character,   in, uplo)
@:args(integer,     in, ipiv(*))
@:args(${TYPE}$, inout, work(*))
@:args(integer,    out, info)
@:args(integer,     in, n, lda, lwork)
end subroutine
#:enddef

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

#:def hegv(NAME,TYPE,KIND)
pure subroutine ${NAME}$(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  inout, a(lda,*))
@:args(${TYPE}$,  inout, b(ldb,*))
@:args(${REAL_TYPE}$, out, w(*))
@:args(integer,   out,   info)
@:args(character, in,    jobz, uplo)
@:args(integer,   in,    n, itype, lda, ldb, lwork)
@:args(${TYPE}$,  inout, work(*))
@:args(${REAL_TYPE}$, in, rwork(*))
end subroutine
#:enddef

#:def heevd(NAME,TYPE,KIND)
pure subroutine ${NAME}$(jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,      inout, a(lda,*))
@:args(${REAL_TYPE}$, out,   w(*))
@:args(integer,       out,   info)
@:args(character,     in,    jobz, uplo)
@:args(integer,       in,    n, lda, lwork, lrwork, liwork)
@:args(${TYPE}$,      inout, work(*))
@:args(${REAL_TYPE}$, inout, rwork(*))
@:args(integer,       inout, iwork(*))
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

#:def potrs(NAME,TYPE,KIND)
pure subroutine ${NAME}$(uplo, n, nrhs, a, lda, b, ldb, info)
    import :: ${KIND}$
@:parameter(integer, wp=${KIND}$)
@:args(${TYPE}$,  in, a(lda,*))
@:args(${TYPE}$,  in, b(ldb,*))
@:args(character, in,  uplo)
@:args(integer,   in,  n, nrhs, lda, ldb)
@:args(integer,   out, info)
end subroutine
#:enddef

#:endmute
module f77_lapack
use iso_fortran_env
implicit none

$:f77_interface('?geqrf',  DEFAULT_TYPES, geqrf_gerqf)
$:f77_interface('?gerqf',  DEFAULT_TYPES, geqrf_gerqf)
$:f77_interface('?getrf',  DEFAULT_TYPES, getrf)
$:f77_interface('?getri',  DEFAULT_TYPES, getri)
$:f77_interface('?getrs',  DEFAULT_TYPES, getrs)
$:f77_interface('?hetrf',  COMPLEX_TYPES, hetrf)
$:f77_interface('?hegv',   COMPLEX_TYPES, hegv)
$:f77_interface('?heevd',  COMPLEX_TYPES, heevd)
$:f77_interface('?gesvd',  DEFAULT_TYPES, gesvd)
$:f77_interface('?potrf',  DEFAULT_TYPES, potrf_potri)
$:f77_interface('?potri',  DEFAULT_TYPES, potrf_potri)
$:f77_interface('?potrs',  DEFAULT_TYPES, potrs)

    interface f77_xerbla
        pure subroutine xerbla(name,info)
            character(*), intent(in) :: name
            integer, intent(in) :: info
        end subroutine
    end interface f77_xerbla

end module

