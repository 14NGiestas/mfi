#:include "common.fpp"

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

$:f77_interface('?potrf',  DEFAULT_TYPES, potrf_potri)
$:f77_interface('?potri',  DEFAULT_TYPES, potrf_potri)

end module

