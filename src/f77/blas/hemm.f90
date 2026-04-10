module f77_blas_hemm
    use iso_fortran_env
    use iso_c_binding
    implicit none

!> Generic old style interface for HEMM.
!> Supports c, z.
!> See also: [[mfi_hemm]], [[chemm]], [[zhemm]].
interface f77_hemm
!> Original interface for CHEMM
!> See also: [[mfi_hemm]], [[hemm]].
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
!> See also: [[mfi_hemm]], [[hemm]].
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
end module

