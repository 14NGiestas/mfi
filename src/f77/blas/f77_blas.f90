module f77_blas
    use iso_fortran_env
    use iso_c_binding
    use f77_blas_copy
    use f77_blas_swap
    use f77_blas_axpy
    use f77_blas_dot
    use f77_blas_dotc
    use f77_blas_dotu
    use f77_blas_asum
    use f77_blas_nrm2
    use f77_blas_rot
    use f77_blas_rotg
    use f77_blas_rotm
    use f77_blas_rotmg
    use f77_blas_scal
    use f77_blas_gbmv
    use f77_blas_gemv
    use f77_blas_ger
    use f77_blas_gerc
    use f77_blas_geru
    use f77_blas_hbmv
    use f77_blas_hemv
    use f77_blas_her
    use f77_blas_her2
    use f77_blas_hpmv
    use f77_blas_hpr
    use f77_blas_hpr2
    use f77_blas_sbmv
    use f77_blas_spmv
    use f77_blas_spr
    use f77_blas_spr2
    use f77_blas_symv
    use f77_blas_syr
    use f77_blas_syr2
    use f77_blas_tbmv
    use f77_blas_tbsv
    use f77_blas_tpmv
    use f77_blas_tpsv
    use f77_blas_trmv
    use f77_blas_trsv
    use f77_blas_gemm
    use f77_blas_hemm
    use f77_blas_herk
    use f77_blas_her2k
    use f77_blas_symm
    use f77_blas_syrk
    use f77_blas_syr2k
    use f77_blas_trmm
    use f77_blas_trsm
    use f77_blas_iamax
    use f77_blas_iamin
    use f77_blas_cublas
    implicit none

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
end module
