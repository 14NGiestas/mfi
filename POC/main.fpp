#:include "../common.fpp"
program main
    use mfi_blas
    implicit none
    integer, parameter :: wp = REAL64
    integer, parameter :: n = 2048*2
    real(wp), allocatable :: A(:,:), B(:,:), C(:,:), ref(:,:)
    allocate(A(n,n))
    allocate(B(n,n))
    allocate(C(n,n))
    allocate(ref(n,n))
    c = 0.0; ref = 0.0

    call random_seed()
    call random_number(A)
    call random_number(B)

    call mfi_init
    @:timeit("MFI:    ", { call mfi_gemm(A,B,C,transa='T') })
    @:timeit("MATMUL: ", { ref = matmul(transpose(A),B)    })

    print*, 'L2', norm2(abs(ref - C)), 'is almost equal? ', all(is_almost_equal(ref,C))
contains
  logical pure elemental function is_almost_equal(x, y)
      real(wp), intent(in) :: x, y
      is_almost_equal = abs(x-y) < 10**6*epsilon(x)
  end function
end program
