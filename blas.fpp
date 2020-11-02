#:include "common.fpp"
module f77_blas
#:for name, types in BLAS_ROUTINES.items()
$:f77_interface(name, types)
#:endfor
end module

module mfi_blas
use f77_blas
use iso_fortran_env
implicit none

#:for name, types in BLAS_ROUTINES.items()
$:mfi_interface(name, types)
#:endfor

contains

#:block mfi_implement('?gemv')
#:contains args
(a,x,y,alpha,beta,trans)
#:contains code
    {T}, intent(in) :: a(:,:)
    {T}, intent(in) :: x(:)
    {T}, intent(inout) :: y(:)
    character(len=1), intent(in), optional :: trans
    {T}, intent(in), optional :: alpha
    {T}, intent(in), optional :: beta
    character(len=1) :: local_trans
    {T} :: local_alpha
    {T} :: local_beta
    integer :: incx
    integer :: incy
    integer :: m
    integer :: n
    integer :: lda
    intrinsic max, present, size
    if(present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1
    endif
    if(present(beta)) then
        local_beta = beta
    else
        local_beta = 0
    endif
    if(present(trans)) then
        local_trans = trans
    else
        local_trans = 'n'
    endif
    incx = 1
    incy = 1
    lda = max(1,size(a,1))
    m = size(a,1)
    n = size(a,2)
    ! <<< call blas77 routine >>>
    call f77_gemv(local_trans,m,n,local_alpha,a,lda,x,incx,local_beta,y,incy)
#:endblock

#:block mfi_implement('?gemm')
#:contains args
(a,b,c,transa,transb,alpha,beta)
#:contains code
    {T}, intent(in) :: a(:,:)
    {T}, intent(in) :: b(:,:)
    {T}, intent(inout) :: c(:,:)
    character(len=1), intent(in), optional :: transa
    character(len=1), intent(in), optional :: transb
    {T}, intent(in), optional :: alpha
    {T}, intent(in), optional :: beta
    character(len=1) :: local_transa
    character(len=1) :: local_transb
    {T} :: local_alpha
    {T} :: local_beta
    integer :: m
    integer :: n
    integer :: k
    integer :: lda
    integer :: ldb
    integer :: ldc
    intrinsic max, present, size
    if(present(alpha)) then
        local_alpha = alpha
    else
        local_alpha = 1
    endif
    if(present(beta)) then
        local_beta = beta
    else
        local_beta = 0
    endif
    if(present(transa)) then
        local_transa = transa
    else
        local_transa = 'n'
    endif
    if(present(transb)) then
        local_transb = transb
    else
        local_transb = 'n'
    endif
    if((local_transa.eq.'n'.or.local_transa.eq.'n')) then
        k = size(a,2)
    else
        k = size(a,1)
    endif
    lda = max(1,size(a,1))
    ldb = max(1,size(b,1))
    ldc = max(1,size(c,1))
    m = size(c,1)
    n = size(c,2)
    call f77_gemm(local_transa,local_transb,m,n,k,local_alpha,a,lda,b,ldb,local_beta,c,ldc)
#:endblock
end module
