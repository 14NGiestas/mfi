#ifdef MFI_CUBLAS_SUPPORT
#include <cuda_runtime.h>
#include <cublas_v2.h>

/* Wrapper: convert C function return value to output parameter.
   Allows Fortran to call as "pure subroutine" with correct ABI. */

void mfi_cuda_malloc(void **devPtr, size_t size, int *stat) {
    *stat = (int)cudaMalloc(devPtr, size);
}

void mfi_cuda_free(void *devPtr, int *stat) {
    *stat = (int)cudaFree(devPtr);
}

void mfi_cublas_create(void **handle, int *stat) {
    *stat = (int)cublasCreate_v2((cublasHandle_t *)handle);
}

void mfi_cublas_destroy(void *handle, int *stat) {
    *stat = (int)cublasDestroy_v2((cublasHandle_t)handle);
}

void mfi_cublas_set_pointer_mode(void *handle, int mode, int *stat) {
    *stat = (int)cublasSetPointerMode_v2((cublasHandle_t)handle,
                                          (cublasPointerMode_t)mode);
}
#endif
