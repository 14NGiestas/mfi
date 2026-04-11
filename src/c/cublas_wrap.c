#ifdef MFI_CUBLAS
#include <omp.h>
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

/* Thread-safe handle lookup — pure because Fortran trusts bind(c) purity.
   The global array mfi_cublas_handles is managed by Fortran module code. */
void mfi_cublas_get_thread_handle(void **handles, int count, void **out_handle) {
    int tid = omp_get_thread_num();
    if (tid >= 0 && tid < count) {
        *out_handle = handles[tid];
    } else {
        *out_handle = NULL;
    }
}
#endif
