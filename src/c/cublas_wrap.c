#ifdef MFI_CUBLAS
#include <omp.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <stdlib.h>

static cublasHandle_t *g_handles = NULL;
static int g_count = 0;

void mfi_cublas_set_count(int count) {
    if (g_handles) {
        for (int i = 0; i < g_count; i++) {
            if (g_handles[i]) cublasDestroy_v2(g_handles[i]);
        }
        free(g_handles);
    }
    g_count = count;
    g_handles = count > 0 ? (cublasHandle_t *)calloc(count, sizeof(cublasHandle_t)) : NULL;
}

void mfi_cublas_init_handles(int *stat) {
    *stat = 0;
    for (int i = 0; i < g_count; i++) {
        if (!g_handles[i]) {
            *stat = (int)cublasCreate_v2(&g_handles[i]);
            if (*stat != 0) return;
            *stat = (int)cublasSetPointerMode_v2(g_handles[i], CUBLAS_POINTER_MODE_HOST);
            if (*stat != 0) return;
        }
    }
}

void mfi_cublas_get_thread_handle(void **out_handle, int *stat) {
    *stat = 0;
    int tid = omp_get_thread_num();
    if (tid >= 0 && tid < g_count) {
        *out_handle = (void *)g_handles[tid];
    } else {
        *out_handle = NULL;
        *stat = -1;
    }
}

void mfi_cublas_finalize_all(int *stat) {
    *stat = 0;
    if (g_handles) {
        for (int i = 0; i < g_count; i++) {
            if (g_handles[i]) {
                *stat = (int)cublasDestroy_v2(g_handles[i]);
                g_handles[i] = NULL;
            }
        }
        free(g_handles);
        g_handles = NULL;
        g_count = 0;
    }
}

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

void mfi_cuda_memcpy(void *dst, const void *src, size_t count, int kind, int *stat) {
    *stat = (int)cudaMemcpy(dst, src, count, (enum cudaMemcpyKind)kind);
}
#endif
