#ifdef MFI_CUBLAS
#include <omp.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static cublasHandle_t *g_handles = NULL;
static int g_count = 0;
static int g_cublas_active = 0;
static int g_env_checked = 0;
static int g_initialized = 0;
static int g_debug = 0;

#define MFI_DPRINTF(fmt, ...) do { \
    if (g_debug) fprintf(stderr, "[MFI] " fmt, ##__VA_ARGS__); \
} while(0)

static void mfi_debug_init(void) {
    if (g_debug) return;
    const char *env = getenv("MFI_DEBUG");
    if (env && atoi(env) != 0) g_debug = 1;
}

void mfi_cublas_set_count(int count) {
    MFI_DPRINTF("cublas_set_count(%d)\n", count);
    if (g_handles) {
        for (int i = 0; i < g_count; i++) {
            if (g_handles[i]) cublasDestroy_v2(g_handles[i]);
        }
        free(g_handles);
        g_handles = NULL;
    }
    g_count = count;
    g_handles = count > 0 ? (cublasHandle_t *)calloc(count, sizeof(cublasHandle_t)) : NULL;
    MFI_DPRINTF("  allocated g_handles=%p count=%d\n", (void*)g_handles, g_count);
}

void mfi_cublas_init_handles(int *stat) {
    MFI_DPRINTF("cublas_init_handles()\n");
    *stat = 0;
    for (int i = 0; i < g_count; i++) {
        if (!g_handles[i]) {
            MFI_DPRINTF("  creating handle %d\n", i);
            *stat = (int)cublasCreate_v2(&g_handles[i]);
            MFI_DPRINTF("  cublasCreate_v2 returned %d, handle=%p\n", *stat, (void*)g_handles[i]);
            if (*stat != 0) {
                /* Try to initialize CUDA runtime first */
                cudaError_t cuda_err = cudaFree(0);
                MFI_DPRINTF("  cudaFree(0) returned %d\n", (int)cuda_err);
                if (cuda_err == cudaSuccess) {
                    *stat = (int)cublasCreate_v2(&g_handles[i]);
                    MFI_DPRINTF("  cublasCreate_v2 retry returned %d, handle=%p\n", *stat, (void*)g_handles[i]);
                }
            }
            if (*stat != 0) return;
            *stat = (int)cublasSetPointerMode_v2(g_handles[i], CUBLAS_POINTER_MODE_HOST);
            if (*stat != 0) return;
        }
    }
}

void mfi_cublas_get_thread_handle(void **out_handle, int *stat) {
    *stat = 0;
    int tid = omp_get_thread_num();
    MFI_DPRINTF("get_thread_handle() tid=%d count=%d g_handles=%p\n", tid, g_count, (void*)g_handles);
    if (tid >= 0 && tid < g_count) {
        MFI_DPRINTF("  returning handle[%d]=%p\n", tid, (void*)g_handles[tid]);
        *out_handle = (void *)g_handles[tid];
    } else {
        MFI_DPRINTF("  tid out of range, returning NULL\n");
        *out_handle = NULL;
        *stat = -1;
    }
}

void mfi_cublas_finalize_all(int *stat) {
    MFI_DPRINTF("cublas_finalize_all()\n");
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
    g_cublas_active = 0;
    g_initialized = 0;
}

int mfi_cublas_is_active(void) { return g_cublas_active; }
void mfi_cublas_set_active(int v) { g_cublas_active = v; }
int mfi_cublas_env_checked(void) { return g_env_checked; }
void mfi_cublas_set_env_checked(int v) { g_env_checked = v; }
int mfi_cublas_is_initialized(void) { return g_initialized; }
void mfi_cublas_set_initialized(int v) { g_initialized = v; }
int mfi_cublas_handle_count(void) { return g_count; }
void mfi_cublas_set_handle_count(int v) { g_count = v; }
void* mfi_cublas_get_handle(int i) { return (void *)(i >= 0 && i < g_count ? g_handles[i] : NULL); }

int mfi_cublas_read_env(void) {
    const char *env = getenv("MFI_USE_CUBLAS");
    int r = env ? atoi(env) : 0;
    MFI_DPRINTF("read_env() -> %d (MFI_USE_CUBLAS=%s)\n", r, env ? env : "(null)");
    return r;
}

int mfi_cublas_read_omp_threads(void) {
    const char *env = getenv("OMP_NUM_THREADS");
    int r = env ? atoi(env) : 0;
    MFI_DPRINTF("read_omp_threads() -> %d (OMP_NUM_THREADS=%s)\n", r, env ? env : "(null)");
    return r;
}

void mfi_cublas_lazy_init(void) {
    mfi_debug_init();
    MFI_DPRINTF("lazy_init() init=%d checked=%d active=%d count=%d\n",
                g_initialized, g_env_checked, g_cublas_active, g_count);
    if (g_initialized) return;
    if (!g_env_checked) {
        g_env_checked = 1;
        g_cublas_active = mfi_cublas_read_env();
    }
    if (g_cublas_active && g_count == 0) {
        int threads = mfi_cublas_read_omp_threads();
        if (threads < 1) threads = 1;
        mfi_cublas_set_count(threads);
        int stat;
        mfi_cublas_init_handles(&stat);
        if (stat != 0) {
            MFI_DPRINTF("  cublas init failed with stat=%d, falling back to CPU\n", stat);
            g_cublas_active = 0;
            /* Don't set g_initialized so it can be retried */
            return;
        }
    }
    g_initialized = 1;
    MFI_DPRINTF("lazy_init() done: active=%d count=%d\n", g_cublas_active, g_count);
}

void mfi_cublas_force_gpu(void) {
    mfi_debug_init();
    MFI_DPRINTF("force_gpu()\n");
    g_cublas_active = 1;
    g_env_checked = 1;
    g_initialized = 0;
    mfi_cublas_lazy_init();
}

void mfi_cublas_force_cpu(void) {
    mfi_debug_init();
    MFI_DPRINTF("force_cpu()\n");
    /* Destroy all existing handles to prevent resource leaks */
    if (g_handles) {
        for (int i = 0; i < g_count; i++) {
            if (g_handles[i]) {
                cublasDestroy_v2(g_handles[i]);
                g_handles[i] = NULL;
            }
        }
        free(g_handles);
        g_handles = NULL;
    }
    g_cublas_active = 0;
    g_env_checked = 1;
    g_initialized = 0;
    g_count = 0;
}

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
