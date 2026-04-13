#ifdef MFI_CUBLAS
#include <omp.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ── State ───────────────────────────────────────────────────────
 *
 *  Handles are kept alive across CPU ↔ GPU switches.
 *  force_cpu() only clears the active flag; handles persist.
 *  force_gpu() just sets the active flag — handles are reused.
 *  finalize() is the only function that actually destroys handles.
 *
 *  g_threads  = desired handle count (set explicitly or from env)
 *  g_active   = 1 → GPU path, 0 → CPU fallback
 *  g_ready    = 1 → handles have been created (may or may not be active)
 * ────────────────────────────────────────────────────────────── */

static cublasHandle_t *g_handles  = NULL;
static int             g_threads  = 0;   /* desired / actual handle count    */
static int             g_active   = 0;   /* is GPU path currently active?    */
static int             g_ready    = 0;   /* have handles been created?       */
static int             g_debug    = 0;

#define MFI_DPRINTF(fmt, ...) do { \
    if (g_debug) fprintf(stderr, "[MFI] " fmt, ##__VA_ARGS__); \
} while(0)

static void mfi_debug_init(void) {
    if (g_debug) return;
    const char *env = getenv("MFI_DEBUG");
    if (env && atoi(env) != 0) g_debug = 1;
}

/* ── Helpers ───────────────────────────────────────────────────── */

static int mfi_read_env_int(const char *name) {
    const char *env = getenv(name);
    return env ? atoi(env) : 0;
}

static void mfi_ensure_handles(void) {
    /* Called internally — creates handles if they don't exist yet.
     * Uses OMP_NUM_THREADS env var (default 1) if not set explicitly. */
    if (g_ready) return;

    if (g_threads < 1) {
        g_threads = mfi_read_env_int("OMP_NUM_THREADS");
        if (g_threads < 1) g_threads = 1;
    }

    g_handles = (cublasHandle_t *)calloc(g_threads, sizeof(cublasHandle_t));
    if (!g_handles) return;

    for (int i = 0; i < g_threads; i++) {
        int stat = (int)cublasCreate_v2(&g_handles[i]);
        if (stat != 0) {
            /* Try waking up the CUDA runtime */
            if (cudaFree(0) == cudaSuccess)
                stat = (int)cublasCreate_v2(&g_handles[i]);
        }
        if (stat != 0) {
            /* Partial init failure — clean up and give up */
            for (int j = 0; j < i; j++)
                if (g_handles[j]) cublasDestroy_v2(g_handles[j]);
            free(g_handles);
            g_handles = NULL;
            g_threads = 0;
            return;
        }
        cublasSetPointerMode_v2(g_handles[i], CUBLAS_POINTER_MODE_HOST);
    }
    g_ready = 1;
    MFI_DPRINTF("handles created: threads=%d\n", g_threads);
}

/* ── Public API ────────────────────────────────────────────────── */

/** Set the number of cuBLAS handles (call before first GPU use, or after finalize).
 *  If called when handles already exist, old handles are destroyed and recreated. */
void mfi_cublas_set_threads(int count) {
    mfi_debug_init();
    if (count < 1) count = 1;
    MFI_DPRINTF("set_threads(%d) [was %d, ready=%d]\n", count, g_threads, g_ready);

    /* Destroy existing handles if count changed */
    if (g_ready && count != g_threads) {
        for (int i = 0; i < g_threads; i++)
            if (g_handles[i]) cublasDestroy_v2(g_handles[i]);
        free(g_handles);
        g_handles = NULL;
        g_ready = 0;
    }

    g_threads = count;

    /* If GPU is currently active, create handles immediately */
    if (g_active) mfi_ensure_handles();
}

/** Switch to GPU mode. Reuses existing handles or creates them on first use. */
void mfi_cublas_force_gpu(void) {
    mfi_debug_init();
    MFI_DPRINTF("force_gpu()\n");
    g_active = 1;

    /* Ensure we know how many handles to create */
    if (g_threads < 1) {
        g_threads = mfi_read_env_int("OMP_NUM_THREADS");
        if (g_threads < 1) g_threads = 1;
    }

    /* Create handles if needed */
    mfi_ensure_handles();
}

/** Switch to CPU mode.  Handles are preserved (warm) for fast re-switch. */
void mfi_cublas_force_cpu(void) {
    mfi_debug_init();
    MFI_DPRINTF("force_cpu()\n");
    g_active = 0;
    /* Handles stay alive — force_gpu will reuse them instantly */
}

/** Destroy all handles and reset all state. */
void mfi_cublas_finalize_all(int *stat) {
    mfi_debug_init();
    MFI_DPRINTF("finalize_all()\n");
    *stat = 0;
    if (g_handles) {
        for (int i = 0; i < g_threads; i++)
            if (g_handles[i]) {
                int s = (int)cublasDestroy_v2(g_handles[i]);
                if (s != 0 && *stat == 0) *stat = s;
                g_handles[i] = NULL;
            }
        free(g_handles);
        g_handles = NULL;
    }
    g_threads = 0;
    g_active  = 0;
    g_ready   = 0;
}

/* ── Query functions ───────────────────────────────────────────── */

int  mfi_cublas_is_active(void)       { return g_active; }
int  mfi_cublas_handle_count(void)    { return g_threads; }
int  mfi_cublas_is_ready(void)        { return g_ready; }
void *mfi_cublas_get_handle(int i)    {
    return (i >= 0 && i < g_threads && g_handles) ? (void *)g_handles[i] : NULL;
}

/** Get handle for current OpenMP thread. Falls back to handle[0] for
 *  non-OpenMP or mismatched thread counts. */
void mfi_cublas_get_thread_handle(void **out_handle, int *stat) {
    *stat = 0;
    if (!g_ready || !g_handles || g_threads < 1) {
        *out_handle = NULL;
        *stat = -1;
        return;
    }

    int tid = omp_get_thread_num();
    if (tid >= 0 && tid < g_threads) {
        *out_handle = (void *)g_handles[tid];
    } else {
        /* Fallback: use handle 0 (safe for single-threaded callers) */
        *out_handle = (void *)g_handles[0];
    }
}

/* ── Legacy / compatibility shims ──────────────────────────────── */

int mfi_cublas_read_env(void)           { return mfi_read_env_int("MFI_USE_CUBLAS"); }
int mfi_cublas_read_omp_threads(void)   { return mfi_read_env_int("OMP_NUM_THREADS"); }
int mfi_cublas_is_initialized(void)     { return g_ready; }
int mfi_cublas_env_checked(void)        { return 1; }  /* always "checked" in new design */
void mfi_cublas_set_active(int v)       { g_active = v; }
void mfi_cublas_set_env_checked(int v)  { (void)v; }    /* no-op in new design */
void mfi_cublas_set_initialized(int v)  { (void)v; }    /* no-op in new design */
void mfi_cublas_set_handle_count(int v) { mfi_cublas_set_threads(v); }

/** Legacy lazy init — kept for callers that expect it.
 *  Now just activates GPU if env var says so and handles aren't ready. */
void mfi_cublas_lazy_init(void) {
    mfi_debug_init();
    if (g_ready) return;
    if (!g_active) {
        g_active = mfi_read_env_int("MFI_USE_CUBLAS");
    }
    if (g_active) mfi_ensure_handles();
}

/* ── Raw CUDA/cuBLAS helpers ───────────────────────────────────── */

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
