/* hip_wrap.c: native HIP and rocBLAS backend for MFI.
 *
 * This file mirrors src/c/cublas_wrap.c with rocBLAS calls. This path uses no
 * CUDA and no ZLUDA. The linker needs librocblas and libamdhip64.
 *
 * The HIP headers need one platform macro at compile time. Set it here so a
 * plain C compiler builds this file without extra flags. */
#ifndef __HIP_PLATFORM_AMD__
#define __HIP_PLATFORM_AMD__
#endif

#ifdef MFI_HIP
#include <omp.h>
#include <hip/hip_runtime.h>
#include <rocblas/rocblas.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* State:
 *
 *  Handles stay alive across CPU and GPU switches.
 *  force_cpu() only clears the active flag. Handles persist.
 *  force_gpu() only sets the active flag. The code reuses old handles.
 *  finalize() is the only function that destroys handles.
 *
 *  g_threads  = desired handle count (set by hand or from the env)
 *  g_active   = 1 means GPU path, 0 means CPU fallback
 *  g_ready    = 1 means handles exist (active or not)
 */

static rocblas_handle *g_handles  = NULL;
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

/* Helpers */

static int mfi_read_env_int(const char *name) {
    const char *env = getenv(name);
    return env ? atoi(env) : 0;
}

static int mfi_hip_ok(rocblas_status stat) {
    return (int)stat == (int)rocblas_status_success;
}

static void mfi_ensure_handles(void) {
    /* Called internally. This creates handles when they do not exist.
     * It uses OMP_NUM_THREADS (default 1) when g_threads is unset. */
    if (g_ready) return;

    if (g_threads < 1) {
        g_threads = mfi_read_env_int("OMP_NUM_THREADS");
        if (g_threads < 1) g_threads = 1;
    }

    g_handles = (rocblas_handle *)calloc(g_threads, sizeof(rocblas_handle));
    if (!g_handles) return;

    for (int i = 0; i < g_threads; i++) {
        rocblas_status stat = rocblas_create_handle(&g_handles[i]);
        if (!mfi_hip_ok(stat)) {
            /* Try starting the HIP runtime, then retry once */
            if ((int)hipInit(0) == (int)hipSuccess)
                stat = rocblas_create_handle(&g_handles[i]);
        }
        if (!mfi_hip_ok(stat)) {
            /* Partial init failure. Clean up and give up. */
            for (int j = 0; j < i; j++)
                if (g_handles[j]) rocblas_destroy_handle(g_handles[j]);
            free(g_handles);
            g_handles = NULL;
            g_threads = 0;
            return;
        }
        /* Scalars alpha and beta live on the host, as in the cuBLAS path. */
        rocblas_set_pointer_mode(g_handles[i], rocblas_pointer_mode_host);
    }
    g_ready = 1;
    MFI_DPRINTF("hip handles created: threads=%d\n", g_threads);
}

/* Public API */

/** Set the number of rocBLAS handles. Call this before the first GPU use.
 *  The code destroys and recreates old handles when the count changes. */
void mfi_hip_set_threads(int count) {
    mfi_debug_init();
    if (count < 1) count = 1;
    MFI_DPRINTF("hip set_threads(%d) [was %d, ready=%d]\n", count, g_threads, g_ready);

    /* Destroy existing handles when the count changed */
    if (g_ready && count != g_threads) {
        for (int i = 0; i < g_threads; i++)
            if (g_handles[i]) rocblas_destroy_handle(g_handles[i]);
        free(g_handles);
        g_handles = NULL;
        g_ready = 0;
    }

    g_threads = count;

    /* Create handles now when the GPU path is active */
    if (g_active) mfi_ensure_handles();
}

/** Switch to GPU mode. This reuses old handles or creates them on first use. */
void mfi_hip_force_gpu(void) {
    mfi_debug_init();
    MFI_DPRINTF("hip force_gpu()\n");
    g_active = 1;

    /* Find out how many handles to create */
    if (g_threads < 1) {
        g_threads = mfi_read_env_int("OMP_NUM_THREADS");
        if (g_threads < 1) g_threads = 1;
    }

    /* Create handles when needed */
    mfi_ensure_handles();
}

/** Switch to CPU mode. The code keeps the handles warm for a fast re-switch. */
void mfi_hip_force_cpu(void) {
    mfi_debug_init();
    MFI_DPRINTF("hip force_cpu()\n");
    g_active = 0;
    /* Handles stay alive. force_gpu reuses them at once. */
}

/** Destroy all handles and reset the state. */
void mfi_hip_finalize_all(int *stat) {
    mfi_debug_init();
    MFI_DPRINTF("hip finalize_all()\n");
    *stat = 0;
    if (g_handles) {
        for (int i = 0; i < g_threads; i++)
            if (g_handles[i]) {
                rocblas_status s = rocblas_destroy_handle(g_handles[i]);
                if (!mfi_hip_ok(s) && *stat == 0) *stat = (int)s;
                g_handles[i] = NULL;
            }
        free(g_handles);
        g_handles = NULL;
    }
    g_threads = 0;
    g_active  = 0;
    g_ready   = 0;
}

/* Query functions */

int  mfi_hip_is_active(void)          { return g_active; }
int  mfi_hip_handle_count(void)       { return g_threads; }
int  mfi_hip_is_ready(void)           { return g_ready; }
void *mfi_hip_get_handle(int i)       {
    return (i >= 0 && i < g_threads && g_handles) ? (void *)g_handles[i] : NULL;
}

/** Get the handle for the current OpenMP thread. This falls back to handle[0]
 *  for a non-OpenMP caller or a mismatched thread count. */
void mfi_hip_get_thread_handle(void **out_handle, int *stat) {
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
        /* Fallback: use handle 0. This is safe for a single-threaded caller. */
        *out_handle = (void *)g_handles[0];
    }
}

/* Legacy and compatibility shims */

int mfi_hip_read_env(void)           { return mfi_read_env_int("MFI_USE_HIP"); }
int mfi_hip_read_omp_threads(void)   { return mfi_read_env_int("OMP_NUM_THREADS"); }
int mfi_hip_is_initialized(void)     { return g_ready; }
int mfi_hip_env_checked(void)        { return 1; }  /* always "checked" here */
void mfi_hip_set_active(int v)       { g_active = v; }
void mfi_hip_set_env_checked(int v)  { (void)v; }    /* no-op here */
void mfi_hip_set_initialized(int v)  { (void)v; }    /* no-op here */
void mfi_hip_set_handle_count(int v) { mfi_hip_set_threads(v); }

/** Legacy lazy init. This stays for callers that expect it.
 *  It activates the GPU path when the env var says so and handles are absent. */
void mfi_hip_lazy_init(void) {
    mfi_debug_init();
    if (g_ready) return;
    if (!g_active) {
        g_active = mfi_read_env_int("MFI_USE_HIP");
    }
    if (g_active) mfi_ensure_handles();
}

/* Raw HIP and rocBLAS helpers */

/** Alias of mfi_hip_get_thread_handle. The short name stays for callers. */
void mfi_hip_handle_get(void **out_handle, int *stat) {
    mfi_hip_get_thread_handle(out_handle, stat);
}

/** Alias of mfi_hip_finalize_all. The short name stays for callers. */
void mfi_hip_finalize(int *stat) {
    mfi_hip_finalize_all(stat);
}

void mfi_hip_malloc(void **devPtr, size_t size, int *stat) {
    *stat = (int)hipMalloc(devPtr, size);
}

void mfi_hip_free(void *devPtr, int *stat) {
    *stat = (int)hipFree(devPtr);
}

void mfi_hip_create(void **handle, int *stat) {
    *stat = (int)rocblas_create_handle((rocblas_handle *)handle);
}

void mfi_hip_destroy(void *handle, int *stat) {
    *stat = (int)rocblas_destroy_handle((rocblas_handle)handle);
}

void mfi_hip_set_pointer_mode(void *handle, int mode, int *stat) {
    *stat = (int)rocblas_set_pointer_mode((rocblas_handle)handle,
                                          (rocblas_pointer_mode)mode);
}

void mfi_hip_memcpy(void *dst, const void *src, size_t count, int kind, int *stat) {
    *stat = (int)hipMemcpy(dst, src, count, (hipMemcpyKind)kind);
}
#endif
