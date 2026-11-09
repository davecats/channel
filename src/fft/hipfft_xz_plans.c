/* hipFFT plans carrying the JIT callbacks of hipfft_xz_callbacks.hip.
 *
 * The SPIR-V of the callbacks is built by CMake and embedded below with
 * .incbin (CHANNEL_XZ_CB_SPV is its path, as a string literal).  Each plan is created with
 * hipfftCreate, gets its callbacks with hipfftXtSetJITCallback and is only
 * then made (hipfftMakePlanMany): JIT callbacks must be set before the plan
 * is initialised.  The plans are returned to Fortran, which executes them
 * with the ordinary hipfftExec* calls (ffts.fypp).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <hip/hip_runtime_api.h>

/* hipfftXt.h is C++ only, so the few hipFFT entry points used here are
 * declared directly (values from hipfft.h / hipfftXt.h). */
typedef struct hipfftHandle_t *hipfftHandle;
typedef int hipfftResult;
typedef int hipfftXtCallbackType;
#define HIPFFT_SUCCESS 0
#define HIPFFT_D2Z 0x6a
#define HIPFFT_Z2D 0x6c
#define HIPFFT_Z2Z 0x69
#define HIPFFT_CB_LD_COMPLEX_DOUBLE 0x1
#define HIPFFT_CB_LD_REAL_DOUBLE 0x3
#define HIPFFT_CB_ST_COMPLEX_DOUBLE 0x5
hipfftResult hipfftCreate(hipfftHandle *plan);
hipfftResult hipfftDestroy(hipfftHandle plan);
hipfftResult hipfftMakePlanMany(hipfftHandle plan, int rank, int *n, int *inembed, int istride, int idist,
                                int *onembed, int ostride, int odist, int type, int batch, size_t *workSize);
hipfftResult hipfftXtSetJITCallback(hipfftHandle plan, const char *symbol_name, const void *bitcode_data,
                                    size_t bitcode_len_bytes, hipfftXtCallbackType cbtype, void **cbdata);

#define CHANNEL_INCBIN(name, path)                                                                        \
  __asm__(".section .rodata\n"                                                                           \
          ".balign 16\n"                                                                                 \
          ".globl " #name "_begin\n" #name "_begin:\n"                                                   \
          ".incbin \"" path "\"\n"                                                                       \
          ".globl " #name "_end\n" #name "_end:\n"                                                       \
          ".byte 0\n"                                                                                    \
          ".previous\n");                                                                                \
  extern const char name##_begin[];                                                                     \
  extern const char name##_end[];
CHANNEL_INCBIN(channel_xz_cb_spv, CHANNEL_XZ_CB_SPV)
CHANNEL_INCBIN(channel_xz_cb_hft_load_spv, CHANNEL_XZ_CB_HFT_LOAD_SPV)

/* Must match channel_xz_cb_info in hipfft_xz_callbacks.hip. */
typedef struct {
  unsigned int nzd, nxB, nzB, nxd1, nxpp, sendcount, nxr1, layout;
  unsigned long long a_off, b_off;
  double factor;
  void *zbuf;
  void *xbuf;
} channel_xz_cb_info;

static channel_xz_cb_info host_info;
static channel_xz_cb_info *dev_info = NULL;
static void *cbdata_array[1];
static int hft_products = 0;

/* 1 when the HFT plan forms the products itself (load callback). */
int channel_xz_cb_hft_products(void) { return hft_products; }

static int check(hipfftResult r, const char *what) {
  if (r != HIPFFT_SUCCESS) fprintf(stderr, "hipfft_xz_plans: %s failed: %d\n", what, (int)r);
  return r != HIPFFT_SUCCESS;
}

static int set_cb(hipfftHandle plan, const char *symbol, hipfftXtCallbackType type) {
  const char *begin = channel_xz_cb_spv_begin, *end = channel_xz_cb_spv_end;
  if (strcmp(symbol, "channel_cb_hft_load_products") == 0) {
    begin = channel_xz_cb_hft_load_spv_begin;
    end = channel_xz_cb_hft_load_spv_end;
  }
  return check(hipfftXtSetJITCallback(plan, symbol, begin, (size_t)(end - begin), type, cbdata_array), symbol);
}

/* Creates the four callback plans (each pointer receives a hipfftHandle, NULL
 * when its creation failed).  Returns the number of failures. */
int channel_xz_cb_create(int nzd, int nxB, int nzB, int nxd, int nyb, void **pIFT, void **pFFT, void **pRFT,
                         void **pHFT) {
  int fails = 0;
  hipfftHandle h;
  size_t ws;
  int npxz = nzd/nzB;
  int n[1], inembed[1], onembed[1];

  host_info.nzd = (unsigned)nzd;
  host_info.nxB = (unsigned)nxB;
  host_info.nzB = (unsigned)nzB;
  host_info.nxd1 = (unsigned)(nxd + 1);
  host_info.nxpp = (unsigned)(npxz*nxB);
  host_info.sendcount = (unsigned)nxB*(unsigned)nzB*(unsigned)nyb;
  host_info.a_off = 0;
  host_info.b_off = 0;
  host_info.factor = 1.0;
  host_info.zbuf = NULL;
  host_info.xbuf = NULL;
  host_info.nxr1 = (unsigned)(2*(nxd + 1));
  host_info.layout = 0;
  {
    const char *lay = getenv("CHANNEL_XZ_CB_LAYOUT");
    if (lay) host_info.layout = (unsigned)atoi(lay);
  }
  if (dev_info == NULL && hipMalloc((void **)&dev_info, sizeof(host_info)) != hipSuccess) return 100;
  if (hipMemcpy(dev_info, &host_info, sizeof(host_info), hipMemcpyHostToDevice) != hipSuccess) return 101;
  cbdata_array[0] = dev_info;

  /* z: complex, length nzd, nxB*nyb lines, contiguous */
  n[0] = nzd; inembed[0] = nzd; onembed[0] = nzd;
  *pIFT = NULL;
  if (!check(hipfftCreate(&h), "create IFT") && !set_cb(h, "channel_cb_ift_store_pack", HIPFFT_CB_ST_COMPLEX_DOUBLE) &&
      !check(hipfftMakePlanMany(h, 1, n, inembed, 1, nzd, onembed, 1, nzd, HIPFFT_Z2Z, nxB*nyb, &ws), "make IFT"))
    *pIFT = h;
  else
    fails++;
  *pFFT = NULL;
  if (!check(hipfftCreate(&h), "create FFT") && !set_cb(h, "channel_cb_fft_load_unpack", HIPFFT_CB_LD_COMPLEX_DOUBLE) &&
      !check(hipfftMakePlanMany(h, 1, n, inembed, 1, nzd, onembed, 1, nzd, HIPFFT_Z2Z, nxB*nyb, &ws), "make FFT"))
    *pFFT = h;
  else
    fails++;

  /* x: real length 2*nxd; complex leading dim nxd+1, real leading dim 2*(nxd+1) */
  n[0] = 2*nxd;
  inembed[0] = nxd + 1; onembed[0] = 2*(nxd + 1);
  *pRFT = NULL;
  if (!check(hipfftCreate(&h), "create RFT") && !set_cb(h, "channel_cb_rft_load_unpack", HIPFFT_CB_LD_COMPLEX_DOUBLE) &&
      !check(hipfftMakePlanMany(h, 1, n, inembed, 1, nxd + 1, onembed, 1, 2*(nxd + 1), HIPFFT_Z2D, nzB*nyb, &ws),
             "make RFT"))
    *pRFT = h;
  else
    fails++;
  inembed[0] = 2*(nxd + 1); onembed[0] = nxd + 1;
  /* HFT: the pack store callback, and with CHANNEL_XZ_CB_HFT=products also
   * the products load callback.  The latter is correct but measured slower
   * (small x4y1 0.0895 vs 0.0846 s/step, big 0.774 vs 0.719): the R2C kernel
   * then loads two factors per element through the callback.  Without it the
   * products come from the build_products kernel. */
  *pHFT = NULL;
  hft_products = 0;
  {
    const char *mode = getenv("CHANNEL_XZ_CB_HFT");
    int want_load = (mode && strcmp(mode, "products") == 0);
    if (want_load && !check(hipfftCreate(&h), "create HFT") &&
        !set_cb(h, "channel_cb_hft_load_products", HIPFFT_CB_LD_REAL_DOUBLE) &&
        !set_cb(h, "channel_cb_hft_store_pack", HIPFFT_CB_ST_COMPLEX_DOUBLE)) {
      if (!check(hipfftMakePlanMany(h, 1, n, inembed, 1, 2*(nxd + 1), onembed, 1, nxd + 1, HIPFFT_D2Z, nzB*nyb, &ws),
                 "make HFT (products + pack)")) {
        *pHFT = h;
        hft_products = 1;
      } else {
        hipfftDestroy(h);
      }
    }
    if (*pHFT == NULL) {
      if (!check(hipfftCreate(&h), "create HFT") &&
          !set_cb(h, "channel_cb_hft_store_pack", HIPFFT_CB_ST_COMPLEX_DOUBLE) &&
          !check(hipfftMakePlanMany(h, 1, n, inembed, 1, 2*(nxd + 1), onembed, 1, nxd + 1, HIPFFT_D2Z, nzB*nyb,
                                    &ws), "make HFT (pack)"))
        *pHFT = h;
      else
        fails++;
    }
  }
  return fails;
}

/* Points the callbacks at transpose buffers for the next execution (NULL
 * keeps the current one). */
int channel_xz_cb_set_buffers(void *zbuf, void *xbuf) {
  int changed = 0;
  if (dev_info == NULL) return 1;
  if (zbuf != NULL && zbuf != host_info.zbuf) { host_info.zbuf = zbuf; changed = 1; }
  if (xbuf != NULL && xbuf != host_info.xbuf) { host_info.xbuf = xbuf; changed = 1; }
  if (!changed) return 0;
  return hipMemcpy(dev_info, &host_info, sizeof(host_info), hipMemcpyHostToDevice) != hipSuccess;
}

/* Selects the two rVVdx components (element offsets from rVVdx(1,1,y0,1)) and
 * the factor the next HFT products callback multiplies. */
int channel_xz_cb_set_products(long long a_off, long long b_off, double factor) {
  if (dev_info == NULL) return 1;
  if ((unsigned long long)a_off == host_info.a_off && (unsigned long long)b_off == host_info.b_off &&
      factor == host_info.factor)
    return 0;
  host_info.a_off = (unsigned long long)a_off;
  host_info.b_off = (unsigned long long)b_off;
  host_info.factor = factor;
  return hipMemcpy(dev_info, &host_info, sizeof(host_info), hipMemcpyHostToDevice) != hipSuccess;
}

void channel_xz_cb_destroy(void *pIFT, void *pFFT, void *pRFT, void *pHFT) {
  if (pIFT) hipfftDestroy((hipfftHandle)pIFT);
  if (pFFT) hipfftDestroy((hipfftHandle)pFFT);
  if (pRFT) hipfftDestroy((hipfftHandle)pRFT);
  if (pHFT) hipfftDestroy((hipfftHandle)pHFT);
}
