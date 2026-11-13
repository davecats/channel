// hipFFT/rocFFT JIT load/store callbacks for the xz transposes.
//
// Compiled to SPIR-V (--offload-arch=amdgcnspirv --offload-device-only) and
// handed to hipfftXtSetJITCallback by hipfft_xz_plans.c; rocFFT links them
// into its runtime-compiled kernels, so they are inlined into the FFT.
//
// The transpose buffers are not the plans' input/output arrays: the plans run
// on VVdz / VVdx as before (so rocFFT never writes outside them) and the
// callbacks read or write zbuf / xbuf from the callback data instead, which
// ffts.fypp points at the right buffer slot before each execution.
//
// The callbacks fold the alltoall pack/unpack kernels (and the build_products
// kernel) into the FFTs.  Every callback gets `offset` in the layout the plan
// was made with (the VVdz layout for the z plans, VVdx / rVVdx for the x
// plans) and computes where that element lives in the transpose buffer.  The
// buffer layouts are those of pack_*/unpack_* in mpi_transpose.f90 (0-based
// here):
//   zTOx block p: p = dest*sendcount + iz + nzB*(ix + nxB*iy)
//   xTOz block p: p = dest*sendcount + ix + nxB*(iz + nzB*iy)
// CHANNEL_XZ_CB_LAYOUT (bit 1: zTOx, bit 2: xTOz) flips a block to the other
// fastest index, moving the strided access from the load to the store side;
// pack and unpack are both callbacks, so they always agree.
// When the transpose is local (npxz = 1) the zTOx "buffer" is VVdz itself and
// the xTOz one is a compact copy of VVdx, so the same formulas apply.
//
// The HFT plan carries two callbacks and rocFFT links each one's bitcode
// into the same kernel, so the products callback is compiled into a blob of
// its own (CHANNEL_CB_HFT_LOAD_ONLY) to keep the symbols unique.
//
// All index arithmetic is 32-bit unsigned: init_MPI guarantees that one
// buffer slot holds fewer than 2^31 elements.
#include <hip/hip_runtime.h>

struct channel_xz_cb_info {
  unsigned int nzd, nxB, nzB, nxd1, nxpp, sendcount, nxr1, layout;
  unsigned long long a_off, b_off;  // HFT products: element offsets of the two factors from rVVdx(1,1,y0,1)
  double factor;                    // HFT products: normalisation
  double2* zbuf;                    // zTOx send (IFT) / xTOz receive (FFT) buffer
  double2* xbuf;                    // zTOx receive (RFT) / xTOz send (HFT) buffer
};

#ifdef CHANNEL_CB_HFT_LOAD_ONLY
// HFT load: the product of two velocity / scalar components, computed from
// rVVdx on the fly (buf is rVVdx(1, 1, y0, 1)).
extern "C" __device__ double channel_cb_hft_load_products(double* buf, size_t offset, void* cbdata, void* shared) {
  const channel_xz_cb_info* c = (const channel_xz_cb_info*)cbdata;
  unsigned int o = (unsigned int)offset;
  double a = buf[(size_t)o + c->a_off];
  double b = buf[(size_t)o + c->b_off];
  return a*b*c->factor;
}

#else
// IFT store: z-pencil element (iz, ix, iy) of VVdz goes straight into the
// zTOx send block of the rank owning iz.
extern "C" __device__ void channel_cb_ift_store_pack(double2* buf, size_t offset, double2 e, void* cbdata,
                                                     void* shared) {
  const channel_xz_cb_info* c = (const channel_xz_cb_info*)cbdata;
  unsigned int o = (unsigned int)offset;
  unsigned int iz = o % c->nzd;
  unsigned int r = o / c->nzd;
  unsigned int ix = r % c->nxB;
  unsigned int iy = r / c->nxB;
  unsigned int dest = iz / c->nzB;
  unsigned int izl = iz - dest*c->nzB;
  if (c->layout & 1)  // zTOx blocks x-fastest: strided store here, contiguous load in the RFT
    c->zbuf[dest*c->sendcount + ix + c->nxB*(izl + c->nzB*iy)] = e;
  else
    c->zbuf[dest*c->sendcount + izl + c->nzB*(ix + c->nxB*iy)] = e;
}

// RFT load: x-pencil element (ix, iz, iy) of VVdx read from the zTOx receive
// buffer; the x padding (ix >= nx+1) is zero.
extern "C" __device__ double2 channel_cb_rft_load_unpack(double2* buf, size_t offset, void* cbdata, void* shared) {
  const channel_xz_cb_info* c = (const channel_xz_cb_info*)cbdata;
  unsigned int o = (unsigned int)offset;
  unsigned int ixf = o % c->nxd1;
  unsigned int r = o / c->nxd1;
  unsigned int iz = r % c->nzB;
  unsigned int iy = r / c->nzB;
  if (ixf >= c->nxpp) return make_double2(0.0, 0.0);
  unsigned int src = ixf / c->nxB;
  unsigned int ix = ixf - src*c->nxB;
  if (c->layout & 1) return c->xbuf[src*c->sendcount + ix + c->nxB*(iz + c->nzB*iy)];
  return c->xbuf[src*c->sendcount + iz + c->nzB*(ix + c->nxB*iy)];
}

// HFT store: x-pencil element (ix, iz, iy) goes into the xTOz send block of
// the rank owning ix; modes ix >= nx+1 are not sent.
extern "C" __device__ void channel_cb_hft_store_pack(double2* buf, size_t offset, double2 e, void* cbdata,
                                                     void* shared) {
  const channel_xz_cb_info* c = (const channel_xz_cb_info*)cbdata;
  unsigned int o = (unsigned int)offset;
  unsigned int ixf = o % c->nxd1;
  unsigned int r = o / c->nxd1;
  unsigned int iz = r % c->nzB;
  unsigned int iy = r / c->nzB;
  if (ixf >= c->nxpp) return;
  unsigned int dest = ixf / c->nxB;
  unsigned int ix = ixf - dest*c->nxB;
  if (c->layout & 2)  // xTOz blocks z-fastest: strided store here, contiguous load in the FFT
    c->xbuf[dest*c->sendcount + iz + c->nzB*(ix + c->nxB*iy)] = e;
  else
    c->xbuf[dest*c->sendcount + ix + c->nxB*(iz + c->nzB*iy)] = e;
}

// FFT load: z-pencil element (iz, ix, iy) of VVdz read from the xTOz receive
// buffer.
extern "C" __device__ double2 channel_cb_fft_load_unpack(double2* buf, size_t offset, void* cbdata, void* shared) {
  const channel_xz_cb_info* c = (const channel_xz_cb_info*)cbdata;
  unsigned int o = (unsigned int)offset;
  unsigned int iz = o % c->nzd;
  unsigned int r = o / c->nzd;
  unsigned int ix = r % c->nxB;
  unsigned int iy = r / c->nxB;
  unsigned int src = iz / c->nzB;
  unsigned int izl = iz - src*c->nzB;
  if (c->layout & 2) return c->zbuf[src*c->sendcount + izl + c->nzB*(ix + c->nxB*iy)];
  return c->zbuf[src*c->sendcount + ix + c->nxB*(izl + c->nzB*iy)];
}
#endif
