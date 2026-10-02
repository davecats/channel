#!/usr/bin/env python3
"""Build a doubled-resolution restart from an existing one, memory-aware.

Source: nx, ny, nz            Target: 2*nx+1, 2*ny, 2*nz
so 1023x500x512 -> 2047x1000x1024, the same physical domain sampled twice as
finely in every direction.

x and z are Fourier and the field is dealiased well below the transform
Nyquist (nxd = nzd = 1536 for nx = nz-modes = 1023/512), so extending the mode
range is exact zero-padding -- no folding to undo.  y is a stretched tanh mesh
in physical space and needs real interpolation; the doubled mesh contains the
source mesh exactly at its even indices, because

    y(iy) = ymin + (ymax-ymin)/2 * (tanh(a*(2*iy/ny - 1))/tanh(a) + 1)

depends on iy/ny, and (2k)/(2*ny) = k/ny.  So half the target points are
copies and only the odd ones are interpolated, by 4-point Lagrange.

MEMORY.  The target field is 376 GiB and is never held.  The file is laid out
Fortran-order (y, kz, kx, component), so one (component, kx) plane is
contiguous: the loop reads one source plane (8 MB), interpolates it (16 MB)
and writes one target plane (33 MB).

The source is read with explicit seeks rather than a memmap.  A memmap works
and is tidier, but every page touched stays resident as clean page cache, so
the first version of this measured 47.5 GB RSS on a 1023x500x512 source -- the
whole file.  Reclaimable, and harmless on a 234 GB node, but it makes the
resident set scale with the input instead of with the plane, which is the
opposite of the point.  With seeks it stays at a few tens of MB.
"""
import argparse, sys, time
import numpy as np


def y_grid(ny, a, ymin, ymax):
    """Node positions iy = -1 .. ny+1, matching case_setup.fypp:178."""
    iy = np.arange(-1, ny + 2, dtype=float)
    return ymin + 0.5 * (ymax - ymin) * (np.tanh(a * (2.0 * iy / ny - 1.0)) / np.tanh(a) + 1.0)


def lagrange4_table(ys, yt):
    """For each target node, four source indices and their Lagrange weights.

    Returns (idx, w) of shape (len(yt), 4).  At a target node that coincides
    with a source node the weights come out exactly (0,1,0,0) -- Lagrange is
    interpolatory -- which is asserted by the caller.
    """
    n = len(ys)
    j = np.clip(np.searchsorted(ys, yt) - 1, 1, n - 3)
    idx = j[:, None] + np.array([-1, 0, 1, 2])
    xs = ys[idx]                                  # (nt, 4)
    w = np.ones_like(xs)
    for k in range(4):
        for m in range(4):
            if k != m:
                w[:, k] *= (yt - xs[:, m]) / (xs[:, k] - xs[:, m])
    return idx, w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source"); ap.add_argument("target")
    ap.add_argument("--ni-new", type=float, required=True,
                    help="kinematic viscosity of the NEW case (the header must match its dns.in)")
    ap.add_argument("--fluct-scale", type=float, default=1.0,
                    help="multiply every mode except (kx=0,kz=0) by this")
    ap.add_argument("--mean-scale", type=float, default=1.0,
                    help="multiply the (kx=0,kz=0) mean profiles by this")
    ap.add_argument("--time", type=float, default=0.0,
                    help="time written to the new header (default 0: a fresh clock)")
    args = ap.parse_args()

    with open(args.source, "rb") as f:
        nx, ny, nz = np.fromfile(f, np.int32, 3)
        alfa0, beta0, ni, a, ymin, ymax, t_src = np.fromfile(f, np.float64, 7)
    nx, ny, nz = int(nx), int(ny), int(nz)
    HDR = 3 * 4 + 7 * 8
    nyf, nzf, nxf = ny + 3, 2 * nz + 1, nx + 1
    import os
    ncomp = round((os.path.getsize(args.source) - HDR) / (16.0 * nyf * nzf * nxf))

    NX, NY, NZ = 2 * nx + 1, 2 * ny, 2 * nz
    NYF, NZF, NXF = NY + 3, 2 * NZ + 1, NX + 1
    print(f"source {nx}x{ny}x{nz} ncomp={ncomp} 1/ni={1/ni:.0f} a={a} t={t_src:.4f}")
    print(f"target {NX}x{NY}x{NZ} ncomp={ncomp} 1/ni={1/args.ni_new:.0f} a={a} t={args.time:.4f}")
    print(f"scaling: mean x{args.mean_scale:.6f}   fluctuations x{args.fluct_scale:.6f}")
    print(f"target file: {(HDR + 16*NYF*NZF*NXF*ncomp)/2**30:.1f} GiB")

    ys, yt = y_grid(ny, a, ymin, ymax), y_grid(NY, a, ymin, ymax)
    assert ys[0] <= yt[0] and yt[-1] <= ys[-1], "target ghost nodes fall outside the source range"
    idx, w = lagrange4_table(ys, yt)
    # Index it holds iy = it-1, so y_t(iy_t) == y_s(iy_s) when iy_t = 2*iy_s, i.e.
    # it = 2*is - 1: the ODD target indices are the interior source nodes.  The
    # source's own ghost rows have no counterpart and are only interpolated from.
    err = np.abs((w * ys[idx]).sum(1) - yt).max()
    coinc = np.abs((w[1::2] * ys[idx[1::2]]).sum(1) - ys[1:-1]).max()
    print(f"y table: reproduces its own abscissae to {err:.2e}; "
          f"odd target nodes match the interior source nodes to {coinc:.2e}")
    assert coinc < 1e-12
    assert err < 1e-12

    plane_elems = nyf * nzf
    plane_bytes = plane_elems * 16
    fsrc = open(args.source, "rb")

    def read_source_plane(c, ix):
        """One (component, kx) plane, (nyf, nzf), Fortran order on disk."""
        fsrc.seek(HDR + plane_bytes * (ix + nxf * c))
        return np.fromfile(fsrc, np.complex128, plane_elems).reshape((nyf, nzf), order="F")

    with open(args.target, "wb") as out:
        out.write(np.array([NX, NY, NZ], dtype=np.int32).tobytes())
        out.write(np.array([alfa0, beta0, args.ni_new, a, ymin, ymax, args.time],
                           dtype=np.float64).tobytes())
        kz_off = (NZF - nzf) // 2                       # same kz -> shifted index
        plane = np.empty((NYF, NZF), dtype=np.complex128, order="F")
        t0 = time.time(); done = 0; total = ncomp * NXF
        for c in range(ncomp):
            for ix in range(NXF):
                plane[:] = 0.0
                if ix < nxf:
                    s = read_source_plane(c, ix)                       # (nyf, nzf)
                    # y-interpolate: (NYF,4) weights against gathered source rows
                    interp = np.einsum('ij,ijk->ik', w, s[idx, :])     # (NYF, nzf)
                    # exact at the coincident nodes, so no interpolation error is
                    # introduced on half the mesh
                    interp *= args.fluct_scale
                    if ix == 0:
                        interp[:, nz] = np.einsum('ij,ij->i', w, s[idx, nz]) * args.mean_scale
                    plane[:, kz_off:kz_off + nzf] = interp
                out.write(np.asfortranarray(plane).tobytes(order="F"))
                done += 1
                if done % 512 == 0 or done == total:
                    el = time.time() - t0
                    print(f"  {done}/{total} planes  {el:6.0f}s  eta {el*(total/done-1):6.0f}s", flush=True)
    fsrc.close()
    print("done")


if __name__ == "__main__":
    main()
