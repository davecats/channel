#!/usr/bin/env python3
"""Mean scalar profiles out of restart snapshots -- the cheap way to watch a run settle.

A restart holds the field in Fourier space, `(ny+3, 2*nz+1, nx+1, 3+nPhi)`
complex doubles in Fortran order, so the plane-averaged profile of scalar i is
one contiguous y-line at `(kx = 0, kz = 0)`: **503 complex numbers out of a 50 GB
file**, at a computable offset. Reading every snapshot of a production case costs
a few hundred kilobytes and a second.

What it is for. Under an isoflux (`bc = neumann`) wall the mean wall gradient is
pinned by the deck and the bulk is pinned by `meantb`, so neither says anything
about whether the mean has settled -- and neither does `corrtx`, which the wall
flux alone determines (`corrtx = 2*alpha*g` to four figures in the isothermal
campaign). What is free is the *shape* of the profile, and the wall value is its
most visible number. `Runtimedata.phi` reports the wall value every step but
nothing about the shape; this reports the shape.

The drift column is the one to read: `|d<phi>/dt|` between consecutive snapshots,
normalised by the wall-to-bulk difference, in units of inverse time. It falls
towards the level set by turbulent sampling noise; that is the run settling.

    ./scalar_mean_profiles.py Dati.cart.*.out
    ./scalar_mean_profiles.py --csv profiles.csv Dati.cart.*.out
    ./scalar_mean_profiles.py --profile 3 Dati.cart.33.out   # print <phi>(y) itself
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

HEADER_BYTES = 3 * np.dtype(np.int32).itemsize + 7 * np.dtype(np.float64).itemsize
COMPLEX_BYTES = np.dtype(np.complex128).itemsize


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", type=Path, nargs="+", help="restart snapshots (Dati.cart.*.out)")
    parser.add_argument("--csv", type=Path, default=None, help="also write one row per snapshot here")
    parser.add_argument("--profile", type=int, default=None, metavar="N",
                        help="print <phi>(y) for scalar N (1-based) of each file instead of the summary")
    return parser.parse_args()


def read_header(path: Path) -> tuple[dict[str, float], int]:
    with path.open("rb") as handle:
        ints = np.fromfile(handle, dtype=np.int32, count=3)
        reals = np.fromfile(handle, dtype=np.float64, count=7)
    if ints.size != 3 or reals.size != 7:
        raise SystemExit(f"{path}: too short to hold a restart header")
    h = dict(zip(("nx", "ny", "nz"), (int(v) for v in ints)))
    h.update(zip(("alfa0", "beta0", "ni", "a", "ymin", "ymax", "time"), (float(v) for v in reals)))
    elements = (h["ny"] + 3) * (2 * h["nz"] + 1) * (h["nx"] + 1)
    payload = path.stat().st_size - HEADER_BYTES
    if payload <= 0 or payload % (elements * COMPLEX_BYTES) != 0:
        raise SystemExit(f"{path}: payload is not a whole number of components")
    return h, payload // (elements * COMPLEX_BYTES)


def mean_line(path: Path, h: dict[str, float], component: int) -> np.ndarray:
    """The (kx = 0, kz = 0) y-line of one component, as nodal values on y(-1:ny+1)."""
    nyf, nzf = h["ny"] + 3, 2 * h["nz"] + 1
    # Fortran order (iy, iz, ix, component): the line at iz = nz (kz = 0), ix = 0.
    offset = HEADER_BYTES + ((component * (h["nx"] + 1) + 0) * nzf + h["nz"]) * nyf * COMPLEX_BYTES
    with path.open("rb") as handle:
        handle.seek(offset)
        line = np.fromfile(handle, dtype=np.complex128, count=nyf)
    if line.size != nyf:
        raise SystemExit(f"{path}: short read on the mean line of component {component}")
    # The solver carries the mean mode's value in the real part; see yintegr and
    # the at_*_wall_real macros, which both take dreal().
    return line.real


def y_grid(h: dict[str, float]) -> np.ndarray:
    iy = np.arange(-1, h["ny"] + 2, dtype=np.float64)
    return h["ymin"] + 0.5 * (h["ymax"] - h["ymin"]) * (
        np.tanh(h["a"] * (2.0 * iy / h["ny"] - 1.0)) / np.tanh(h["a"]) + 1.0)


def integrate(f: np.ndarray, y: np.ndarray, ny: int) -> float:
    """yintegr from channel_equations.fypp: three-point Newton-Cotes over node pairs."""
    total = 0.0
    for iy in range(1, ny + 1, 2):
        yp1 = y[iy + 2] - y[iy + 1]          # the arrays are offset by one: index 0 is y(-1)
        ym1 = y[iy] - y[iy + 1]
        a1 = -ym1 / 3.0 + yp1 / 6.0 + yp1 * yp1 / (6.0 * ym1)
        a3 = +yp1 / 3.0 - ym1 / 6.0 - ym1 * ym1 / (6.0 * yp1)
        a2 = yp1 - ym1 - a1 - a3
        total += a1 * f[iy] + a2 * f[iy + 1] + a3 * f[iy + 2]
    return total


def wall_gradient(f: np.ndarray, y: np.ndarray, lower: bool, ny: int) -> float:
    """d<phi>/dy at a wall, from the five nodes the solver's own row uses."""
    if lower:
        nodes, at = slice(0, 5), y[1]            # y(-1 .. 3), evaluated at y(0)
    else:
        nodes, at = slice(ny - 2, ny + 3), y[ny + 1]   # y(ny-3 .. ny+1), at y(ny)
    # Exact first derivative of the quartic through those five points, which is
    # what d140 and d14n are.
    return float(np.polyval(np.polyder(np.polyfit(y[nodes] - at, f[nodes], 4)), 0.0))


def main() -> None:
    args = parse_args()
    files = sorted(args.files, key=lambda p: read_header(p)[0]["time"])

    rows = []
    for path in files:
        h, n_components = read_header(path)
        nphi = n_components - 3
        y = y_grid(h)
        profiles = [mean_line(path, h, 3 + i) for i in range(nphi)]
        rows.append((h, profiles, path))

        if args.profile is not None:
            i = args.profile - 1
            if not 0 <= i < nphi:
                raise SystemExit(f"--profile {args.profile}: the file has {nphi} scalar(s)")
            print(f"# {path.name}  t = {h['time']:.4f}  scalar {args.profile}")
            print("# y  <phi>")
            for yy, pp in zip(y, profiles[i]):
                print(f"{yy:.10e} {pp:.10e}")
            continue

    if args.profile is not None:
        return

    nphi = len(rows[0][1])
    print(f"{len(rows)} snapshot(s), nPhi = {nphi}")
    header = f"{'time':>12} {'dt':>9}"
    for i in range(nphi):
        header += f" | {'<phi>_w0':>10} {'<phi>_wN':>10} {'INT phi':>9} {'g0':>9} {'|drift|':>9}"
    print(header)

    if args.csv:
        with args.csv.open("w") as fh:
            fh.write(",".join(["time", "dt"] + [f"{name}_{i + 1}" for i in range(nphi)
                                                for name in ("w0", "wN", "int", "g0", "drift")]) + "\n")

    prev = None
    for h, profiles, path in rows:
        ny = h["ny"]
        y = y_grid(h)
        dt = h["time"] - prev[0]["time"] if prev else float("nan")
        line = f"{h['time']:12.4f} {dt:9.2f}"
        csv_row = [h["time"], dt]
        for i, f in enumerate(profiles):
            w0, wn = f[1], f[ny + 1]
            bulk = integrate(f, y, int(ny))
            g0 = wall_gradient(f, y, True, int(ny))
            if prev is not None:
                scale = max(abs(w0 - bulk / (h["ymax"] - h["ymin"])), 1e-30)
                drift = np.max(np.abs(f - prev[1][i])) / (abs(dt) * scale)
            else:
                drift = float("nan")
            line += f" | {w0:10.5f} {wn:10.5f} {bulk:9.5f} {g0:9.4f} {drift:9.2e}"
            csv_row += [w0, wn, bulk, g0, drift]
        print(line)
        prev = (h, profiles)
        if args.csv:
            with args.csv.open("a") as fh:
                fh.write(",".join(f"{v:.10e}" for v in csv_row) + "\n")

    print("\n|drift| is max|d<phi>/dt| over the profile, divided by the wall-to-bulk")
    print("difference -- an inverse time.  It falling towards a floor is the run settling;")
    print("the floor is turbulent sampling noise, not a residual transient.")


if __name__ == "__main__":
    main()
