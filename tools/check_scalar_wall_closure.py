#!/usr/bin/env python3
"""Solvability of the scalar wall-normal system under Neumann rows at both walls.

Rebuilds the compact stencils of `case_setup.fypp:setup_derivatives` on a given
mesh, assembles the per-substep scalar system exactly as
`assemble_compact_helmholtz_system` plus the four closure rows of
`channel_bcs.f90` do, and reports:

  1. the four pivots the wall closure divides by (wall_closure.fypph:48/56,
     compact_component_solve.fypp:364/386, compact_line_solvers.fypp:162),
  2. rank and condition number at k = 0 and at the extreme modes,
  3. whether a solve actually carries the prescribed wall gradients,
  4. the integral of the mean-correction field, which `linsolve_scalar`
     divides by to form corrtx.

The short answer it exists to support: with lambda = RK_rai(1,i)/deltat > 0 the
operator is a Helmholtz, not a Laplacian, so Neumann-at-both-walls is
nonsingular.  The familiar pure-Neumann null space appears only at lambda = 0,
which the scalars never solve.  Run with --lambda 0 to see it.

Usage:  python3 tools/check_scalar_wall_closure.py [--ny 500] [--stretching 1.66]
                                                   [--ni 20000] [--pr 0.025]
                                                   [--deltat 4e-3] [--lambda L]
"""
import argparse

import numpy as np
from numpy.linalg import cond, matrix_rank, solve, svd


def build(ny, a, ymin, ymax):
    """The mesh and the compact rows, as setup_derivatives computes them."""
    idx = np.arange(-1, ny + 2)
    Y = ymin + 0.5 * (ymax - ymin) * (np.tanh(a * (2.0 * idx / ny - 1)) / np.tanh(a) + 1)
    y = lambda i: Y[i + 1]

    der = np.zeros((ny + 1, 4, 5))
    for iy in range(1, ny):
        d = np.array([y(n) - y(iy) for n in range(iy - 2, iy + 3)])
        M = np.array([[d[j] ** (4.0 - i) for j in range(5)] for i in range(5)])
        t = np.zeros(5)
        t[0] = 24.0
        der[iy, 3] = np.linalg.solve(M, t)                                  # D^4
        M2 = np.array([[(5.0 - i) * (6.0 - i) * (7.0 - i) * (8.0 - i) * d[j] ** (4.0 - i)
                        for j in range(5)] for i in range(5)])
        t = np.array([np.sum(der[iy, 3] * d ** (8.0 - i)) for i in range(5)])
        der[iy, 0] = np.linalg.solve(M2, t)                                 # D^0
        t = np.zeros(5)
        for i in range(3):
            t[i] = np.sum(der[iy, 0] * (4.0 - i) * (3.0 - i) * d ** (2.0 - i))
        der[iy, 2] = np.linalg.solve(M, t)                                  # D^2

    def one_sided(nodes, anchor, order):
        d = np.array([y(n) - anchor for n in nodes])
        M = np.array([[d[j] ** (4.0 - i) for j in range(5)] for i in range(5)])
        t = np.zeros(5)
        t[4 - order] = 1.0
        return np.linalg.solve(M, t)

    rows = dict(
        d140=one_sided([-1, 0, 1, 2, 3], y(0), 1),
        d040=np.array([0, 1, 0, 0, 0], float),
        d14n=one_sided([ny - 3, ny - 2, ny - 1, ny, ny + 1], y(ny), 1),
        d04n=np.array([0, 0, 0, 1, 0], float),
    )
    return Y, der, rows


def assemble(ny, der, lam, alpha, k2, lower, upper):
    """Square system on nodes -1 .. ny+1: two closure rows, the interior, two more."""
    n = ny + 3
    A = np.zeros((n, n))
    col = lambda node: node + 1
    (wl, gl), (wu, gu) = lower, upper
    for j, c in enumerate(range(-1, 4)):
        A[0, col(c)] = gl[j]
        A[1, col(c)] = wl[j]
    for iy in range(1, ny):
        for j, c in enumerate(range(iy - 2, iy + 3)):
            A[iy + 1, col(c)] = lam * der[iy, 0, j] - alpha * (der[iy, 2, j] - k2 * der[iy, 0, j])
    for j, c in enumerate(range(ny - 3, ny + 2)):
        A[ny + 1, col(c)] = wu[j]
        A[ny + 2, col(c)] = gu[j]
    return A


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ny", type=int, default=500)
    p.add_argument("--stretching", type=float, default=1.66)
    p.add_argument("--ymin", type=float, default=0.0)
    p.add_argument("--ymax", type=float, default=2.0)
    p.add_argument("--ni", type=float, default=20000.0, help="as in the deck, i.e. Re")
    p.add_argument("--pr", type=float, default=0.025)
    p.add_argument("--deltat", type=float, default=4.0e-3)
    p.add_argument("--lambda", dest="lam", type=float, default=None,
                   help="override lambda; 0 exposes the pure-Neumann null space")
    args = p.parse_args()

    ny = args.ny
    Y, der, r = build(ny, args.stretching, args.ymin, args.ymax)
    alpha = (1.0 / args.ni) / args.pr
    lam = args.lam if args.lam is not None else 3.75 / args.deltat  # min RK_rai(1,i)/deltat

    lower = {"dirichlet": (r["d040"], der[1, 3]), "neumann": (r["d140"], der[1, 3])}
    upper = {"dirichlet": (r["d04n"], der[ny - 1, 3]),
             "neumann (fixed)": (r["d14n"], der[ny - 1, 3]),
             "neumann (in tree)": (r["d14n"], r["d04n"])}

    print(f"mesh ny={ny} stretching={args.stretching}  dy_1={Y[2]-Y[1]:.4e}")
    print(f"alpha={alpha:.4e}  lambda={lam:.4e}\n")

    print("closure pivots -- every one of these is a divisor in the solver")
    for k, (w, g) in lower.items():
        fold = w[1:] - g[1:] * w[0] / g[0]
        print(f"  lower {k:18s} eqm1(-2)={g[0]:+.4e}  folded eq0(-1)={fold[0]:+.4e}")
    for k, (w, g) in upper.items():
        if g[4] == 0.0:
            print(f"  upper {k:18s} eqnp1(+2)={g[4]:+.4e}  folded eqn(+1)=DIVISION BY ZERO")
        else:
            fold = w[:4] - g[:4] * w[4] / g[4]
            print(f"  upper {k:18s} eqnp1(+2)={g[4]:+.4e}  folded eqn(+1)={fold[3]:+.4e}")

    eqlb = lambda A: A / np.abs(A).max(axis=1, keepdims=True)
    print("\nconditioning at k=0 (raw cond is dominated by the D^4 ghost-row scaling;")
    print("the row-equilibrated number is the meaningful one)")
    for nm, lo, up in [("Dirichlet/Dirichlet", lower["dirichlet"], upper["dirichlet"]),
                       ("Neumann/Neumann    ", lower["neumann"], upper["neumann (fixed)"])]:
        A = assemble(ny, der, lam, alpha, 0.0, lo, up)
        # rank on the equilibrated matrix: matrix_rank's tolerance is relative to
        # the largest singular value, which the 1e12-scale D^4 ghost rows dominate.
        print(f"  {nm}  rank={matrix_rank(eqlb(A))}/{ny+3}  raw={cond(A):.2e}  equilibrated={cond(eqlb(A)):.2e}")

    if lam == 0.0:
        A = assemble(ny, der, 0.0, alpha, 0.0, lower["neumann"], upper["neumann (fixed)"])
        s = svd(eqlb(A), compute_uv=False)
        print(f"\n  lambda=0: smallest singular values {s[-3:]} -- one null direction, the constant")
        return

    lo, up = lower["neumann"], upper["neumann (fixed)"]
    for k2, label in [(0.0, "k=0"), (0.0625, "k=alfa0"), (262144.0, "k=k_max")]:
        A = assemble(ny, der, lam, alpha, k2, lo, up)
        print(f"  Neumann/Neumann  {label:8s} rank={matrix_rank(eqlb(A))}/{ny+3}")

    g_bot, g_top = +37.0, -37.0
    A = assemble(ny, der, lam, alpha, 0.0, lo, up)
    b = np.zeros(ny + 3)
    b[1], b[ny + 1] = g_bot, g_top
    b[2:ny + 1] = 1.0
    phi = solve(A, b)
    print("\nrecovered wall gradients (the row is only real if these come back)")
    print(f"  bottom asked {g_bot:+.6f} got {np.dot(r['d140'], phi[0:5]):+.15f}")
    print(f"  top    asked {g_top:+.6f} got {np.dot(r['d14n'], phi[ny-2:ny+3]):+.15f}")

    tcor = solve(A, np.concatenate(([0, 0], np.ones(ny - 1), [0, 0])))
    w = np.zeros(ny + 3)
    w[1:ny + 2] = np.concatenate(([(Y[2] - Y[1]) / 2],
                                  (Y[3:ny + 2] - Y[1:ny]) / 2,
                                  [(Y[ny + 1] - Y[ny]) / 2]))
    print(f"\nmean correction: integral of tcor = {np.dot(w, tcor):.6e}"
          f"  (corrtx divides by this; ~2/lambda = {2/lam:.3e})")


if __name__ == "__main__":
    main()
