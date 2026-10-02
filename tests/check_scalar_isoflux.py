#!/usr/bin/env python3
"""Checks the isoflux scalar wall condition in a whole run, next to an isothermal one.

One run, two scalars at the same Prandtl number over the same velocity field,
`bc = dirichlet neumann`.  Three properties, each of which the unit test on the
line solver cannot see:

  1. the Neumann scalar's wall gradient is the value the deck asked for, at
     every step and at both walls, to round-off -- the condition is enforced on
     the instantaneous field and survives the explicit half of the time step
  2. the Dirichlet scalar's wall gradient is free and does move, so the two
     scalars really took different rows in the same run
  3. the global balance closes: d/dt INT(phi) dy = alpha*(phi'(2) - phi'(0)),
     which is what says the imposed flux is the only flux

Property 3 is the one that catches a wall condition that is enforced but wrong,
e.g. a sign or a factor in the row.  It is checked to 1%, which is the spatial
and temporal truncation error of the deck below; refining either tightens it.

`meantb = 0`: the bulk-pinning correction would add to INT(phi) deliberately,
and property 3 is about what the wall does.
"""

import os
import shutil
import subprocess
import sys
import tempfile

NI = 1000.0          # the deck's [velocity] ni, i.e. 1/ni in code units
PR = 1.0
ALPHA = (1.0 / NI) / PR
G = 1.0              # imposed dphi/dy: +G at the lower wall, -G at the upper

DECK = f"""\
[mesh]
nx = 7
ny = 64
nz = 4
alfa0 = 1.0
beta0 = 2.0
stretching = 1.6
ymin = 0.0
ymax = 2.0

[velocity]
ni = {NI}
meanpx = 0.0
meanpz = 0.0
meanflowx = 2.0
meanflowz = 0.0
u0 = 0.0
uN = 0.0
perturbation_amplitude = 1.0e-3
seed = 7

[scalars]
nPhi = 2
meantx = 0.0
meantb = 0.0
pr = {PR} {PR}
bc = dirichlet neumann
t0 = 0.0 {G}
tN = 0.0 {-G}

[timestepping]
deltat = 0.01
cflmax = 0.0
time = 0.0
dt_field = 1000.0
dt_save = -1.0
t_max = 1000.0
time_from_restart = false
nstep = 200
"""


def run(channel, mpiexec):
    """Run channel from a seeded cold start; return the rows of Runtimedata.phi."""
    workdir = tempfile.mkdtemp(prefix="scalar_isoflux_")
    try:
        with open(os.path.join(workdir, "dns.in"), "w") as fh:
            fh.write(DECK)

        proc = subprocess.run(
            [mpiexec, "-np", "1", channel],
            cwd=workdir, capture_output=True, text=True, timeout=1200,
        )
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout + proc.stderr)
            raise SystemExit(f"channel failed (exit {proc.returncode})")

        path = os.path.join(workdir, "Runtimedata.phi")
        if not os.path.exists(path):
            sys.stderr.write(proc.stdout + proc.stderr)
            raise SystemExit("no Runtimedata.phi written")
        with open(path) as fh:
            rows = [[float(w) for w in line.split()] for line in fh if line.strip()]
        # Columns: time | dphi/dy lower (nPhi) | upper (nPhi) | bulk (nPhi) |
        # corrtx (nPhi).  The first row is the start field, written before the
        # first solve, so its scalar still carries the initial profile's gradient.
        if len(rows) < 3:
            raise SystemExit(f"Runtimedata.phi has only {len(rows)} rows")
        return rows[1:]
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_scalar_isoflux.py <channel-binary> <mpiexec>")
    channel, mpiexec = os.path.abspath(sys.argv[1]), sys.argv[2]

    rows = run(channel, mpiexec)
    # nPhi = 2, so: time, g0(1), g0(2), gN(1), gN(2), bulk(1), bulk(2), ...
    time = [r[0] for r in rows]
    dirichlet_g0, neumann_g0 = [r[1] for r in rows], [r[2] for r in rows]
    dirichlet_gN, neumann_gN = [r[3] for r in rows], [r[4] for r in rows]
    neumann_bulk = [r[6] for r in rows]

    failures = []

    worst = max(max(abs(g - G) for g in neumann_g0),
                max(abs(g + G) for g in neumann_gN))
    if not worst <= 1.0e-12:
        failures.append(f"Neumann scalar: wall gradient off the imposed value by {worst:.3e} "
                        f"(want {G} at y = 0 and {-G} at y = 2)")

    spread = max(max(dirichlet_g0) - min(dirichlet_g0),
                 max(dirichlet_gN) - min(dirichlet_gN))
    if not spread > 1.0e-8:
        failures.append(f"Dirichlet scalar: wall gradient varies by only {spread:.3e} over the run "
                        "(both scalars appear to be taking the same rows)")

    measured = (neumann_bulk[-1] - neumann_bulk[0]) / (time[-1] - time[0])
    expected = ALPHA * (-G - G)
    if not abs(measured - expected) <= 0.01 * abs(expected):
        failures.append(f"Neumann scalar: d/dt INT(phi) dy = {measured:.6e}, "
                        f"but alpha*(phi'(2) - phi'(0)) = {expected:.6e}")

    if failures:
        for f in failures:
            print("FAIL:", f)
        raise SystemExit(1)
    print(f"scalar isoflux: wall gradient held to {worst:.2e}, "
          f"global balance closes to {abs(measured / expected - 1.0):.2%}, "
          f"isothermal scalar's gradient free over {spread:.2e}")


if __name__ == "__main__":
    main()
