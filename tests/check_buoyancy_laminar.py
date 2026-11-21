#!/usr/bin/env python3
"""Checks the active scalar (Boussinesq buoyancy along x, [scalars] ri) and the
velocity-weighted sink ([scalars] sink = velocity) against the laminar fully
developed solution, and the Cf / Nu / mixing-cup diagnostics with it.

Laminar mixed convection between Neumann walls at a fixed flow rate Q has an
x-independent steady state that solves two coupled ODEs,

    nu U'' + ri phi + P = 0,      U(0) = U(L) = 0,       INT U dy = Q
    kappa phi'' = gamma U,        phi'(0) = t0, phi'(L) = tN,
    gamma = kappa (tN - t0)/Q,

which are solved here by Chebyshev collocation and compared with the end of a
run started from the generated field.  Two cases:

  1. both walls heated (t0 = -1, tN = +1), aiding buoyancy ri = +0.5
  2. lower wall heated, upper wall insulated (tN = 0), opposing buoyancy
     ri = -0.2

Compared: dU/dy at both walls and Cf, the wall values and the mixing-cup bulk
of phi, and Nu at both walls (0 at the insulated one), to 1e-4 relative --
the mesh and time-step error of the deck is 1e-6.  Also: the bulk INT phi dy
must stay constant, since the sink balances the wall heat input exactly (to
1e-4 of that input; the discrete residue is a few 1e-6).
At ri = 0 the two Nusselt numbers are the textbook 140/17 and 70/13; the
reference solver reproduces both, so the definitions are pinned as well.

Only the mean mode is tested here; the fluctuating buoyancy terms are not.
"""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

RE, PR, Q, L = 10.0, 1.0, 4.0 / 3.0, 2.0

DECK = """\
[mesh]
nx = 2
ny = 48
nz = 2
alfa0 = 0.5
beta0 = 1.0
stretching = 1.0
ymin = 0.0
ymax = {L}

[velocity]
ni = {RE}
meanpx = 0.0
meanpz = 0.0
meanflowx = {Q!r}
meanflowz = 0.0
u0 = 0.0
uN = 0.0
perturbation_amplitude = 1.0e-6
seed = 7

[scalars]
nPhi = 1
meantx = 0.0
meantb = 0.0
pr = {PR}
bc = neumann
sink = velocity
t0 = {t0}
tN = {tN}
ri = {ri}

[timestepping]
deltat = 0.04
cflmax = 0.0
time = 0.0
dt_field = 100000.0
dt_save = -1.0
t_max = 60.0
time_from_restart = false
nstep = 100000
"""


def read_list_directed(text, row_length):
    """Rows of a file of Fortran list-directed records (see check_scalar_isoflux.py)."""
    values = []
    for token in re.split(r"[,\s]+", text.strip()):
        count, star, value = token.rpartition("*")
        values.extend([float(value)] * (int(count) if star else 1))
    if len(values) % row_length:
        raise SystemExit(f"{len(values)} values, not a multiple of the {row_length} per row")
    return [values[i:i + row_length] for i in range(0, len(values), row_length)]


def reference(ri, t0, tN, n=48):
    """The laminar solution by Chebyshev collocation; INT phi dy = 0 (the run's
    own bulk is added back).  Unknowns U, phi (n+1 points each) and P.  The
    Neumann problem for phi is singular by one and consistent (gamma balances
    the wall fluxes), so one interior phi row makes way for the gauge row."""
    nu, kappa = 1.0 / RE, 1.0 / (RE * PR)
    gamma = kappa * (tN - t0) / Q
    m = n + 1
    x = np.cos(np.pi * np.arange(m) / n)            # y = L (1 - x)/2
    c = np.r_[2.0, np.ones(n - 1), 2.0] * (-1.0) ** np.arange(m)
    d = np.outer(c, 1.0 / c) / (x[:, None] - x[None, :] + np.eye(m))
    d1 = -2.0 / L * (d - np.diag(d.sum(axis=1)))
    d2 = d1 @ d1
    moments = np.zeros(m)                           # Clenshaw-Curtis weights
    moments[::2] = 2.0 / (1.0 - np.arange(0, m, 2) ** 2)
    w = L / 2.0 * np.linalg.solve(np.polynomial.chebyshev.chebvander(x, n).T, moments)
    eye, zero = np.eye(m), np.zeros((m, m))
    interior = list(range(1, m - 1))
    a = [np.c_[nu * d2, ri * eye, np.ones(m)][interior],
         np.c_[-gamma * eye, kappa * d2, np.zeros(m)][[i for i in interior if i != m // 2]],
         np.c_[eye, zero, np.zeros(m)][[0, m - 1]],
         np.c_[zero, d1, np.zeros(m)][[0, m - 1]],
         np.r_[w, np.zeros(m), 0.0][None, :], np.r_[np.zeros(m), w, 0.0][None, :]]
    rhs = np.zeros(2 * m + 1)
    rhs[-4:] = t0, tN, Q, 0.0
    sol = np.linalg.solve(np.vstack(a), rhs)
    u, phi = sol[:m], sol[m:2 * m]
    du0, duN = d1[0] @ u, d1[-1] @ u
    phi_m = w @ (u * phi) / Q
    ub = Q / L
    return dict(dudy0=du0, dudyN=-duN, cf0=2 * nu * du0 / ub**2, cfN=-2 * nu * duN / ub**2,
                phi0=phi[0], phiN=phi[-1], phim=phi_m,
                nu0=-t0 * 2 * L / (phi[0] - phi_m), nuN=tN * 2 * L / (phi[-1] - phi_m))


def run(channel, mpiexec, ri, t0, tN, env_extra):
    workdir = tempfile.mkdtemp(prefix="buoyancy_laminar_")
    try:
        with open(os.path.join(workdir, "dns.in"), "w") as fh:
            fh.write(DECK.format(L=L, RE=RE, Q=Q, PR=PR, t0=t0, tN=tN, ri=ri))
        env = dict(os.environ, CHANNEL_DISABLE_RESTART_WRITE="1", **env_extra)
        proc = subprocess.run([mpiexec, "-np", "1", channel], cwd=workdir, env=env,
                              capture_output=True, text=True, timeout=1200)
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout + proc.stderr)
            raise SystemExit(f"channel failed (exit {proc.returncode})")
        with open(os.path.join(workdir, "Runtimedata")) as fh:
            vel = read_list_directed(fh.read(), 13)
        with open(os.path.join(workdir, "Runtimedata.phi")) as fh:
            sca = read_list_directed(fh.read(), 12)
        return vel, sca
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def check(name, channel, mpiexec, ri, t0, tN, env_extra, failures):
    vel, sca = run(channel, mpiexec, ri, t0, tN, env_extra)
    v, s = vel[-1], sca[-1]
    # Runtimedata: time, dU/dy(0), -dU/dy(L), dW/dy x2, flow, px, flow z, pz, cfl, dt, Cf0, CfN
    # Runtimedata.phi: time, dphi/dy(0), dphi/dy(L), INT phi, corr, phi(0), phi(L),
    #                  var(0), var(L), Nu(0), Nu(L), phi_m
    ref = reference(ri, t0, tN)
    bulk = s[3]
    got = dict(dudy0=v[1], dudyN=v[2], cf0=v[11], cfN=v[12],
               phi0=s[5] - bulk / L, phiN=s[6] - bulk / L, phim=s[11] - bulk / L, nu0=s[9], nuN=s[10])
    worst = 0.0
    for key, want in ref.items():
        scale = max(abs(want), 1.0)
        err = abs(got[key] - want) / scale
        worst = max(worst, err)
        if not err <= 1.0e-4:
            failures.append(f"{name}: {key} = {got[key]:.8f}, laminar solution {want:.8f}")
    # Over the second half of the run, the change of INT phi dy against the heat
    # put in through the walls over the same time.  Analytically zero; the
    # discrete residue (quadrature vs. wall rows) is a few 1e-6 of the input,
    # while a sink with a wrong factor or sign would be O(1).
    mid = sca[len(sca) // 2]
    drift = abs(sca[-1][3] - mid[3])
    heat_in = (1.0 / RE / PR) * abs(tN - t0) * (sca[-1][0] - mid[0])
    if not drift <= 1.0e-4 * heat_in:
        failures.append(f"{name}: INT phi dy changed by {drift:.3e} while {heat_in:.3e} came in "
                        "through the walls; the velocity-weighted sink should balance it")
    drift /= heat_in
    print(f"{name}: worst relative deviation {worst:.2e}, Nu = {got['nu0']:.4f} / {got['nuN']:.4f}, "
          f"Cf = {got['cf0']:.5f} / {got['cfN']:.5f}, bulk drift {drift:.1e} of the heat input")


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_buoyancy_laminar.py <channel-binary> <mpiexec>")
    channel, mpiexec = os.path.abspath(sys.argv[1]), sys.argv[2]

    # The reference must reproduce the textbook laminar values at ri = 0.
    for t0, tN, want in ((-1.0, 1.0, 140.0 / 17.0), (-1.0, 0.0, 70.0 / 13.0)):
        got = reference(0.0, t0, tN)["nu0"]
        if abs(got - want) > 1.0e-8:
            raise SystemExit(f"reference solver: Nu = {got} at ri = 0, textbook {want}")

    failures = []
    check("both walls heated, ri = +0.5", channel, mpiexec, 0.5, -1.0, 1.0, {}, failures)
    check("one wall heated, ri = -0.2", channel, mpiexec, -0.2, -1.0, 0.0, {}, failures)
    if failures:
        for f in failures:
            print("FAIL:", f)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
