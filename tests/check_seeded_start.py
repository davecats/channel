#!/usr/bin/env python3
"""Checks that a seeded restartless start field is reproducible.

Three properties, each of which has a way of being quietly false:

  1. the same seed gives the same field twice          (reproducible at all)
  2. the same seed gives the same field at np=1 and np=2, and at npy=2
                                                       (decomposition-independent)
  3. a different seed gives a different field          (the seed is really used)

Property 2 is the one that needs a counter-based generator rather than a seeded
intrinsic one, and is the reason the test drives real runs instead of unit
testing the mixer.

Each run uses nstep = 0, so the time loop never executes and the restart file
the run writes is exactly the generated start field.  Runs happen in a private
temporary directory because channel writes fixed output names.
"""

import os
import shutil
import subprocess
import sys
import tempfile

DECK = """\
[mesh]
nx = 7
ny = 12
nz = 4
alfa0 = 1.0
beta0 = 2.0
stretching = 1.6
ymin = 0.0
ymax = 2.0

[velocity]
ni = 1000.0
meanpx = 0.0
meanpz = 0.0
meanflowx = 2.0
meanflowz = 0.0
u0 = 0.0
uN = 0.0
perturbation_amplitude = 1.0e-3
{seed_line}

[scalars]
nPhi = 0
meantx = 0.0
meantb = 0.0
t0 = 0.0
tN = 0.0

[timestepping]
deltat = 0.01
cflmax = 0.0
time = 0.0
dt_field = 1000.0
dt_save = -1.0
t_max = 1000.0
time_from_restart = false
nstep = 0
"""


def run(channel, mpiexec, nprocs, seed, env_extra=None):
    """Run channel with no restart file present; return the field it writes."""
    workdir = tempfile.mkdtemp(prefix="seeded_start_")
    try:
        seed_line = "" if seed is None else f"seed = {seed}"
        with open(os.path.join(workdir, "dns.in"), "w") as fh:
            fh.write(DECK.format(seed_line=seed_line))

        env = dict(os.environ)
        if env_extra:
            env.update(env_extra)

        proc = subprocess.run(
            [mpiexec, "-np", str(nprocs), channel],
            cwd=workdir, env=env, capture_output=True, text=True, timeout=600,
        )
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout + proc.stderr)
            raise SystemExit(f"channel failed at np={nprocs} (exit {proc.returncode})")

        out = os.path.join(workdir, "Dati.cart.out")
        if not os.path.exists(out):
            sys.stderr.write(proc.stdout + proc.stderr)
            raise SystemExit(f"no Dati.cart.out written at np={nprocs}")
        with open(out, "rb") as fh:
            return fh.read()
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_seeded_start.py <channel-binary> <mpiexec>")
    # Each run happens in its own temporary directory, so a relative path to
    # the binary would not resolve there.
    channel, mpiexec = os.path.abspath(sys.argv[1]), sys.argv[2]

    failures = []

    a = run(channel, mpiexec, 1, 4242)
    b = run(channel, mpiexec, 1, 4242)
    if a != b:
        failures.append("same seed, two runs at np=1: fields differ")

    c = run(channel, mpiexec, 2, 4242)
    if a != c:
        failures.append("same seed at np=1 and np=2: fields differ "
                        "(start field depends on the decomposition)")

    d = run(channel, mpiexec, 2, 4242, env_extra={"CHANNEL_NPY": "2"})
    if a != d:
        failures.append("same seed at np=1 and npy=2: fields differ "
                        "(start field depends on the wall-normal split)")

    e = run(channel, mpiexec, 1, 4243)
    if a == e:
        failures.append("seeds 4242 and 4243 give the same field "
                        "(the seed is not reaching the perturbation)")

    if failures:
        for f in failures:
            print("FAIL:", f)
        raise SystemExit(1)
    print("seeded start field: reproducible, decomposition-independent, "
          "and sensitive to the seed")


if __name__ == "__main__":
    main()
