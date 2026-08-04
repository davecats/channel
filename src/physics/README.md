# `src/physics` — the layer you change

Everything in the rest of `src/` is machinery: transforms, transposes,
wall-normal solvers, I/O. This directory is the problem being solved. If you
are adapting `channel` to a different flow, these are the files, and for most
changes there is exactly one.

| I want to change… | Edit |
| --- | --- |
| Reynolds number, box size, mesh, scalars, flow rate, wall velocities | nothing here — `dns.in`, see the README's *Input File* section |
| what the wall condition **is** (values), including in time | `channel_bcs.f90` → `apply_wall_values` |
| which nodes a wall condition **constrains** (rows) | `channel_bcs.f90` → `setup_boundary_conditions` |
| the field a restartless run starts from | `initial_condition.f90` |
| the viscous operators | `channel_operators.fypph` |
| the nonlinear terms, the right-hand side, the mean-flow correction | `channel_equations.fypp` |
| half channel, or Neumann scalars at the wall | neither — CMake options, see below |

`channel_state.f90` is the run state itself (fields, physical parameters,
clock). You read it constantly; you rarely change it.

## The equations

Incompressible Navier–Stokes in velocity–vorticity form. Wall-parallel
derivatives are algebraic in each Fourier mode, so each mode is a pair of
wall-normal problems,

```
v   :  lambda*(D^2 - k^2) v  -  ni*(D^2 - k^2)^2 v  =  explicit terms
eta :  lambda*eta           -  ni*(D^2 - k^2) eta   =  explicit terms
```

plus one Helmholtz problem per passive scalar. Each operator is written **once**
in `channel_operators.fypph` and used by both halves of its equation — the
implicit half assembled into the banded matrix as a slice `j = -2:2`, the
explicit half applied at a single offset. That is the point of the file: these
expressions were once written twice, hundreds of lines apart and in different
notations, so a change to the viscous model could be applied to one half and
silently not the other. **Change the operator there and both halves follow.**

## Compile-time options

Two choices are made at build time, not in the deck:

```bash
cmake -S . -B build -DCHANNEL_HALF_CHANNEL=ON   # half channel
cmake -S . -B build -DCHANNEL_PHI_NEUMANN=ON    # Neumann scalars at the wall
```

They stay compile-time because each selects which stencil rows close the wall;
making them runtime would put a branch in kernels that execute for every mode
on every substep. The `physics_options_compile` test keeps the guarded branches
building, but **nothing checks that their numbers are right** — if you rely on
either, validate it against a case you know.

## Two rules that are not style

Both come from how this code reaches the GPU, and breaking either produces a
build failure or a dead kernel rather than a warning.

1. **Shared code that runs inside `!$omp target` must be a fypp macro, not a
   procedure.** A `declare target` procedure called across a module boundary
   does not survive `nvlink`. This is why `channel_operators.fypph` and
   `stencil_macros.fypph` are macro files rather than modules of functions.

2. **A module that declares `declare target` variables must not `use` anything.**
   `nvfortran` fails to resolve such variables inside downstream target regions
   when the declaring module also has `use` statements. `channel_state.f90`,
   `channel_grid.f90` and `stencil_coefficients.f90` are use-free for this
   reason. If you add state, add it to one of those, not to a module that
   computes something.

A corollary worth knowing when you change a hot path: a macro should take *the
text of each call site's own expressions* rather than fixed variable names. Then
neither call site changes its storage, the generated Fortran for the kernel is
unchanged, and you can prove it by diffing the expanded output instead of
having to benchmark.

## After you change something

The suite never runs the `channel` binary — every test drives its own `test_*`
program — so a green `ctest` is necessary and not sufficient. Run the binary and
check its exit status too:

```bash
ctest --test-dir build -j8
mpirun -np 2 build/channel        # from a directory holding dns.in
```

If the change is meant to preserve behaviour, the cheap proof is to compare
output against the parent commit with both trees built using `FFTW_ESTIMATE`
instead of `FFTW_PATIENT` (`src/fft/ffts.fypp`) — `FFTW_PATIENT` re-plans per
run, so CPU output is not otherwise reproducible at ULP level.

Note that a run with **no** restart file cannot be compared with anything:
`initial_condition.f90` calls `RANDOM_NUMBER` without ever seeding it from the
deck, so two runs of the same binary start from different fields.
