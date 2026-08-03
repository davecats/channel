# Prompt for the next session

Copy the block below as the opening message.

---

Continue simplifying `channel` for readability, without changing behaviour.

Read `.claude/projects/-home-ws-xt8786-Codes-channel/memory/MEMORY.md` first —
it points at the verification gates (including which machines have a usable
GPU), the NVHPC constraints on module structure, and the state of the layering
plan. Those were learned the hard way and will save you several wasted build
cycles.

## Where the previous session left it

The wall-normal solve is done. Every such problem — v, eta, scalars, d/dy,
d²/dy², Poisson p and dp/dy — is a pentadiagonal system plus four closure rows
(a ghost row and a wall row at each end), and closing it is always fold →
eliminate → solve → reconstruct. Each of those now has one definition:

- `wall_closure.fypph` (linsolve) — fold and eliminate, shared by the batched
  device kernels in `channel_equations` and the single-line host solver in
  `compact_line_solvers`.
- `channel_boundaries.fypph` (core) — the target-region scaffolding the four
  boundary-assembly callbacks used to each carry.
- `channel_stencils.fypph` — `at_lower_wall` / `at_upper_wall`, which own the
  span the one-sided wall rows must be applied to.
- `compact_assembly_context` now names its *outputs* too, so callbacks no
  longer write into public `y_line_solvers` workspace. `pressure_output` does
  not reference that module at all any more.

Two things were deliberately left duplicated, both recorded in comments at the
sites: `reconstruct` is spelled `a - sum(b)` on the host and
`a - b1 - b2 - b3` on the device, which do not round alike; and
`assemble_compact_derivative_interior` keeps the interleaved `ys_gpsv_ds`
view, which is genuinely a different layout from the owner-indexed one.

Four commits, each independently bit-identity verified on CPU and GPU.

`.fypp` sources are now `dnsdata`, `channel_equations`, `byte_workspace`,
`pressure_output`, `statistics`, `compact_line_solvers`.

## Candidates, roughly in value order

1. **`convvelo.f90` (1190 lines).** Never touched beyond mechanical changes,
   and now the largest thing that is not a solver. It mixes statistics
   accumulation, workspace management and MPI-IO in one module. Start by
   reading it and reporting what it actually does before changing anything.

2. **`y_line_solvers.f90` (2440 lines).** Still the biggest module: workspace
   management, the cuSPARSE/hipSPARSE wrapper, the endpoint-Schur solver, the
   pipelined-LU solver (MPI and NCCL variants written twice, ~460 lines) and
   the pentadiagonal kernels. The kernels **cannot** be extracted — that was
   tried and reverted, see the memory notes. The MPI/NCCL pipeline pair is
   genuinely different code (one non-blocking with overlap, one blocking), so
   unifying it would change behaviour; read both before deciding.

3. **The halo pack/unpack quartet in `y_line_solvers`** (`ys_pack_lower/upper_
   solution_halo`, `ys_unpack_lower/upper_solution_halo`). Four near-identical
   target regions differing only in a base index and a destination array. This
   is *index* dispatch, so a single parameterised routine is a better fit than
   fypp. It touches target regions, so gate on a GPU build.

4. **`FFT`/`IFT` and `RFT`/`HFT` in `ffts.f90`.** Each pair differs only by
   plan and direction across three `#ifdef` backends. Worth collapsing, but the
   backend switching makes it fiddly — judge whether the result is actually
   more readable before committing to it.

5. **`compute_cfl` and the transform orchestration still in `dnsdata`.** These
   are numerics rather than configuration; `dnsdata` would read better as pure
   setup if they moved.

## Two open questions that are the user's call, not yours

- **`FFTW_PATIENT`** (`ffts.f90`) makes CPU runs non-reproducible at ULP level
  run to run. Any bit-identity check needs a scratch build patched to
  `FFTW_ESTIMATE`. Switching the default is a real performance tradeoff —
  ask, don't decide.
- **`-Dbodyforce`** compiles but nothing builds or tests it. It hid two
  separate mistakes in an earlier session. Either give it a smoke test or
  delete it. The same goes for `-Dhalfchannel` and `-DphiNeumann`: a latent
  double-application of the wall fold survived in both for as long as it did
  precisely because nothing exercises them.

## How to work

Small, independently verified commits — one logical change each, so a
regression bisects to one peel. For every change: CPU `ctest` (30 tests), an
NVHPC GPU build, and bit-identity against the previous commit using a
`FFTW_ESTIMATE` scratch build. A green CPU suite has repeatedly been
insufficient: three separate module changes passed all 30 tests and broke the
GPU build.

Two techniques worth reusing:

- **Diff the generated Fortran, not just the output.** When a change is fypp
  expansion, expanding both trees and diffing the result proves the compiled
  program is unchanged — stronger than a bit-identical run, and it retires the
  performance question without needing an idle GPU. Shaping macros to take the
  *text of each call site's own expressions* rather than fixed variable names
  is what makes that achievable in hot kernels.
- **Compare against the checked-in references too**, not only base-vs-work. A
  bad test invocation once left stale files behind and reported a difference
  that did not exist.

Report honestly when something is reverted and why. Several extractions have
been abandoned across these sessions for good reasons, and recording that has
been more useful than forcing them through.
