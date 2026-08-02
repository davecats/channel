# Prompt for the next session

Copy the block below as the opening message.

---

Continue simplifying `channel` for readability, without changing behaviour.

Read `.claude/projects/-home-ws-xt8786-Codes-channel/memory/MEMORY.md` first — it
points at the verification gates, the NVHPC constraints on module structure,
and the state of the layering plan. Those were learned the hard way and will
save you several wasted build cycles.

## Where the previous session left it

`dnsdata` went from 1827 to ~530 lines. The physics now lives in
`channel_equations.fypp`, with the differential operators defined once in
`channel_operators.fypph` and used by both the implicit and explicit halves.
Supporting modules: `channel_grid` (mesh + decomposition), `channel_state`
(fields, BC planes, parameters, clock), `compact_stencils`, `channel_bcs`,
`statistics`, `restart_io`, `env_options`. fypp is wired into CMake.

Everything was verified bit-identical to the pre-refactor tree.

## Candidates, roughly in value order

1. **`convvelo.f90` (1192 lines) and `pressure_output.f90` (831).** Neither
   was touched beyond mechanical changes. `convvelo` in particular mixes
   statistics accumulation, workspace management and MPI-IO in one module, and
   is the largest remaining file after `channel_equations`. Start by reading it
   and reporting what it actually does before changing anything.

2. **`y_line_solvers.f90` (~2400 lines).** Still the biggest module. It holds
   workspace management, the cuSPARSE/hipSPARSE wrapper, the endpoint-Schur
   solver, the pipelined-LU solver (MPI and NCCL variants written twice, ~460
   lines) and the pentadiagonal kernels. The kernels **cannot** be extracted —
   that was tried and reverted, see the memory notes. The MPI/NCCL pipeline
   pair is genuinely different code (one non-blocking with overlap, one
   blocking), so unifying it would change behaviour; read both before deciding.

3. **The halo pack/unpack quartet in `y_line_solvers`** (`ys_pack_lower/upper_
   solution_halo`, `ys_unpack_lower/upper_solution_halo`). Four near-identical
   target regions differing only in a base index and a destination array. This
   is *index* dispatch, so a single parameterised routine is a better fit than
   fypp. Left undone deliberately; it touches target regions, so gate on a GPU
   build.

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
  separate mistakes last session. Either give it a smoke test or delete it.

## How to work

Small, independently verified commits — one logical change each, so a
regression bisects to one peel. For every change: CPU `ctest` (30 tests), an
NVHPC GPU build, and bit-identity against the previous commit using a
`FFTW_ESTIMATE` scratch build. A green CPU suite has repeatedly been
insufficient: three separate module changes passed all 30 tests and broke the
GPU build.

Report honestly when something is reverted and why. Two extractions were
abandoned last session for good reasons, and recording that was more useful
than forcing them through.
