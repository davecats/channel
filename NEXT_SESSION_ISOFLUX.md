# Prompt for the next session: the isoflux (H2) production run

This is a **separate workstream** from `NEXT_SESSION.md`, which is still
outstanding and is about readability refactoring of the pipelined-LU send side.
Do not mix the two in one commit.

**The code is done and the production case is built.** Stages 0 to 4 of
`ISOFLUX_BC_PLAN.md` are implemented, tested and committed — see §9 of that file
for what landed and what the verification established. The four configuration
choices were settled on 2026-10-02 and `newcases_direct/re1000_Pr_0.025_1_isoflux`
is prepared from them. What is left is to submit the chain, and to run the GPU
suite.

Copy the block below as the opening message.

---

Continue the **local constant wall heat flux (isoflux, H2) scalar campaign**.
The boundary condition is implemented and verified; read `ISOFLUX_BC_PLAN.md`
§9 first for the state, and §1 and §5 for the physics and the configuration
question. `src/physics/README.md` and the "How to work" section of
`NEXT_SESSION.md` apply unchanged.

The condition is a deck setting, per scalar — the production case uses three
Neumann scalars, but one run can carry both conditions over the same velocity
field:

```ini
[scalars]
nPhi = 6
pr   = 0.025 0.4 1 0.025 0.4 1
bc   = dirichlet dirichlet dirichlet neumann neumann neumann
t0   = 0 0 0   5.361  30.717  51.701
tN   = 0 0 0  -5.361 -30.717 -51.701
meantb = 2.0
meantx = 0.0
```

`bc`, `t0` and `tn` take one value or one per scalar. Under `neumann`, `t0`/`tN`
are `dphi/dy` at the wall in the global `+y` sense, so cooled walls with a
volumetric source are `+g` and `-g`.

---

## Settled — do not re-derive, and do not re-verify

The physics is in `ISOFLUX_BC_PLAN.md` §1–§5 and has not changed. On top of
that, these are now measured facts rather than expectations:

- **The row is exact on the instantaneous field.** Every Neumann run reports the
  imposed gradient at both walls to within 1.1e-15 at every step.
  `scalar_wall_closure` pins the whole closure, ghost reconstruction included,
  against an exact quartic to 8.2e-15 at the production `lambda`.
- **The global balance closes** to 0.024% at ny = 64, and the residual converges
  at the scheme's own order in both h and dt — so there is no leak at the wall.
  `scalar_isoflux_balance` asserts this in 6 s as part of `ctest`.
- **`meantb` adds no wall flux.** With it on, the bulk holds and the imposed
  gradient stays exact. This was the open worry in the risk register; it is
  closed.
- **Both conditions in one run work**, at the same Prandtl number over the same
  velocity field.
- **`Runtimedata.phi` carries the settling monitors**: the mean wall value and
  the wall temperature variance, per scalar, per wall, every step. The variance
  is the one to read — it starts at exactly zero on a restart from an isothermal
  field and has to grow to a plateau. The gradient, bulk and `corrtx` columns are
  all pinned by the condition and say nothing about settling; §9 of the plan
  tabulates which signal does what, and why `corrtx` in particular is a trap.
- **`CHANNEL_PHI_NEUMANN` is gone.** Configuring with it is a FATAL_ERROR that
  names the deck keys. Do not look for it.

## The four decisions — settled 2026-10-02, and the case is built

1. **Symmetric.** `t0 = +g`, `tN = -g`, `meantb = 2.0` left on. A pure H1 → H2
   swap against the existing database.
2. **Three scalars, a separate run.** No restart in the workspace holds six
   scalars at this Re and Pr set — the only many-scalar cases are
   `re180_Pr_0.025_1` and `re500_Pr_0.025_1`, eight scalars at
   `Pr = 0.005 … 1`, a different Pr set and a different Re. So the
   same-realisation comparison was not available, and the run keeps the
   campaign's own Pr/Re. Nothing to convert: `nPhi` stays 3.
3. **`g_i` measured**: 5.361, 30.717, 51.701 at Pr = 0.025, 0.4, 1, from the
   isothermal run's own `Runtimedata.phi` over t = 3722 .. 3825, symmetrised
   between the walls (the last 20% instead moves them by under 0.4%). The tidy
   `g = Re_tau*Pr = 25, 400, 1000` is 4.7x, 13x and 19.3x those, so the mean
   profile would have to grow by that factor before anything could be sampled.
   Rescale offline by `phi_tau_i = 0.214, 0.0768, 0.0517` for Hetsroni's
   `theta+`.
4. **5 eddy turnovers of warm-up, then the isothermal run's own statistics
   length.** `T_eddy = Re_b/Re_tau = 20` time units; the restart is at
   t = 3825.152, so `[convvelo] t_start = 3925.2` and the warm-up is simply the
   part of the run convvelo does not see. The isothermal run covered 1646.0 time
   units = 83 T_eddy in 8 links of 90000 steps (205.7 each), so **9 links** give
   1752 sampled = 85 T_eddy.

The case is **built and not submitted**:

    /hkfs/work/workspace/scratch/.../newcases_direct/re1000_Pr_0.025_1_isoflux/

Its `README.md` carries the reasoning, the column map of `Runtimedata.phi`, what
to watch and what the signals mean. Four deck lines differ from
`re1000_Pr_0.025_1` and nothing else does. `CASE_BIN=../bin_isoflux/channel`
pins a binary built from this commit — the campaign's `bin/channel` would *not*
reject this deck, it would ignore `bc` and read only the first number of `t0`
and run a plausible-looking isothermal case, so `run_case.sh` now fails a link
whose deck asks for `neumann` and whose banner does not report one. To launch:

    cd .../newcases_direct && ./chain.sh re1000_Pr_0.025_1_isoflux 9

about 10650 GPU-h and 3.7 TB.

## Still outstanding

- **Nothing has run on a GPU yet.** The NVHPC build compiles all 20 targets with
  `-mp=gpu`; `apply_wall_values` generates its two kernels and
  `wall_plane_variance` its one, so the per-scalar `t0s(iPhi)` and the wall
  reduction inside target regions are fine at compile time. But this box is a
  login node with no GPU, so gate leg 2 is open. A production-scale smoke
  (`./smoke.sh re1000_Pr_0.025_1_isoflux 60`, 8 nodes, 60 steps, restart read
  through a symlink with `CHANNEL_DISABLE_RESTART_WRITE=1`) is the first thing to
  look at: the gradient columns flat at 5.361 / 30.717 / 51.701 and the variance
  columns growing from zero is the whole check. Establish the failing set at HEAD
  on a GPU node before reading anything into a `ctest` result, as
  `HPC_SESSION.md` says.
- **Stage 5 is untouched and still optional**: the exact `-u_x/u_B`
  Kasagi/Tiselj source. The wall physics does not depend on it. Only needed to
  compare one-for-one with Tiselj's numbers rather than self-consistently with
  your own isothermal database.
- **`NEXT_SESSION.md` is a different, still-open subject** — the pipelined-LU
  send side.
