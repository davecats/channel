# Prompt for the next session: local constant wall heat flux (H2) scalars

This is a **separate workstream** from `NEXT_SESSION.md`, which is still
outstanding and is about readability refactoring of the pipelined-LU send side.
Do not mix the two in one commit.

Copy the block below as the opening message.

---

Continue the work on the **local constant wall heat flux boundary condition**
for the passive scalars. The physics, the reading of the three source papers,
the audit of the current code and the staged implementation plan are already
done and live in `ISOFLUX_BC_PLAN.md` at the repo root — read that file first,
in full, before touching anything. This file is only the handoff: what is
settled, what the first move is, and what not to re-litigate.

The goal is a second pass over the existing Re_tau = 1000 / Pr = {0.025, 0.4, 1}
scalar campaign in which the wall imposes a uniform *instantaneous* heat flux
and the wall temperature is free to fluctuate — Kong, Choi & Lee (2000) Eq. (24)
and Hetsroni et al. (2004) Eq. (6) — instead of the isothermal wall used so far.

Also read `src/physics/README.md` and the "How to work" section of
`NEXT_SESSION.md`; the verification gate there applies unchanged, and leg 3 of
it is the one that matters most here (see below).

---

## Settled, do not re-derive

These came out of reading the three papers against the code. They are
established; the plan file carries the evidence.

- **The boundary condition is a Neumann row on the total instantaneous field**,
  uniform in x, z and t. Hetsroni's `<theta+>(wall) = 0` is not a second
  boundary condition, it pins the free additive constant of a pure-Neumann
  scalar.
- **Kasagi's "constant average heat flux" is not this.** It is a Dirichlet wall
  plus a source term, and gives a fluctuating *flux* at a fixed-temperature
  wall — the opposite. Your current production deck
  (`newcases/re1000_Pr_0.025_1/dns.in`: `t0 = tN = 0`, `meantb = 2.0`) is
  already a variant of that device, via the `corrtx` uniform-source correction
  in `linsolve_scalar`. Do not mistake it for the thing being built.
- **Everything in Kong beyond Eq. (24)** — the rescaling/recycling inflow
  generator, the enthalpy-thickness Newton iteration, Eqs. (4)–(16) — is inflow
  generation for a spatially developing boundary layer and does not apply to a
  periodic channel.
- **Explicit is not an option.** On the production mesh `dy_1+ = 0.97`, so an
  explicit Neumann closure gives `alpha*dt/dy_1^2 = 2.1 … 21` at Pr = 0.025.
  Unconditionally unstable there, marginal at Pr = 0.4. The row must stay in the
  implicit banded system, which is where the code already puts it.
- **No change is needed in `channel_equations.fypp`.**
  `solve_wall_normal_component` already reconstructs `phi(0)`, `phi(-1)`,
  `phi(ny)`, `phi(ny+1)` after every solve, and `build_products` already loops
  `y_first-2 .. y_last+2`, so the now-nonzero wall value enters the advection
  products correctly on its own.
- **`convvelo` already accumulates the wall plane** (`ny0-2 .. nyN+2`), so
  wall-temperature statistics and Hetsroni's `U_cq+` come out of the existing
  machinery with no change. The Python post-processing is BC-agnostic; only the
  normalisation `phi_tau_i = g_i/(Re_tau*Pr_i)` is new.

## Do this first: the top-wall ghost row

`src/physics/channel_bcs.f90`, inside the `#ifdef phiNeumann` block for the top
wall:

```fortran
phinbc = d14n; phinp1bc = d04n      ! <-- the second assignment is the bug
```

`d04n` is `d04n = 0; d04n(1) = 1` (`case_setup.fypp:324`) — a value row at node
`ny` with a **zero coefficient on the ghost node `ny+1`**. Four places divide by
exactly that entry:

- `src/linsolve/wall_closure.fypph:56` — `fold_ghost_into_upper_wall`
- `src/linsolve/wall_closure.fypph:97` — `eliminate_upper_last_row`
- `src/numerics/compact_component_solve.fypp:386` — reconstruction of `phi(ny+1)`
- `src/linsolve/compact_line_solvers.fypp:162` — the single-line host path used
  by `solve_mean_correction_line`

So `-DCHANNEL_PHI_NEUMANN=ON` NaNs on the first substep, on CPU and GPU alike.
It has never been run: `cmake/check_physics_options.cmake` only checks that the
guarded branches compile, and `NEXT_SESSION.md`'s gate leg 3 says in as many
words that nothing builds `-DphiNeumann`.

The bottom wall is already right — `phi0m1bc = der(1,3,:)`, the compact
fourth-derivative relation at `iy = 1`, which is BC-independent. The fix is to
make the top mirror it, i.e. simply delete the override so `phinp1bc` keeps
`der(ny-1, 3, :)`:

```fortran
#ifdef phiNeumann
  phinbc = d14n
#endif
```

One deleted assignment. Then build `-DCHANNEL_PHI_NEUMANN=ON` and run the
binary — gate leg 3 exists precisely because a green `ctest` would not have
caught this.

## Then, in order

Full detail in `ISOFLUX_BC_PLAN.md` §6. In brief:

1. **Stage 1 — validate before burning node-hours.** Add the unit test that
   would have caught the bug: solve the 1-D Helmholtz through
   `solve_full_line_compact` with the Neumann rows and inhomogeneous wall values,
   compare against the closed-form cosh/sinh solution, and assert the recovered
   `phi'(0)`, `phi'(2)`. Register it in `cmake/tests.cmake`. Then a short
   `nPhi = 1` conservation run (`t0 = tN = g`, `meantb = 0`) asserting
   `d/dt ∫phi dy ≈ 0` and flat wall-gradient columns in `Runtimedata.phi`.
2. **Stage 2 — per-scalar, runtime-selected BC.** This is the one worth doing
   before the campaign, not after. It lets a single run carry three Dirichlet
   and three Neumann scalars over the *same velocity realisation*, which is what
   makes Kong's Fig. 12 comparison possible and removes any doubt that a
   difference comes from the BC rather than from sampling. Host-side row
   assembly only, no kernel changes, a few dozen lines. The file-by-file table
   is in the plan.
3. **Stage 3 — diagnostics.** Under Neumann the wall-gradient columns of
   `Runtimedata.phi` become constants (useful only as an assertion); add the
   mean wall *value* `<phi_w>` instead. Two array reads in
   `runtime_diagnostics.fypp`, the mean line is already gathered there.
4. **Stage 4 — restart conversion and the production run.** See the trap below.
5. **Stage 5 (optional)** — the exact `-u_x/u_B` Kasagi/Tiselj source, only if
   you want to compare one-for-one with Tiselj's numbers rather than
   self-consistently with your own isothermal database.

## Decisions that are the user's call, not yours

- **Which global-balance variant.** The plan recommends symmetric
  (`t0 = +g`, `tN = -g`, `meantb` left on) because it is a pure H1 → H2 swap
  against the existing database — same Re, Pr set, box, mesh, driving and bulk.
  Flux-through (`t0 = tN = +g`, no source) is the alternative. Do not switch
  without asking.
- **One 6-scalar run or two 3-scalar runs.** Six scalars costs ~1.5x in field
  memory, restart size (the current `Dati.cart.out` is 50 GB) and scalar solve
  time, and buys the same-realisation instantaneous comparison. Two separate
  runs from the same restart are statistically equivalent — the scalars are
  passive, so the velocity path is untouched — and cost nothing extra.
- **The value of `g_i`.** `g = Re_tau*Pr` makes `phi_tau = 1` and `phi` literally
  Hetsroni's `theta+`. But the better choice is `g_i` = the mean wall gradient
  the matching isothermal run already shows, so the isoflux run starts in
  near-equilibrium instead of relaxing its whole mean profile — which otherwise
  costs a very long transient, badly so at Pr = 0.025. Read it from the existing
  `Runtimedata.phi`, columns `2 .. 1+nPhi` (lower wall), averaged in time.

## Traps

- **The restart header does not record `nPhi`** (`restart_io.f90:47`), and the
  MPI-IO subarray is built for `3 + nPhi` components. Going from 3 to 6 scalars
  reads past the end of the existing 50 GB file **with no diagnostic**. Write
  the converter (append copies of the three existing scalar planes as scalars
  4–6, which also gives the Neumann scalars a near-equilibrium start) before
  changing `nPhi`, and consider adding `nPhi` to the header validation.
- **`t0`/`tN` are single deck reals shared by all scalars** today. Until Stage 2
  you can only drive every scalar at one shared `g`.
- **The bulk-pinning correction is not a wall flux.** `tcor` solves the
  Helmholtz problem with *homogeneous* boundary rows, which under Neumann means
  zero gradient — so it contributes no wall flux and the imposed `g` stays
  exact. Assert this from `Runtimedata.phi` rather than assuming it.
- **`-DCHANNEL_PHI_NEUMANN=ON` is all-or-nothing and also couples to the
  half-channel branch.** Stage 2 retires that coupling; until then do not build
  the two options together for a production run.
- **Default build must stay bit-identical.** Stage 0 and Stage 2 touch guarded
  branches and host-side row assembly only. Gate leg 4 (`FFTW_ESTIMATE`
  bit-identity against the parent commit) applies.

## Acceptance criteria

The H2 signature, from Kong §IV.A (Figs. 5b, 7, 9) and Hetsroni Eqs. (9)–(10),
against the H1 scalars of the same campaign:

| quantity | H1 | H2 |
|---|---|---|
| `theta_rms` at the wall | 0 | finite (≈ 2.0 in wall units at Pr ≈ 0.7) |
| `theta'u'` as y → 0 | `~ y^2` | `~ y` |
| `theta'v'` as y → 0 | `~ y^3` | `~ y^2` |
| `alpha_t` as y → 0 | `~ y^3` | `~ y^2` |
| `Pr_t` as y → 0 | const ≈ 1.1 | `~ y` → 0 |
| thermal streak spacing | ≈ 100 wall units | ≈ 140 wall units |

These are local to the wall and converge in a few eddy turnovers, long before
the mean profile does — so they are the right early check on a short run.
