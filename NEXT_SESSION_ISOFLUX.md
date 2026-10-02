# Prompt for the next session: the isoflux (H2) production run

This is a **separate workstream** from `NEXT_SESSION.md`, which is still
outstanding and is about readability refactoring of the pipelined-LU send side.
Do not mix the two in one commit.

**The code is done.** Stages 0 to 4 of `ISOFLUX_BC_PLAN.md` are implemented,
tested and committed — see §9 of that file for what landed and what the
verification established. What is left is the run, and the four choices below
that are yours rather than the code's.

Copy the block below as the opening message.

---

Continue the **local constant wall heat flux (isoflux, H2) scalar campaign**.
The boundary condition is implemented and verified; read `ISOFLUX_BC_PLAN.md`
§9 first for the state, and §1 and §5 for the physics and the configuration
question. `src/physics/README.md` and the "How to work" section of
`NEXT_SESSION.md` apply unchanged.

The condition is now a deck setting, per scalar:

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
- **`CHANNEL_PHI_NEUMANN` is gone.** Configuring with it is a FATAL_ERROR that
  names the deck keys. Do not look for it.

## The four decisions, which are yours

1. **Which global-balance variant.** The plan recommends symmetric
   (`t0 = +g`, `tN = -g`, `meantb = 2.0` left on), because it is a pure H1 → H2
   swap against the existing database — same Re, Pr set, box, mesh, driving and
   bulk. Flux-through (`t0 = tN = +g`, no source) is the alternative, and is
   also verified to work. Nothing in the code prefers either.

2. **One 6-scalar run, or a second 3-scalar run from the same restart.** Six
   scalars costs ~1.5x in field memory, restart size (the current
   `Dati.cart.out` is 50 GB) and scalar solve time, and buys the
   same-realisation instantaneous comparison that Kong's Fig. 12 needs. Two runs
   are statistically equivalent and cost nothing extra, but are two samples.

3. **The `g_i`.** Use the measured equilibrium gradients of the matching
   isothermal run, so the isoflux scalars start near equilibrium instead of
   relaxing a whole mean profile. From
   `newcases_direct/re1000_Pr_0.025_1/Runtimedata.phi`, time-averaged over
   t = 3722 .. 3825 and symmetrised between the walls:

   | Pr | g_i | `Re_tau*Pr` | `phi_tau = g/(Re_tau*Pr)` |
   |---|---|---|---|
   | 0.025 | 5.361 | 25 | 0.214 |
   | 0.4 | 30.717 | 400 | 0.0768 |
   | 1 | 51.701 | 1000 | 0.0517 |

   The tidy alternative `g = Re_tau*Pr` would make `phi` literally Hetsroni's
   `theta+`, but it is 4.7x, 13x and 19.3x the equilibrium gradient, so the mean
   profile would have to grow by that factor before anything could be sampled —
   badly so at Pr = 0.025. Rescale offline by `phi_tau_i` instead if that
   normalisation is wanted.

4. **How long, and when to start sampling.** The near-wall H2 signature
   equilibrates in a few eddy turnovers; the mean profile does not. Choosing
   `g_i` as above is what keeps that affordable.

## Doing the run

1. **If going to 6 scalars, convert the restart first.** The reader now refuses
   a file with the wrong component count rather than reading past its end, so
   this is a hard stop rather than a silent corruption — but convert, do not
   fight it:

   ```bash
   postpro/add_restart_scalars.py Dati.cart.out Dati.cart.6phi.out --scalars 1 2 3 1 2 3
   ```

   The new scalars 4–6 start as copies of 1–3, which is the near-equilibrium
   start that makes decision 3 worth making. ~100 GB out, streamed.
2. **Short run first.** `nPhi = 1`, `pr = 1`, `bc = neumann`, `t0 = 51.701`,
   `tn = -51.701` on the production mesh from the existing restart, a few
   hundred steps. Confirm from `Runtimedata.phi` that the gradient columns are
   flat at `g` and that the new wall-value columns move. Then commit the
   campaign.
3. **Check the acceptance criteria early.** They are local to the wall and
   converge long before the mean profile:

   | quantity | H1 | H2 |
   |---|---|---|
   | `theta_rms` at the wall | 0 | finite (≈ 2.0 in wall units at Pr ≈ 0.7) |
   | `theta'u'` as y → 0 | `~ y^2` | `~ y` |
   | `theta'v'` as y → 0 | `~ y^3` | `~ y^2` |
   | `alpha_t` as y → 0 | `~ y^3` | `~ y^2` |
   | `Pr_t` as y → 0 | const ≈ 1.1 | `~ y` → 0 |
   | thermal streak spacing | ≈ 100 wall units | ≈ 140 wall units |

   `convvelo` already accumulates the wall plane (`ny0-2 .. nyN+2`), so wall
   temperature statistics and Hetsroni's `U_cq+` come out of the existing
   machinery. The Python post-processing is BC-agnostic; the one new constant is
   `phi_tau_i = g_i/(Re_tau*Pr_i)`.

## Still outstanding

- **Gate leg 2 has not been run.** The NVHPC build compiles all 20 targets with
  `-mp=gpu` and `apply_wall_values` still generates its two GPU kernels, so the
  per-scalar `t0s(iPhi)` inside the target region is fine — but this box is a
  login node with no GPU, so nothing has been *run* on one. Establish the
  failing set at HEAD on a GPU node before reading anything into a result, as
  `HPC_SESSION.md` says.
- **Stage 5 is untouched and still optional**: the exact `-u_x/u_B`
  Kasagi/Tiselj source. The wall physics does not depend on it. Only needed to
  compare one-for-one with Tiselj's numbers rather than self-consistently with
  your own isothermal database.
- **`NEXT_SESSION.md` is a different, still-open subject** — the pipelined-LU
  send side.
