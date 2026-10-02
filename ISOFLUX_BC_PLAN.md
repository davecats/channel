# Local constant wall heat flux (H2) for the passive scalars — implementation plan

Status: plan only, no code changed.
Target: run the existing Re_tau = 1000 / Pr = {0.025, 0.4, 1} scalar campaign a
second time with a **local, instantaneous** constant wall heat flux, so the wall
temperature fluctuates — the "isoflux" / "H2" condition of Kong, Choi & Lee
(2000) and Tiselj/Hetsroni (2004), rather than the isothermal ("H1") condition
used so far.

---

## 1. What the three papers actually prescribe

**Kong, Choi & Lee (2000), Phys. Fluids 12, 2555** — thermal boundary layer,
isothermal vs isoflux. Their Eq. (24) is the whole of it:

```
theta = 0           at the wall   (isothermal, H1)
d(theta)/dy = 1     at the wall   (isoflux,   H2)
```

That is a **Neumann condition on the full instantaneous temperature**, uniform
in x, z and t. Nothing is done to the mean: the wall temperature is free and
fluctuates. Everything else in that paper (the rescaling/recycling inflow
generator, the enthalpy-thickness Newton iteration, Eqs. 4–16) is about
generating turbulent inflow for a *spatially developing* boundary layer and is
irrelevant to us — our channel is streamwise periodic.

**Hetsroni, Tiselj, Bergant, Mosyak & Pogrebnyak (2004), J. Heat Transfer 126,
843** — turbulent flume, their Eqs. (5)–(6):

```
<theta+>_{x,z,t}(wall) = 0        (gauge: fixes the free additive constant)
d(theta)/dy(wall)      = Re_tau * Pr
```

Eq. (6) is again a uniform Neumann value on the total field. `Re_tau * Pr` is
simply what `d(theta+)/d(y/h) = 1` becomes in their non-dimensionalisation
(`theta+ = (T_w - T)/T_tau`, `T_tau = q_w/(u_tau rho c_p)`, `y` scaled by `h`).
Eq. (5) is not a second boundary condition — a pure-Neumann scalar has an
undetermined additive constant and they pin it.

**The trap, which the task statement already flags.** Hetsroni's Eq. (3) also
carries a source term `u_x+/u_B+`. That term is *not* the boundary condition.
It exists only because, in a streamwise-periodic box with a net wall heat input,
the mean temperature grows linearly in `x`; the term subtracts that growth so a
statistically stationary field exists. Kasagi et al. (1992) use the same device
with a **Dirichlet** wall (`theta = 0`), and that combination — "constant
*average* heat flux" — produces an isothermal wall with a fluctuating local flux,
i.e. the exact opposite of what we want. The two must not be conflated:

| | wall row | wall T | local wall flux | global balance closed by |
|---|---|---|---|---|
| H1 isothermal (current runs) | Dirichlet | fixed | fluctuates | volumetric source |
| Kasagi "const. average flux" | Dirichlet | fixed | fluctuates | `u/u_B` source |
| **H2 isoflux (wanted)** | **Neumann** | **fluctuates** | **exactly constant** | volumetric source |

The physics we are buying is in the last column of row 3 only. What closes the
global balance is a separate, secondary choice (section 5).

**Signature to validate against** (Kong §IV.A, Fig. 5b/7/9; Hetsroni Eqs. 9–10):

| quantity | H1 (isothermal) | H2 (isoflux) |
|---|---|---|
| `theta_rms` at the wall | 0 | finite (≈ 2.0 in wall units at Pr ≈ 0.7) |
| `theta'u'` as y → 0 | `~ y^2` | `~ y` |
| `theta'v'` as y → 0 | `~ y^3` | `~ y^2` |
| `alpha_t` as y → 0 | `~ y^3` | `~ y^2` |
| `Pr_t` as y → 0 | const ≈ 1.1 | `~ y` → 0 |
| thermal streak spacing | ≈ 100 wall units | ≈ 140 wall units |

---

## 2. What the code does today

The scalar wall rows live in exactly one place, `src/physics/channel_bcs.f90`:

```fortran
phi0bc = d040; phi0m1bc = der(1, 3, :)        ! Dirichlet, bottom
#ifdef phiNeumann
phi0bc = d140; phi0m1bc = der(1, 3, :)        ! Neumann,   bottom
#endif
phinbc = d04n; phinp1bc = der(ny - 1, 3, :)   ! Dirichlet, top
#ifdef phiNeumann
phinbc = d14n; phinp1bc = d04n                ! Neumann,   top   <-- broken
#endif
```

So a Neumann scalar BC is **already present** as the compile-time option
`-DCHANNEL_PHI_NEUMANN=ON` (`CMakeLists.txt:179`), and it is already *implicit*:
the row goes into the banded wall-normal system assembled in
`assemble_selected_compact_component_boundaries` and solved by
`solve_wall_normal_component`. That is the right architecture. Three things are
nevertheless missing or wrong.

The value side, `apply_wall_values`, sets `bc0(0,0,5+iPhi) = t0` on the mean mode
and `0` on every other mode — which, read through a Neumann row, is exactly
`d(phi)/dy = t0` uniformly in x and z. That is the condition we want, with no
further change.

**The current production configuration** (`newcases/re1000_Pr_0.025_1/dns.in`)
is `t0 = tN = 0`, `meantb = 2.0`, `nPhi = 3`, `pr = 0.025 0.4 1`, `ni = 20000`
(⇒ `Re_tau = 1000`, `u_tau = 0.05`). `meantb` drives the correction in
`linsolve_scalar` (`channel_equations.fypp:250`), which adds `corrtx(iPhi)`
times the solution of the Helmholtz problem with a **uniform unit source** and
homogeneous walls, chosen each substep so that the bulk stays at `meantb`. In
other words today's runs are *uniform volumetric heating with isothermal walls*,
with the source strength adjusted to hold the bulk temperature — itself a
"constant average heat flux" device, which is precisely why it must not be
mistaken for the isoflux condition.

---

## 3. The one real blocker: the top-wall ghost row is broken

`phinp1bc = d04n` under `phiNeumann` is a bug, and it is fatal rather than
inaccurate.

`d04n` is set in `case_setup.fypp:324` as `d04n = 0; d04n(1) = 1` — a *value*
row at node `ny`, with a **zero coefficient on the ghost node `ny+1`** (band
index `+2`). But every consumer of the upper ghost row divides by that entry:

- `src/linsolve/wall_closure.fypph:56` — `fold_ghost_into_upper_wall` divides by
  `eqnp1(2)`;
- `wall_closure.fypph:97` — `eliminate_upper_last_row` divides by `ghost_diag =
  ys_eqnp1_owner(2,...)`;
- `compact_component_solve.fypp:386` — the reconstruction of `phi(ny+1)` divides
  by `ys_eqnp1_owner(2,...)`;
- `compact_line_solvers.fypp:162` — the single-line host path, used by
  `solve_mean_correction_line`, does `x(ny+1) = (...)/upper_ghost_bc(2)`.

All four are a division by exactly zero. Building `-DCHANNEL_PHI_NEUMANN=ON`
today and running it produces NaN on the first substep, on CPU and GPU alike.
It has never been exercised: `cmake/check_physics_options.cmake` only checks that
the guarded branches *compile*, and `src/physics/README.md` says so in as many
words ("nothing checks that their numbers are right").

The bottom wall is correct: `phi0m1bc = der(1,3,:)` keeps the compact fourth-
derivative relation at `iy = 1` as the ghost closure, exactly as the Dirichlet
case does, and that relation is BC-independent. The fix is to make the top wall
mirror it — i.e. simply not to override `phinp1bc` at all:

```fortran
#ifdef phiNeumann
  phinbc = d14n                               ! leave phinp1bc = der(ny-1, 3, :)
#endif
```

That is the entire Stage 0 change: one deleted assignment.

### Is the system still solvable with Neumann at both walls? Yes.

Worth settling explicitly, because "pure Neumann is singular" is a true
statement about the wrong operator. What the scalars solve each substep is a
**Helmholtz**, not a Laplacian:

```
lambda*phi - alpha*(D^2 - k^2) phi = rhs ,   lambda = RK_rai(1,i)/deltat
```

and `RK_rai(1,i)` is `{3.75, 15, 6}`, so `lambda` is `O(1/deltat)` — never zero.
Multiplying the homogeneous problem by `phi` and integrating gives
`lambda*INT(phi^2) + alpha*INT(phi'^2) = 0` under homogeneous Neumann rows,
hence `phi = 0`: the operator is coercive and the solve is unique at **every**
mode, `k = 0` included. The familiar null space is the `lambda = 0` limit, which
the scalars never reach. (`apply_dy`/`apply_d2y` do pass `lambda = 0`, but they
use a different assembly and their own boundary callback, not these rows.)

Verified numerically on the production mesh by
`tools/check_scalar_wall_closure.py`, which rebuilds the stencils and assembles
the same system:

| | rank at k=0 | cond (row-equilibrated) |
|---|---|---|
| Dirichlet / Dirichlet (today) | 503/503 | 13.3 |
| Neumann / Neumann (fixed rows) | 503/503 | 15.9 |
| Neumann / Neumann at `lambda = 0` | 502/503 | 1.4e16 |

So Neumann is no worse conditioned than what is running today, and the
`lambda = 0` row is there only to show where the null space *would* be — its
smallest singular value is `1.7e-16` with exactly one null direction, the
constant. (Compare the raw condition numbers at your peril: they are `2e13` and
`3e10` and are dominated by the `1e12`-scale `D^4` ghost rows, which is a row
scaling, not conditioning.)

Three further things the same tool checks, because solvability of the operator
is not the same as solvability of what the code actually factorises:

- **All four closure pivots stay nonzero.** The fold/eliminate/reconstruct
  sequence divides by `eqm1(-2)`, the folded `eq0(-1)`, `eqnp1(+2)` and the
  folded `eqn(+1)`. Under Neumann the folded wall pivots are `-1.89e3` and
  `+1.89e3` where Dirichlet has `1.0` — different magnitude, nowhere near zero.
  Only the in-tree top ghost row gives exactly `0`, which is §3.
- **The prescribed gradients come back.** Solving with `t0 = +37`, `tN = -37`
  and re-applying `d140`/`d14n` to the result returns them to `1.4e-14`. The row
  is genuinely enforced, not merely present.
- **`corrtx` stays well defined.** `linsolve_scalar` divides by
  `yintegr(tcor)`, and under homogeneous *Neumann* rows `tcor` is
  `approx 1/lambda` across the channel rather than a bump vanishing at the
  walls, giving `INT(tcor) = 2.13e-3 approx 2/lambda` — nonzero, and in fact a
  fatter pivot than the Dirichlet case.

What Neumann at both walls does *not* determine is the long-time **mean level**,
since `d/dt INT(phi) dy = alpha*(phi'(2) - phi'(0))`. That is a stationarity and
gauge question, not a solvability one, and section 5 is where it is settled.

---

## 4. Why it must stay implicit (the "explicit will not run" warning, quantified)

The tempting shortcut is to keep the Dirichlet row and, once per step, set the
wall value by one-sided extrapolation of the previous field so the gradient comes
out right — i.e. do it inside `apply_wall_values` and touch nothing else. That is
an explicit Neumann closure and its stability limit is `alpha*dt/dy_1^2 < O(1)`.

On the production mesh (`ny = 500`, `stretching = 1.66`, `ymax = 2`) the first
spacing is `dy_1 = 9.674e-4`, i.e. `dy_1+ = 0.97`. With `alpha = ni/Pr`:

| Pr | `alpha*dt/dy_1^2` at dt = 1e-3 | at dt = 4e-3 | at dt = 1e-2 |
|---|---|---|---|
| 0.025 | 2.14 | 8.55 | 21.4 |
| 0.4 | 0.134 | 0.534 | 1.34 |
| 1 | 0.053 | 0.214 | 0.534 |

The Pr = 0.025 scalar is unconditionally unstable at any timestep the CFL
condition would ever allow, and Pr = 0.4 is marginal. The implicit route — the
Neumann row inside the banded system, which is what the code already does — has
no such limit.

Note also that the *existing* implicit path is already correct in the respect
that matters here: `solve_wall_normal_component` reconstructs `phi(0)`, `phi(-1)`,
`phi(ny)`, `phi(ny+1)` after every solve
(`compact_component_solve.fypp:356–395`), so the explicit (Crank–Nicolson /
RK) half built in `buildrhs_prepare` sees valid wall and ghost planes, and
`build_products` already loops `y_first-2 .. y_last+2`, so the now-nonzero wall
value participates in the advection products as it should. No change is needed
anywhere in `channel_equations.fypp`.

---

## 5. Choosing the configuration: what closes the global balance

With Neumann at both walls the mean level is a free constant and the integrated
balance is `d/dt ∫phi dy = alpha*(phi'(2) - phi'(0))`. Two well-posed choices:

**(A) Flux-through.** `t0 = tN = +g`: the same uniform flux enters one wall and
leaves the other. Globally balanced, no source needed, mean profile monotonic
across the channel. Both walls still see a genuine local isoflux condition.

**(B) Symmetric, matching the existing campaign — recommended.**
`t0 = +g`, `tN = -g`, with the existing `meantb` correction left on. The walls
remove heat uniformly and the volumetric source puts it back, exactly as in
today's isothermal runs; only the wall row changes. Crucially, the correction
field `tcor` is the solution of the Helmholtz problem with **homogeneous**
boundary rows — under Neumann those are zero-gradient rows, so the correction
contributes *no wall flux at all* and acts as a near-uniform offset. The imposed
gradient stays exactly `g`; the bulk pin is pure gauge-fixing. In steady state
the source self-adjusts to `S = alpha*g`.

Variant (B) is the right one because it is a pure H1 → H2 swap against the
database you already have: same Re, same Pr set, same box, same mesh, same
volumetric driving, same bulk value. Nothing moves but the wall.

### Sign convention

The rows are `d(phi)/dy` at `y = 0` and at `y = 2` in the global `+y` sense
(`d140` and `d14n` are both plain first-derivative rows). Heat leaves the fluid
at both walls — matching today's source-heated, wall-cooled runs — for
`t0 = +g`, `tN = -g`. Flip both signs for hot walls.

### Choosing g, and why it should be per scalar

In code units `alpha = ni/Pr`, `u_tau = Re_tau*ni`, so the friction temperature
of a scalar driven at gradient `g` is

```
phi_tau = alpha*g/u_tau = g/(Re_tau*Pr)
```

So `g = Re_tau*Pr` gives `phi_tau = 1` and makes `phi` literally Hetsroni's
`theta+` — their Eq. (6) is this identity. At Re_tau = 1000 that is
`g = 25, 400, 1000` for Pr = 0.025, 0.4, 1.

But there is a better criterion than tidiness: **pick `g_i` equal to the mean
wall gradient the matching isothermal run already has.** Then the isoflux run
starts in near-equilibrium instead of relaxing its entire mean profile, which
otherwise costs a very long transient (the global thermal relaxation time is
`h^2/alpha_t ~ O(200)` time units, and far longer at Pr = 0.025). Read `g_i`
directly from the existing `Runtimedata.phi`, whose columns are
(`runtime_diagnostics.fypp:60`):

```
time | lower-wall mean dphi/dy (nPhi cols) | upper-wall (nPhi cols) | bulk (nPhi) | corrtx (nPhi)
```

Take the time-averaged column `1+i` of the isothermal run and set
`t0_i = +g_i`, `tN_i = -g_i`.

Because `g_i` differs per scalar, and because `t0`/`tN` are today single deck
reals shared by all scalars, this is the main argument for Stage 2 below. With
Stage 0 alone you can only run all scalars at one shared `g` (e.g. `t0 = 1`,
`tN = -1`, rescaling by `phi_tau_i = 1/(Re_tau*Pr_i)` offline) and accept a long
mean-profile transient on two of the three.

---

## 6. Staged implementation

### Stage 0 — unblock the existing option (1 line)

`src/physics/channel_bcs.f90`: drop `phinp1bc = d04n` from the `phiNeumann`
branch so the top keeps `der(ny-1, 3, :)`, mirroring the bottom.

Build: `cmake -S . -B build-h2 -DCHANNEL_PHI_NEUMANN=ON && cmake --build build-h2 -j`.

Deliverable: an all-Neumann binary. Enough to validate the numerics and to run a
single-BC campaign.

### Stage 1 — validate before burning node-hours

1. **Unit test (would have caught the bug).** Extend
   `src/test/test_mean_correction.f90`, or add `test_scalar_neumann.f90` beside
   it: solve `lambda*phi - alpha*(D^2 - k^2)phi = f` on the real mesh at `k = 0`
   through `solve_full_line_compact` with `phi0bc = d140`,
   `phi0m1bc = der(1,3,:)`, `phinbc = d14n`, `phinp1bc = der(ny-1,3,:)` and
   inhomogeneous wall values `a`, `b`; compare against the closed-form
   cosh/sinh solution. Assert the recovered `phi'(0)`, `phi'(2)` equal `a`, `b`
   to round-off. Register it in `cmake/tests.cmake`.
2. **Conservation check, end to end.** `nPhi = 1`, Neumann, `t0 = tN = g`
   (variant A), `meantb = 0`, short run from a turbulent restart. Assert
   `d/dt ∫phi dy ≈ 0` and that the wall-gradient columns of `Runtimedata.phi`
   reproduce `g` to round-off — this is the cheapest direct proof that the row
   is being enforced on the instantaneous field, not just on the mean.
3. **Regression on the unchanged path.** `ctest` with the default build must be
   bit-identical: Stage 0 touches only a `#ifdef phiNeumann` branch.
4. **Physics acceptance.** Short Pr = 1 run; check the near-wall exponents of
   section 1 (`theta'u' ~ y`, `theta'v' ~ y^2`, `Pr_t ~ y`, finite wall
   `theta_rms`). These are local to the wall and converge quickly, long before
   the mean profile does.

### Stage 2 — per-scalar, runtime-selected BC (recommended)

The compile-time, all-or-nothing flag is the real limitation. Making the choice
per scalar lets **one run carry both conditions** — e.g. `nPhi = 6`,
`pr = 0.025 0.4 1 0.025 0.4 1`, three Dirichlet and three Neumann — over the
*same velocity realisation*. That is what makes Kong's Fig. 12 comparison
(instantaneous `theta'` fields for H1 and H2 at the same instant, same flow)
possible at all, and it removes any doubt that a difference comes from the BC
rather than from sampling.

Changes, all small and confined:

| file | change |
|---|---|
| `channel_state.f90` | `phi0bc`/`phi0m1bc`/`phinbc`/`phinp1bc` become allocatable rank-2 `(-2:2, nPhi)`; add `t0s(:)`, `tNs(:)`, `phi_bc_kind(:)`. Keep the module `use`-free (the nvfortran constraint). |
| `case_setup.fypp` | parse `[scalars] bc`, `t0`, `tn` as vectors via the existing `require_real_vector` / a new `require_string_vector`; default to `dirichlet` and broadcast a scalar `t0`/`tn` for backward compatibility. Map `t0s`/`tNs` to the device the way `pra` is (`target enter data map(to:)`). |
| `channel_bcs.f90` | `setup_boundary_conditions`: loop `iPhi`, fill Dirichlet or Neumann rows per scalar (with the Stage 0 fix). `apply_wall_values`: `bc0(0,0,5+iPhi) = t0s(iPhi)`, `bcn(...) = tNs(iPhi)`. |
| `channel_equations.fypp` | `linsolve_scalar` and its `solve_mean_correction_line` call index the rows with `iPhi`. |
| `CMakeLists.txt` | keep `CHANNEL_PHI_NEUMANN` as a deprecated alias that sets the default kind, or retire it. |

No kernel changes: the rows are copied host-side into the assembly context
(`assemble_selected_compact_component_boundaries`) before any `target` region,
so the generated device code is unchanged.

Deck after Stage 2:

```ini
[scalars]
nPhi = 6
pr   = 0.025 0.4 1 0.025 0.4 1
bc   = dirichlet dirichlet dirichlet neumann neumann neumann
t0   = 0 0 0   <g1>  <g2>  <g3>
tn   = 0 0 0  -<g1> -<g2> -<g3>
meantb = 2.0
meantx = 0.0
```

### Stage 3 — diagnostics

Under Neumann the wall-gradient columns in `Runtimedata.phi` become constants
(useful only as an assertion) and the interesting quantity — the **mean wall
value** `<phi_w>`, which gives `Nu` and the wall-to-bulk difference — is not
printed. Add `mean_line_scalar(0)` and `mean_line_scalar(ny)` as extra columns in
`runtime_diagnostics.fypp` (the full mean line is already gathered there, so this
is two array reads).

Good news on statistics: `convvelo` already accumulates over
`ny0-2 .. nyN+2`, i.e. **including the wall plane** `iy = 0`. So wall-temperature
rms/skewness/flatness and, in particular, Hetsroni's `U_cq+` — the convection
velocity of the wall temperature fluctuation, the paper's central result — come
out of the existing machinery with no change. The Python side (`derive.py`,
`postpro/`) is BC-agnostic; only the normalisation changes, and
`phi_tau_i = g_i/(Re_tau*Pr_i)` is the one constant to add.

### Stage 4 — restart files and the production run

- The restart header stores `nx, ny, nz, alfa0, beta0, ni, a, ymin, ymax, time`
  but **not `nPhi`** (`restart_io.f90:47`), and the MPI-IO subarray is built for
  `3 + nPhi` components. Going from 3 to 6 scalars will therefore read past the
  end of the existing 50 GB `Dati.cart.out` with no diagnostic. Write a small
  converter (alongside `postpro/interpolate_restart.py`) that appends copies of
  the three existing scalar planes as scalars 4–6 — which also gives the Neumann
  scalars a near-equilibrium start. Consider adding `nPhi` to the header check
  while you are there.
- Cost: `nPhi` 3 → 6 grows the field and the restart by ~1.5x and adds three more
  Helmholtz solves per substep. If that is unacceptable, run a **separate**
  3-scalar Neumann job from the same restart: the scalars are passive, so the
  velocity path is untouched and the two runs stay statistically identical —
  you only lose the instantaneous, same-realisation comparison.
- Initial condition: `initial_condition.f90:83` seeds scalars with
  `1.5*y*(2-y)` (bulk 2), which is a Dirichlet profile. Irrelevant when
  restarting, but if you ever cold-start an H2 case, give it the offset profile.
- Allow a long transient before sampling. The near-wall H2 signature equilibrates
  in a few eddy turnovers; the mean profile does not. Choosing `g_i` from the
  isothermal run (section 5) is what keeps this affordable, especially at
  Pr = 0.025.

### Stage 5 (optional) — the exact Kasagi/Tiselj source

Our global balance is closed by a *uniform* volumetric source (the existing
`meantb` mechanism), whereas Tiselj/Hetsroni and Kasagi use `-u_x/u_B`, which
arises from a mean temperature growing linearly in `x`. The wall physics — all of
section 1's signature — is set by the wall row and is unaffected; the difference
shows up in the outer-layer mean profile and weakly in `theta'u'`, because a
`u`-weighted source correlates with the streamwise velocity. Since the existing
isothermal database is already driven by the uniform source, keeping it is the
*self-consistent* choice and Stage 5 is only needed to compare one-for-one with
Tiselj's numbers. If wanted: add a term `-u/u_B * phi_x` to the scalar RHS in
`buildrhs` (physical space, alongside the existing advection products), and the
dormant `meantx` input — currently parsed and used only in a diagnostic print —
is the natural place to hang its coefficient.

---

## 7. Risk register

| risk | mitigation |
|---|---|
| Top ghost row divides by zero (§3) | Stage 0, plus the Stage 1 unit test that pins the recovered gradients |
| Explicit Neumann attempted as a shortcut | §4 — unstable at Pr = 0.025 for any usable dt; the implicit row is already in place |
| Confusing the bulk-pinning source with the BC | Under Neumann `tcor` has homogeneous **gradient** rows ⇒ contributes zero wall flux; the imposed `g` is exact. Assert it via `Runtimedata.phi` |
| Very long mean-profile transient | Choose `g_i` from the isothermal run's measured wall gradient (§5) |
| Silent restart corruption when `nPhi` changes | Converter + add `nPhi` to the header validation (§Stage 4) |
| `-DCHANNEL_PHI_NEUMANN=ON` also changes half-channel behaviour | Not used here; the Stage 2 runtime switch retires the coupling entirely |
| Default build regression | Stage 0/2 touch guarded branches and host-side row assembly only; `ctest` must stay bit-identical |

---

## 8. Shortest path to a first result

1. Delete `phinp1bc = d04n` from `channel_bcs.f90`.
2. Build with `-DCHANNEL_PHI_NEUMANN=ON`.
3. `nPhi = 1`, `pr = 1`, `t0 = <g>`, `tn = -<g>`, `meantb = 2.0`, short run from
   the existing restart, where `<g>` is the time-averaged lower-wall gradient of
   the Pr = 1 scalar in the current `Runtimedata.phi`.
4. Confirm: `Runtimedata.phi` wall-gradient columns are flat at `<g>`;
   `theta_rms` at `iy = 0` is finite and non-zero; `theta'v' ~ y^2` near the wall.

Then do Stage 2 before committing the full campaign — the per-scalar switch is
what buys the H1/H2 comparison on one velocity realisation, and it is a few
dozen lines.
