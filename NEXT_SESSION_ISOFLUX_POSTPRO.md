# Prompt for the next session: post-processing the isoflux (H2) campaign

The simulations are **finished**. Six cases, 29 scalars, 33 links, 23 717 time
units sampled, completed 2026-10-24. `ISOFLUX_BC_PLAN.md` §9 is the
implementation record and `NEXT_SESSION_ISOFLUX.md` the run handoff; this file
is the handoff for turning the output into results.

The data lives on **HoreKa**, 13.7 TB of it, so this is a session that works
over `ssh` on HoreKa or pulls a subset — not one that copies the campaign to a
workstation. What is small enough to pull is listed under *Working off HoreKa*.

Copy the block below as the opening message.

---

Process the **isoflux (H2) scalar campaign** and compare it with its isothermal
(H1) twin. Both campaigns are complete and verified; read
`NEXT_SESSION_ISOFLUX_POSTPRO.md` in full first, then
`postpro/CONVVELO_DIRECT.md` for what the convection-velocity estimator stores
and why. The physics is in `ISOFLUX_BC_PLAN.md` §1.

The deliverable is the H1 → H2 comparison the campaign was run for: Hetsroni's
convection velocity of the wall temperature fluctuation, and the near-wall
signature of the isoflux condition, at Re_tau = 180, 500 and 1000 over Prandtl
0.005 to 30.

**Do not re-derive the physics or re-verify the runs.** Both are settled and
recorded. What is open is in *What to produce* and *Caveats*, and the caveats
are not optional reading — three of them change numbers.

---

## Where the data is

Everything is on HoreKa under

    /hkfs/work/workspace/scratch/xt8786-xt8786-scalars_highre-1787706663/newcases_direct

with each isoflux case beside its isothermal twin, same name plus `_isoflux`:

| case | H1 twin | H2 | nPhi | Pr |
|---|---|---|---|---|
| re180_Pr_0.025_1 | ✓ | `_isoflux` | 8 | 0.005 0.025 0.05 0.1 0.2 0.4 0.7 1 |
| re180_Pr_1_7 | ✓ | `_isoflux` | 4 | 1 2 4 7 |
| re180_Pr_7_30 | ✓ | `_isoflux` | 3 | 7 15 30 |
| re500_Pr_0.025_1 | ✓ | `_isoflux` | 8 | 0.005 … 1 |
| re500_Pr_1_2_4 | ✓ | `_isoflux` | 3 | 1 2 4 |
| re1000_Pr_0.025_1 | ✓ | `_isoflux` | 3 | 0.025 0.4 1 |

Sizes: 0.2, 1.5, 3.3, 2.3, 2.9 and 3.5 TB for the six isoflux cases. The H1
campaign is another 14 TB. `re2000_Pr_0.025_1` has no isoflux counterpart — its
H1 case was never run.

Read `newcases_direct/README.md` first: it carries the campaign table, the
measured `g_i`, the known warts and the convergence notes for **both**
campaigns. Each case directory then has its own `README.md` (its `g_i`,
`phi_tau`, diffusive times, link arithmetic, column map) and `LAUNCH.txt` (job
ids, start and final time, smoke results).

## The database

A case directory holds four kinds of output.

**`convvelo.<i>.bin`** — the accumulated spectral statistics, one file per
window, with two text sidecars each:

- `convvelo.bin.fields` names what is in them. Every case here is
  direct-only (`output_mode = none, direct = true`), so that is
  `mean_u mean_v mean_w mean_t` / `u_cross_u v_cross_v w_cross_w u_cross_dtu
  v_cross_dtv w_cross_dtw` / `t_cross_t t_cross_dtt` — the dt numerators and
  the energies they divide by, and nothing else.
- `convvelo.<i>.bin.timing` carries `n_triplets`, `n_dropped_triggers`, the
  step sizes, and `energy_rate` per component.

The convection velocity per mode is
`c(kx, kz, y) = -Im<u* du/dt> / (kx <u* u>)`, undefined at `kx = 0`.
`Im<u* du/dt>` is stored in the **imaginary** slot of the dt fields; read
`.imag` from them. `postpro/CONVVELO_DIRECT.md` has the estimator, its
frequency response and why `deltat` is not frozen across a triplet.

**`Dati.cart.<i>.out`** — full spectral snapshots, 606 of them across the six
cases, `(ny+3, 2nz+1, nx+1, 3+nPhi)` complex doubles in Fortran order after a
68-byte header (`nx ny nz` as int32, then `alfa0 beta0 ni a ymin ymax time` as
float64). The component index is slowest, so one component is a contiguous
block and the mean profile of scalar *i* is one contiguous y-line at
`(kx, kz) = (0, 0)` — 503 numbers out of a 50 GB file.
`postpro/scalar_mean_profiles.py` does exactly that and is the cheap way to
look at mean profiles.

**`Runtimedata.phi`** — per step, `1 + 8·nPhi` values per record:

    1              time
    2 .. 1+n       dphi/dy at y0        imposed constant under H2
    2+n .. 1+2n    dphi/dy at yN        imposed constant under H2
    2+2n .. 1+3n   INT(phi) dy          pinned at 2.0 by meantb
    2+3n .. 1+4n   corrtx + meantx      pinned at 2*alpha*g by the wall flux
    2+4n .. 1+5n   <phi> at y0          free under H2
    2+5n .. 1+6n   <phi> at yN          free under H2
    2+6n .. 1+7n   <phi'^2> at y0       free under H2, identically 0 under H1
    2+7n .. 1+8n   <phi'^2> at yN       free under H2, identically 0 under H1

**Parse it by token, not by line.** It is written list-directed, so one record
wraps over several physical lines, and `awk '{print $20}'` silently reads the
wrong column. **It retains only the newest link**, not the whole run — the
whole-run signal is in the snapshots and the convvelo windows.

**`Runtimedata`** — the velocity side, 11 values per record;
`u_tau = sqrt(ni·|du/dy|_wall)` from columns 2 and 3.

### The pipeline

`postpro/build_dataset.py` turns a case's convvelo files into `convvelo.zarr`,
and `postpro/evaluate_convvelo.py <case-dir>` reads that and inverts the
estimator. `aggregate_convection_velocity` defaults to `weighting="leastsquares"`
(`c_LS = Σ kx² E c / Σ kx² E`) rather than energy weighting, because the
deconvolution amplifies noise as `1/kx` and the low-`kx` modes are the noisiest.
That choice is argued in the function's own docstring; keep it.

A run of this on HoreKa needs a long job, not a login node: the source
campaign's own processing used `submit.sh` with a 48 h `large`-partition job.

## Constants you will need

`phi_tau_i = g_i/(Re_tau·Pr_i)` is the one new constant H2 needs — multiply by
it to get Hetsroni's `theta+`. The `g_i` are each case's `[scalars] t0`, and the
measured `Re_tau` are 178.7–178.9, 499.6–500.1 and 999.2, giving
`T_eddy = h/u_tau` of 15.65–15.67, 18.17–18.18 and 20.02.

| case | sampled TU | T_eddy | window starts at `t_start` |
|---|---|---|---|
| re180_Pr_0.025_1 | 9233 | 590 | 16863.3 |
| re180_Pr_1_7 | 4722 | 301 | 12069.3 |
| re180_Pr_7_30 | 2199 | 140 | 4189.1 |
| re500_Pr_0.025_1 | 4210 | 232 | 8442.4 |
| re500_Pr_1_2_4 | 1599 | 88 | 2920.4 |
| re1000_Pr_0.025_1 | 1755 | 88 | 3925.2 |

The warm-up (5 eddy turnovers) is already excluded: convvelo never saw it.

## Caveats

**1. Weight windows by `n_triplets`.** Each `convvelo.<i>.bin` is written as an
*average over its own samples*, not a raw sum — `average_convvelo_field` divides
by `n_field_samples` before writing. So combining windows with equal weight is
wrong. It matters here because **every case has exactly one window holding a
single triplet**: the periodic-write trigger is measured from `t_start`, and
with `t_start > 0` mid-run it fires on the first sample instead of arming from
it. Six such files across the campaign. Either weight by `n_triplets` or drop
any window whose sidecar reports `n_triplets = 1`. No samples are lost either
way. `build_dataset.py` already n-weights `mean_h` and `energy_rate` in
`_read_convvelo_timing_sidecars`; the spectral fields are stacked and the
weighting is left to whoever averages them.

**2. The near-wall exponents are not in the convvelo output.** Direct-only keeps
`t_cross_t` and `t_cross_dtt` but **no u–theta cross terms**, so
`theta'u'`, `theta'v'`, `alpha_t` and `Pr_t` — the rows of the H2 signature that
change exponent — have to come from the 606 field snapshots. Collecting them
online would have meant `output_mode = minimal`, a memory decision deliberately
not taken. `theta_rms` at the wall *is* available, both from `t_cross_t` (which
accumulates the wall plane, `ny0-2 .. nyN+2`) and per step from
`Runtimedata.phi`.

**3. Statistical uncertainty is 0.5–1.3% on wall statistics, and the signal
wanders.** The wall variance moves on roughly **60 eddy turnovers** with ±1–4%
excursions, the same pattern in all scalars of a case at once while `u_tau`
moves ±0.1% — it is the flow realisation, not per-scalar physics. Consequences:
a batch-means error estimate only converges when the batch exceeds that, so a
window shorter than ~60 T_eddy cannot give a trustworthy error bar; and an
autocorrelation estimate truncated at the first zero crossing is optimistic.
Both campaigns are statistically steady — tested on `energy_rate`, which is
zero in a steady state, across 12 case-campaigns with max |t| = 1.60, all below
2 — so this is sampling noise, not a transient.

**4. Per-mode spectra carry real error at the lowest `kx`.** From the H1
campaign's own convergence section: by wavenumber the first-half/second-half
error is **largest at the smallest kx** and falls monotonically — 6.6% median at
kx = 0.2 against 0.1% at kx = 72; by energy, modes below 1e-4 of peak have a p90
of 23% with a tail to 400%. Aggregates are converged (c_LS shifts 0.01–0.06%
between all data and half); raw `c(kx, kz)` plots should be masked or binned,
not plotted raw. Expect the same for H2.

**5. Two pair cross-checks are still open, and they are the campaign's free
acceptance test.** Three (Re_tau, Pr) pairs appear in two independent cases
each — Re_tau 180 at Pr 1 and Pr 7, Re_tau 500 at Pr 1. The two members differ
only in which other scalars ride along, and passive scalars do not interact, so
converged statistics must agree. Measured on the per-step record they are 0.04%,
0.17% and 1.39% apart, but only the first is a tight test: the others compare
windows of 13–39 T_eddy, below one wander time. **Close them on the full
accumulated convvelo statistics**, where the whole window is available. Expect
agreement at the 0.5–1.3% level; anything beyond that is a bug, not physics.

**6. `re1000_Pr_0.025_1`'s imposed gradients are rounded to three decimals**
(5.361 / 30.717 / 51.701) where the five later decks carry four, so its `g` sits
0.0013% off its twin's measured mean. Immaterial physically, but it is why that
case alone looks wrong if you check the enforced gradient against the *measured*
values rather than against its own deck. Against its deck it holds to 4.3e-14
like every other case.

**7. One link was lost to a GPU fault and retried.** `re180_Pr_7_30`'s job
5184333 died mid-run with `CUDA_ERROR_UNKNOWN`; `dt_save = -1` means no restart
was written, so nothing is corrupted and the next link repeated the interval.
The case has its full seven advances. Nothing to do; noted so the gap in the
job ids does not look like missing data.

## Working off HoreKa

What is small enough to pull to a workstation:

- every `*.timing` sidecar and `*.bin.fields` — kilobytes, and enough for the
  steady-state test and the window bookkeeping;
- `Runtimedata` and `Runtimedata.phi` — tens of MB per case, the newest link
  only;
- mean profiles via `postpro/scalar_mean_profiles.py` — 8 KB per scalar per
  snapshot, so a whole case's mean-profile history is a few hundred KB;
- the case `README.md` and `LAUNCH.txt` files.

What is not: the `convvelo.*.bin` files (2.1 to 101.4 GB each, 77 of them) and the
`Dati.cart.*.out` snapshots (1.1 to 50.7 GB each, 606 of them). Those are
processed where they sit. `sync_to_horeka.sh` in the parent directory pushes
cases *to* HoreKa and excludes `*.zarr/`; there is no pull script.

## What to produce

1. `convvelo.zarr` and the aggregate convection velocity for all six isoflux
   cases, next to the H1 values already computed.
2. **Hetsroni's `U_cq+`** — the convection velocity of the wall temperature
   fluctuation — which is the result the campaign exists for, and which has no
   H1 counterpart because under a Dirichlet wall that fluctuation is zero.
3. The near-wall signature from the snapshots, against H1 at the same Re and Pr:
   `theta'u' ~ y` against `y²`, `theta'v' ~ y²` against `y³`, `Pr_t → 0`
   against a constant ≈ 1.1, and the thermal streak spacing ≈ 140 wall units
   against ≈ 100.
4. The pair cross-checks of caveat 5, closed properly.

For reference, `theta_rms` at the wall from the per-step record is already known
and is the campaign's headline check: 1.96 in wall units at Pr = 0.7, against
the ≈ 2.0 Kong and Hetsroni report, and identically zero under H1. It runs
monotonically from 0.010 at Pr = 0.005 to 20.4 at Pr = 30.
