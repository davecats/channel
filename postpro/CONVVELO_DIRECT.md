# The direct convection velocity

Reference for the online estimator added in `8dfb13f` and the post-processing
added in `1bb3f3e`. Enabled with `[convvelo] direct = true`; see the README for
the input keys.

Two independent routes to the same quantity live in this tree:

- **reconstruction** — `du/dt` rebuilt term by term from the momentum equation,
  from ~30 accumulated cross-spectra. Temporally exact, but model-dependent:
  it hardwires `t_forcing = 0`, `minimal` mode drops the mean-crossflow terms,
  and almost every term needs a y-derivative.
- **direct** — `du/dt` from a three-point centred difference over consecutive
  timesteps. Assumes nothing about the right-hand side and needs no spatial
  operator, but is band-limited in frequency.

They are complementary, and where both are valid they agree: measured 0.68% for
`y+ >= 2` on `re1000` (Re_tau = 1009, 1023x500x512).

## What is computed

Per mode, for u, v, w and every scalar:

    c(kx, kz, y) = -Im<u* du/dt> / (kx <u* u>)

`kx = 0` carries no streamwise phase and is undefined there; the denominator is
masked so the result is NaN rather than the `+-inf` an unmasked division gives
(the numerator does not vanish at `kx = 0`).

`output_mode = none, direct = true` stores exactly the two ingredients —
`u_cross_dtu` and `u_cross_u` per component — and nothing else: 22 spectral
planes per rank against `minimal`'s 49 and `full`'s 64.

## The estimator

With `h1 = t_n - t_(n-1)` and `h2 = t_(n+1) - t_n`, the second-order derivative
at the middle node is the exact Lagrange one:

    du/dt|n = a u_(n-1) + b u_n + c u_(n+1)
    a = -h2/(h1(h1+h2))    b = (h2-h1)/(h1 h2)    c = h1/(h2(h1+h2))

collapsing to `(u_(n+1) - u_(n-1))/2h` when `h1 = h2`.

**`deltat` is deliberately not frozen across a triplet.** The quadratic moment
`a h1^2 + c h2^2` vanishes by construction, so the operator stays second order
at `h1 != h2`. Expanding the response for `u = A exp(-i omega t)`:

    G = -i omega [1 - omega^2 h1 h2/6 + ...] + omega^4 h1 h2 (h2-h1)/24 + ...
         |__ imaginary: this is c __|          |__ real: only this feels h1 != h2 __|

The asymmetry first enters at the *fourth* moment and only in the real part,
where `c` does not live; the imaginary part tracks `sinc(omega sqrt(h1 h2))` to
~1e-4 even at 20% asymmetry. Measured asymmetry under CFL control is 0.37%
median. Freezing `deltat` would let a diagnostic perturb the timestepping and
would cost the bit-exactness gate, for a fourth-order correction to a quantity
that is not stored.

Only `Im<u* du/dt>` is kept — it is the whole numerator — so it is accumulated
into a real array and written into the **imaginary** slot of a complex field.
Post-processing therefore reads `.imag` from a dt field exactly as it would from
a full complex accumulator, and the 24-byte header and `.fields` layout are
unchanged.

## Band limit, and how to undo it

For a travelling mode the difference returns, exactly,

    du/dt|est = -i sin(omega h)/h u,    omega = kx c

so the measured velocity is `c sinc(omega h)`. Since `kx h c_meas = sin(omega h)`
this inverts in closed form,

    c = arcsin(kx h c_meas) / (kx h)

and it is self-diagnosing: where `|kx h c_meas| > 1` there is no solution, the
mode is temporally aliased, and `evaluate_convvelo` masks it rather than
correcting it (equivalently `omega h > pi/2`). Exact for a monochromatic mode;
real turbulence carries a band of frequencies per `(kx, kz, y)`, so this is a
leading-order correction. The direction of the bias and the location of the
failure are robust; the corrected value is not exact.

At production parameters the correction is small. Measured over 1.1e6 modes of
`re1000` at `deltat = 0.00233`: every mode resolved (`omega h` median 0.038, max
0.662 against the `pi/2` limit), aggregate bias 0.001% energy-weighted and
0.145% least-squares-weighted. The worst single high-`kx` mode is ~7% low
uncorrected. Use `*_uc_direct_corrected`; it costs nothing.

## Aggregating over modes

`aggregate_convection_velocity(res, component, weighting=...)` gives the two
standard collapses, both computable from what a direct-only file stores:

    energy         c_E  = sum(E c)/sum(E)         = sum(-Im<u* du/dt>/kx) / sum(<u* u>)
    leastsquares   c_LS = sum(kx^2 E c)/sum(kx^2 E) = sum(-kx Im<u* du/dt>) / sum(kx^2 <u* u>)

`c_LS` is the aggregate form of `-<du/dt du/dx>/<(du/dx)^2>` and is the default
because it is better conditioned: `c` carries a `1/kx`, so the low-`kx` modes are
the noisiest and energy weighting gives exactly those the most weight. The
`kx^2` cancels that. `kx = 0` is excluded from both sums. Deconvolve per mode
*before* weighting — each mode has its own `sinc`, and one factor cannot stand
for a mixture — which is why `corrected=True` is the default.

## Near the wall

This is the practical reason to prefer the direct estimator. The reconstruction
needs y-derivatives for nearly every term (`u_cross_dyyu` viscous,
`u_cross_dyuv` nonlinear, `v_cross_dpdy` pressure, plus the mean gradients).
Near the wall those terms are individually large and nearly cancel, so the
residual is a small difference of big numbers taken with one-sided stencils on a
stretched mesh. The direct estimator differences in time at each `(kx, kz, y)`
independently: no y-derivative, no coupling between y levels, no cancellation.

Measured, `re1000`, same kx/kz sample and the same least-squares aggregation
applied to both:

| iy | y+ | reconstruction | direct |
| --- | --- | --- | --- |
| 0 | 0.00 (wall) | **8.9e13** | 0.0063 |
| 1 | 0.98 | 0.4037 | 0.4584 |
| 2 | 1.97 | 0.4548 | 0.4510 |
| 5 | 5.01 | 0.4597 | 0.4561 |
| 49 | 64.93 | 0.7417 | 0.7365 |
| 249 | 1002.19 | 1.1031 | 1.1062 |

The reconstruction diverges at the wall, where `u = 0` makes its denominator
vanish. At the first point off the wall it is 11.9% low and the wrong shape: its
first differences over `y+ 0.98 -> 5` run `+0.051, -0.005, +0.005, +0.006` —
two sign changes, a kink. The direct profile runs `-0.007, -0.002, +0.002,
+0.005` — one sign change, a single shallow minimum near `y+ ~ 3`.

So **the direct near-wall data needs no y-averaging, no near-wall exclusion and
no smoothing.** The temporal resolution also improves toward the wall, because
`omega h = kx c deltat` and `c` is smallest there — 0.27 at `y+ = 1` against
0.66 at the centreline.

Two rows must still be dropped, and neither is a near-wall exclusion:

- `iy = -1` and `iy = ny+1` are **ghost rows outside the channel**.
  `build_dataset.py` includes them in the `y` coordinate, and they return
  *plausible* values (0.49 on re1000), not obvious garbage — which makes them
  dangerous in a y-integral.
- `iy = 0` and `iy = ny` are the **walls**, where `c` is `0/0`.

Usable range is `iy = 1 .. ny-1`, i.e. `.isel(y=slice(2, -2))`.

## Sampling and statistical convergence

Each `dt_compute` trigger becomes a three-step burst. The trigger step is taken
as `n-1`, so the sample centres one `deltat` later than the nominal trigger —
negligible against `dt_compute`, and it avoids predicting a trigger one step
ahead, which adaptive `deltat` makes unreliable. The reconstruction statistics
move to the centre step so both estimates share a time and a denominator. A
`dt_write` landing mid-triplet is deferred until the triplet closes.

`n_dropped_triggers` in the `.timing` sidecar counts triggers that fired while a
triplet was open. At production settings (`dt_compute/deltat ~ 80` against the 3
a triplet needs) these are **duplicate sightings, not lost samples**: the
interval detector tests `floor((t +- deltat/2)/period)`, so when `deltat` shrinks
between steps the two windows overlap and one boundary is seen twice. The second
is correctly ignored. Only worry if the count approaches the triplet count.

Snapshots are **independent averaging windows** carrying their own sample count,
combined by a weighted mean (`_average_convvelo_steps`). Windows are contiguous
and non-overlapping. **Never average snapshots with equal weight** — a case
holds a mix of full windows (~1000 triplets) and short run-end remainders (~10).

From the 2026-10 campaign, splitting each case's snapshots in half:

- **aggregates converge fast**: least-squares `c` at mid-channel shifts
  0.01-0.06% between half the data and all of it.
- **per-mode spectra converge slowly at both ends**. By energy: above 1e-2 of
  peak, ~1% median / 4% p90; below 1e-4 of peak, 23% p90 with a long tail. By
  wavenumber the error is *largest at the smallest kx* and falls monotonically —
  6.6% median at `kx = 0.2`, 2.2% at `kx = 1`, 0.08% at `kx = 256` — because the
  longest structures decorrelate slowest. That limit is set by simulated time,
  not by snapshot count.

## Validation

- `test_convvelo_direct` (1 and 2 ranks) drives the estimator with
  `u = A exp(-i omega t)` at `h1 != h2` and checks the accumulated numerator
  against the exact `|A|^2 Im(a e^(i w h1) + b + c e^(-i w h2))`.
- The direct path only *reads* the field, so `direct = true` must leave
  `Dati.cart.out` byte-identical to `direct = false`. That is the right gate
  here, and it holds.
