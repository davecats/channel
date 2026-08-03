# Prompt for the next session

Copy the block below as the opening message.

---

Continue simplifying `channel` for readability, without changing behaviour.

Read `.claude/projects/-home-ws-xt8786-Codes-channel/memory/MEMORY.md` first.
It points at the verification gates (including which machines have a usable
GPU and how to reach them), the NVHPC constraints on module structure, a
separate NVHPC trap about device calls, and the state of the layering plan.
All of it was learned the hard way and will save you several wasted build
cycles.

## Where the previous session left it

Three commits, each independently verified.

- **The convvelo MPI-IO writer is out.** `src/core/convvelo_io.f90` holds the
  raw-statistics writer, the profile write, the field averaging and the
  `.fields` layout file. It takes its state as arguments rather than `use`ing
  convvelo -- the shape `restart_io` already had -- so the dependency runs one
  way and it needs no field catalogue, only the names its caller passes it.
  `convvelo.fypp` 1127 -> 913 lines. `write_convvelo_raw_stats` is still the
  public entry point (two tests and the snapshot call it); it now does the two
  device updates and hands the host arrays over.
  - `normalize_convvelo_field_for_output` and `copy_convvelo_field_average`
    were the same loop; both now call `average_convvelo_field`.
  - The mean profiles were three copy-pasted blocks plus a loop over the
    scalars; the displacement arithmetic works out to one loop over 3 + nPhi.
  - `test_convvelo_utils` had its own copy of the header size and now uses
    convvelo_io's.
- **The convvelo transform routines asked `fft_transpose_is_local` three times
  in a row.** The wait and the unpack only apply to the alltoall the first
  `if` just posted, so they moved into that branch.
- **The blocking xz transpose is one routine.** pack / post / wait / unpack --
  or repack in place -- was written out at six sites. `mpi_transpose` now has
  `transpose_zTOx` / `transpose_xTOz`, taking the two arrays and a profiler
  label; the buffers stay behind the boundary. convvelo and pressure_output
  each went from a fifteen-name `use mpi_transpose` to two names.
  `transform_to_physical` keeps its steps apart on purpose (it overlaps the
  alltoall with a whole component) and `mpi_autotune` was left alone because
  its copy is inside the loop it times.

## Do this first

Nothing in convvelo is as clean a peel as the writer was. Judgement calls,
roughly in value order:

1. **The test-support block in `convvelo`** (`fill_convvelo_synthetic_state_
   for_test`, `synthetic_profile_value`, `synthetic_field_value`, ~55 lines,
   public, shipping in the production module). The previous session looked at
   moving it and **did not**, because every shape it tried was worse: the
   filler writes `n_field_samples`, `n_mean_samples` and the averaging window,
   all private, so a test-side module means making them public; and the two
   `!$omp target update to(...)` calls should stay with the declarations (see
   [[gpu-build-constrains-module-structure]] -- a module mapping *another*
   module's variables has failed here before). The least-bad version is a
   public setter plus host writes through the already-public
   `convvelo_stats`/`component_means`, which trades one test-only public
   routine for another. Worth doing only if you find a shape that shrinks the
   public surface rather than reshuffling it. The synthetic pattern is
   duplicated in `tests/convvelo/verify_synthetic_raw.py` either way -- that
   copy is across languages and cannot be removed.
2. **`convvelo` and `pressure_output` now have the same three routines.**
   `spectral_field_to_real_x` and `real_x_to_spectral_field` are, after the
   transpose helper, near-identical in the two files -- pressure_output's pair
   differs only by hoisting `ny0`/`nyN` into `y_first`/`y_last` locals and by
   `default(none)` on the target regions. `load_*_field_to_zbuf` differs more.
   The obstacle is `VVdx`/`VVdz`: they come from `ffts` on GPU and from
   `dnsdata` on CPU, so a shared module would need both spellings behind the
   same `#if`, which is how both files already open. Worth a look; the
   `default(none)` difference has to be resolved deliberately, not silently.
   Note the transform round trip in `channel_transforms` is *not* a third copy
   -- it is the pipelined version (`to`/`from` indices, `requests(m)`, a
   component-ahead loop) and unifying it would change behaviour.
3. **The MPI/NCCL pipeline pair in `y_line_solvers`** (~460 lines written
   twice). Genuinely different code -- one non-blocking with overlap, one
   blocking -- so unifying it would change behaviour. Read both before
   deciding, and be willing to report that it should stay.
4. **`FFT`/`IFT` and `RFT`/`HFT` in `ffts.f90`.** Each pair differs only by
   plan and direction across three `#ifdef` backends. Judge whether the result
   is actually more readable before committing to it.

## Open items that are the user's call, not yours

- **`FFTW_PATIENT`** (`ffts.f90`) makes CPU runs non-reproducible at ULP level.
  Decided: keep it, and keep patching a scratch tree to `FFTW_ESTIMATE` for
  bit-identity checks. Do not re-litigate.
- **`-Dbodyforce` is dead code and is being kept anyway.** `F` is allocated,
  zeroed, and never assigned, so a `-Dbodyforce` build is byte-identical to a
  default one. The user decided to keep it rather than delete it. **It still
  needs validating**: a source for `F` (deck section or file read) plus a test
  that a known forcing produces the expected response. Until then the forcing
  terms in `buildrhs_prepare` have never been exercised with a non-zero `F` --
  do not trust them. See [[bodyforce-is-dead-code]].
- **`-Dhalfchannel` and `-DphiNeumann`** do take effect (both change the
  output), but nothing checks their numbers against a reference.
- **17 of the 30 GPU tests now fail on the local workstation**, at HEAD and
  before it -- a UCX/CMA crash in multi-rank runs that share the one device,
  not a code bug. It was 5 in the previous session, so the environment moved.
  Every single-rank test and every convvelo test still passes, including the
  4-rank `convvelo_mpi_io_npxz2_npy2`. Establish the failing set at HEAD
  before reading anything into a GPU run here, or use `istmcetus`.
- **The autotuner tries NCCL even when ranks share one GPU.** Arguably it
  should skip NCCL when it sees a shared device -- but that is a behaviour
  change, so ask.

## How to work

Small, independently verified commits -- one logical change each, so a
regression bisects to one peel. For every change: CPU `ctest` (30 tests), an
NVHPC GPU build, `-Dbodyforce` compiles, **run the `channel` binary and check
its exit status**, and bit-identity against the previous commit using a
`FFTW_ESTIMATE` scratch build. A green CPU suite has repeatedly been
insufficient.

Techniques worth reusing, in order of how much they save:

- **Diff the preprocessed or generated Fortran, not just the output.** If a
  change is fypp expansion or conditional removal, expanding both trees and
  diffing proves the compiled program is unchanged. That retires the
  performance question without needing an idle GPU. Ignore blank lines: `cpp`
  leaves one where it removes a directive.
- **Shape macros to take the *text* of each call site's own expressions**
  rather than fixed variable names. That is what makes the generated code
  identical to what was there before -- and, as the GPU fix showed, it is also
  the way to keep device code free of dummy-argument indirection.
- **Ask which branch your run actually took.** A 2-rank run of the convvelo
  deck autotunes to npxz=1, where `fft_transpose_is_local` is true and the
  alltoall path is never entered -- so the first bit-identity check of the
  branch merge above proved nothing about the code it changed. `CHANNEL_NPXZ=2
  CHANNEL_NPY=1` forces the other side. Same lesson as `CHANNEL_OVERLAPPING=1`
  and `CHANNEL_YS_FORCE_CUSTOM_GPSV`: "30/30 green" is a statement about the
  30.
- **Check that a flag you are relying on actually took effect.** Two claims
  went wrong this way. `CHANNEL_OVERLAPPING=1` is bit-identical by design, so
  output cannot confirm it -- the memory line does (buffers double). And a
  comparison meant to isolate one flag had a second flag varying with it.
- **Compare against the checked-in references too**, not only base-vs-work.

Report honestly when something is reverted and why. Several extractions have
been abandoned across these sessions for good reasons, and recording that has
been more useful than forcing them through.
