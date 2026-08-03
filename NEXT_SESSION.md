# Prompt for the next session

Copy the block below as the opening message.

---

Continue simplifying `channel` for readability, without changing behaviour.
The next step is peeling the MPI-IO writer out of `convvelo`.

Read `.claude/projects/-home-ws-xt8786-Codes-channel/memory/MEMORY.md` first.
It points at the verification gates (including which machines have a usable
GPU and how to reach them), the NVHPC constraints on module structure, a
separate NVHPC trap about device calls, and the state of the layering plan.
All of it was learned the hard way and will save you several wasted build
cycles.

## Where the previous session left it

Seven commits, each independently verified.

- **The halo pack pair** in `y_line_solvers` is one routine. The *unpack* pair
  stays written twice on purpose -- they differ by destination array, not by an
  index -- with the reason recorded at the site. The handoff had this filed as
  a quartet; only half of it was index dispatch.
- **`channel` no longer segfaults on exit.** `driver`'s `finalize` was calling
  the old `mpi_finalize_` binding with no `ierror` argument, on every run, at
  every rank count. `ctest` never runs the `channel` binary, which is why 30
  green tests hid it. Gate leg 5 now says to run it and check `$?`.
- **NOMPI is gone.** It had not compiled for some time; MPI is now required and
  all 88 `HAVE_MPI` conditionals are deleted (288 lines). Verified by showing
  every touched file preprocesses identically before and after, under three
  backend flag sets -- stronger than a bit-identical run.
- **The convvelo field catalogue is defined once** in `convvelo_fields.fypph`.
  It had been written out twice, and `acc_convvelo_stats` looked every field up
  by string at runtime: 43 linear scans per accumulation, all of literals.
- **`dnsdata` is pure setup** (528 -> 350 lines). The transform pass moved to
  `channel_transforms`.
- **The custom pentadiagonal solve works on GPU.** It had been dying with
  `CUDA_ERROR_ILLEGAL_ADDRESS` since before the session. Root cause: nvfortran
  cannot forward an assumed-shape *dummy* array from inside a target region
  into a `declare target` routine. `y_line_solvers.f90` -> `.fypp`; the solve
  is expanded once per array set. See [[nvfortran-dummy-array-device-call]].

## Do this first

**`convvelo` still mixes six concerns.** The ~180-line MPI-IO writer
(`write_convvelo_raw_stats` and its helpers, from roughly the four
`MPI_Type_create_subarray` calls to `write_convvelo_field_layout`) is the
cleanest remaining peel: self-contained, narrow interface, and four dedicated
tests plus the regression pair cover it. `convvelo` has no `declare target` of
its own, so a new module that also declares none is the safe shape.

Two smaller things in the same file, worth folding in or doing next:

- `normalize_convvelo_field_for_output` and `copy_convvelo_field_average` are
  the same loop -- divide a field by its sample count -- differing only in an
  early-out and one device-update call.
- The test-support block (`fill_convvelo_synthetic_state_for_test`,
  `synthetic_profile_value`, `synthetic_field_value`) is public and ships in
  the production module.

## After that, roughly in value order

1. **The MPI/NCCL pipeline pair in `y_line_solvers`** (~460 lines written
   twice). Genuinely different code -- one non-blocking with overlap, one
   blocking -- so unifying it would change behaviour. Read both before
   deciding, and be willing to report that it should stay.
2. **`FFT`/`IFT` and `RFT`/`HFT` in `ffts.f90`.** Each pair differs only by
   plan and direction across three `#ifdef` backends. Judge whether the result
   is actually more readable before committing to it.
3. **The transform round trip is written twice** -- `channel_transforms` and a
   private copy in `convvelo` (`spectral_field_to_real_x`,
   `real_x_to_spectral_field`). Now that the first has its own module,
   unifying them is worth a look. The convvelo copy also spells
   `if (.not. fft_transpose_is_local)` three times in a row where one block
   would do.

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
- **The autotuner tries NCCL even when ranks share one GPU**, which makes 5
  tests fail on any single-GPU box. `CHANNEL_COMM=mpi` works around it.
  Arguably it should skip NCCL when it sees a shared device -- but that is a
  behaviour change, so ask.

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
- **Check that a flag you are relying on actually took effect.** Two claims
  went wrong this way. `CHANNEL_OVERLAPPING=1` is bit-identical by design, so
  output cannot confirm it -- the memory line does (buffers double). And a
  comparison meant to isolate one flag had a second flag varying with it.
- **Compare against the checked-in references too**, not only base-vs-work.
- **When a run passes, ask what it did not cover.** `CHANNEL_OVERLAPPING` and
  `CHANNEL_YS_FORCE_CUSTOM_GPSV` are set by no test; the second was hiding a
  GPU crash. "30/30 green" is a statement about the 30.

Report honestly when something is reverted and why. Several extractions have
been abandoned across these sessions for good reasons, and recording that has
been more useful than forcing them through.
