# Prompt for the next session

Copy the block below as the opening message.

---

Continue simplifying `channel` for readability, without changing behaviour.

Read `.claude/projects/-home-ws-xt8786-Codes-channel/memory/MEMORY.md` first.
It points at the verification gates (including which machines have a usable
GPU, how to reach them, and what is currently broken on the local one), the
NVHPC constraints on module structure, a separate NVHPC trap about device
calls, and the state of the layering plan. All of it was learned the hard way
and will save you several wasted build cycles.

## Where the previous session left it

Four code commits, each independently verified against the full gate.

- **The convvelo MPI-IO writer is out.** `src/core/convvelo_io.f90` (229
  lines, plain `.f90`) holds the raw-statistics writer, the profile write, the
  field averaging and the `.fields` layout file. It takes its state as
  arguments rather than `use`ing convvelo -- the shape `restart_io` already
  had -- so the dependency runs one way and it needs no field catalogue, only
  the names its caller passes it. `convvelo.fypp` 1127 -> 913 lines.
  `write_convvelo_raw_stats` is still the public entry point (two tests and
  the snapshot call it); it now does the two device updates and hands the host
  arrays over.
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
  label; `sendbuf`/`recvbuf` stay behind the boundary. `transform_to_physical`
  keeps its four steps apart on purpose (it overlaps the alltoall with a whole
  component), and `mpi_autotune` was left alone because its copy sits inside
  the loop it is timing.
- **The single-field xz round trip is shared.**
  `spectral_field_to_real_x` / `real_x_to_spectral_field` are now one copy in
  `channel_transforms`, sitting next to the pipelined round trip they are the
  non-overlapped sibling of; the profiler label is an argument. The
  `pressure_output` form is the one that survived -- `default(none)` kept, and
  convvelo gained the `y_first`/`y_last` hoist that satisfies it -- rather
  than dropping the clause from the stricter side. The x-padding zero turned
  out to be a third copy of `zero_vvdx_hft` and now calls it.
  - **The `VVdx`/`VVdz` obstacle the last handoff warned about was stale.**
    They have lived in `ffts` on *both* the GPU and the FFTW path since they
    stopped being threaded through `init_fft`; there is a comment saying so at
    `ffts.f90:67`. So the shared routines need no `#if` and both callers lost
    theirs, and convvelo and pressure_output no longer `use mpi_transpose` at
    all.
  - `load_convvelo_field_to_zbuf` vs `load_pressure_field_to_zbuf` differ much
    more than these did; they were left alone.

## Do this first, roughly in value order

1. **The MPI/NCCL pipeline pair in `y_line_solvers`** (~460 lines written
   twice). Genuinely different code -- one non-blocking with overlap, one
   blocking -- so unifying it would change behaviour. Read both before
   deciding, and be willing to report that it should stay.
2. **`FFT`/`IFT` and `RFT`/`HFT` in `ffts.f90`.** Each pair differs only by
   plan and direction across three `#ifdef` backends. Judge whether the result
   is actually more readable before committing to it.

## Looked at and deliberately left alone

Do not redo these without new information; the reasoning is in the commits.

- **The test-support block in `convvelo`** (`fill_convvelo_synthetic_state_
  for_test`, `synthetic_profile_value`, `synthetic_field_value`, ~55 lines,
  public, shipping in the production module). Every shape tried was worse: the
  filler writes `n_field_samples`, `n_mean_samples` and the averaging window,
  all private, so a test-side module means making them public; and its two
  `!$omp target update to(...)` calls should stay with the declarations (see
  [[gpu-build-constrains-module-structure]] -- a module mapping *another*
  module's variables has failed here before). The least-bad version trades one
  test-only public routine for another. Worth doing only with a shape that
  shrinks the public surface rather than reshuffling it. The synthetic pattern
  is duplicated in `tests/convvelo/verify_synthetic_raw.py` either way -- that
  copy is across languages and cannot be removed.
- **`mpi_autotune`'s copy of the blocking transpose.** Converting it would put
  two extra roctx markers inside the loop the autotuner times, and that
  measurement feeds the decomposition choice.
- **The halo unpack pair in `y_line_solvers`** -- the two differ by
  destination array, not by an index; the reason is a comment at the site. The
  *pack* pair was unified.
- **Reconstruct (step 3 of the wall closure)** stays written twice: the host
  spells it `a - sum(b)`, the batched one `a - b1 - b2 - b3`, and those do not
  round alike.

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
- **The autotuner tries NCCL even when ranks share one GPU.** Arguably it
  should skip NCCL when it sees a shared device -- but that is a behaviour
  change, so ask.

## How to work

Small, independently verified commits -- one logical change each, so a
regression bisects to one peel. The gate, none of which is optional:

1. CPU `ctest` (30 tests) -- necessary, never sufficient.
2. An NVHPC GPU build, and the GPU suite compared against **HEAD built in a
   scratch tree**, not against an expectation. As of Aug 2026 the local box
   fails 17 of 30 at HEAD -- a UCX/CMA crash in multi-rank runs sharing the
   one device, up from 5 the session before, so the environment moves
   underneath you. `module load toolkits/nvhpc/25.9` silently does nothing
   there; a fresh configure needs
   `PATH=$NV/compilers/bin:$NV/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH`
   with `NV=/opt/Nvidia/nvhpc/Linux_x86_64/25.9`, plus
   `-DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_Fortran_COMPILER=$NV/compilers/bin/
   nvfortran` -- that recipe works, and the same 17 fail by name in the
   scratch tree, so the comparison is worth the one extra build.
   `istmcetus` is the clean machine.
3. `-Dbodyforce` compiles (uncomment it in `src/core/header.h`; nothing builds
   it, so it hides breakage).
4. **Run the `channel` binary and check `$?`.** `ctest` never runs it -- every
   test drives its own `test_*` binary -- which is how it segfaulted in
   `MPI_Finalize` on every run with 30 green tests.
5. Bit-identity against the parent commit, both trees built with
   `FFTW_PATIENT` patched to `FFTW_ESTIMATE`. The working recipe: `git archive
   HEAD` into a scratch dir, `sed` the plan type, build, then run
   `mpirun -np 2 .../channel` in a directory holding
   `tests/convvelo/dns_runtime.in` as `dns.in` and
   `tests/data/start_field_scalar.out` as `Dati.cart.out`, and `cmp` the
   `convvelo.*.bin`, `Dati.cart*.out` and `Runtimedata*` outputs. That deck
   drives convvelo, pressure_output and the solver in one run. `run.log`
   always differs -- autotune timings and rank interleaving.

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
  alltoall path is never entered -- so a bit-identity check of the transpose
  work proved nothing about the code it changed until it was rerun with
  `CHANNEL_NPXZ=2 CHANNEL_NPY=1`. Same lesson as `CHANNEL_OVERLAPPING=1` and
  `CHANNEL_YS_FORCE_CUSTOM_GPSV`, which no test sets and which was hiding a
  GPU crash: "30/30 green" is a statement about the 30.
- **Check that a flag you are relying on actually took effect.** Two claims
  went wrong this way. `CHANNEL_OVERLAPPING=1` is bit-identical by design, so
  output cannot confirm it -- the memory line does (buffers double). And a
  comparison meant to isolate one flag had a second flag varying with it.
- **Compare against the checked-in references too**, not only base-vs-work.
- **Re-check the obstacle before working around it.** This handoff said a
  shared transform module would need `VVdx`/`VVdz` behind an `#if` because
  they came from two different modules. They had not for some time -- one
  `grep` retired the whole complication. Notes about *why something is hard*
  age faster than notes about what was done.
- `fprettify` corrupts `.fypp` (it eats the space in `${macro}$ - x`); the
  pre-commit hook is restricted to `\.f90$`. `pre-commit` is not on PATH on
  this box -- run `fprettify -i 2 -w 2` from a throwaway venv instead.

Report honestly when something is reverted and why. Several extractions have
been abandoned across these sessions for good reasons, and recording that has
been more useful than forcing them through.
