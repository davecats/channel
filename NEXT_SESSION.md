# Prompt for the next session

Copy the block below as the opening message.

---

Continue simplifying `channel` for readability, without changing behaviour.
This time the subject is the six xz pack/repack kernels in
`src/mpi/mpi_transpose.f90` -- see "Do this first" below. Deciding which of
them, if any, should share a body is the task, and "these should stay written
out, because ..." is a complete answer. What is *not* acceptable is deciding
without reading all six, or unifying them in a way that changes a loop nest
order.

Read `.claude/projects/-home-ws-xt8786-Codes-channel/memory/MEMORY.md` first,
then this whole file. Between them they carry the verification gate (including
which machines have a usable GPU, how to reach them, and the one environment
variable that decides whether half the GPU suite runs at all), two separate
NVHPC traps about module structure and device calls, and the state of the
layering plan. All of it was learned the hard way and will save you several
wasted build cycles.

Work in small, independently verified commits -- one logical change each. Run
the full gate on every one; none of its legs is optional, and each has caught
something the others missed.

## Where the previous sessions left it

Eight code commits, each independently verified against the full gate.

- **The convvelo MPI-IO writer is out.** `src/core/convvelo_io.f90` holds the
  raw-statistics writer, the profile write, the field averaging and the
  `.fields` layout file, taking its state as arguments rather than `use`ing
  convvelo -- the shape `restart_io` already had. `convvelo.fypp` 1127 -> 913
  lines. Folded in on the way: two identical averaging loops became
  `average_convvelo_field`, three copy-pasted mean-profile writes became one
  loop over 3 + nPhi, and `test_convvelo_utils` stopped keeping its own copy
  of the header size.
- **The convvelo transform routines asked `fft_transpose_is_local` three
  times in a row**; the wait and the unpack moved into the branch that posts
  the alltoall.
- **The blocking xz transpose is one routine.** pack / post / wait / unpack --
  or repack in place -- was written out at six sites. `mpi_transpose` now has
  `transpose_zTOx` / `transpose_xTOz`; `sendbuf`/`recvbuf` stay behind the
  boundary. `transform_to_physical` keeps its four steps apart on purpose (it
  overlaps the alltoall with a whole component), and `mpi_autotune` was left
  alone because its copy sits inside the loop it is timing.
- **The single-field xz round trip is shared.**
  `spectral_field_to_real_x` / `real_x_to_spectral_field` are one copy in
  `channel_transforms`; the profiler label is an argument. The x-padding zero
  turned out to be a third copy of `zero_vvdx_hft` and now calls it. The
  `VVdx`/`VVdz` obstacle an older handoff warned about was **stale** -- one
  `grep` retired the whole complication.
- **The pipelined-LU buffer layout is named once**, and **its instrumentation
  is two fypp macros** (both below).
- **`-Dbodyforce` is gone** (below).
- **`FFT` and `IFT` are written once**; `RFT`/`HFT` and the three `free_fft`
  bodies were looked at and left (below).

## Do this first: the six xz pack kernels in `mpi_transpose.f90`

`src/mpi/mpi_transpose.f90:418-551`. Six small `!$omp target teams distribute
parallel do` kernels move the same data between the z-pencil and x-pencil
layouts:

- `repack_zTOx_local` / `repack_xTOz_local` (19 lines each) -- the
  `fft_transpose_is_local` case, a plain transpose with no buffer.
- `pack_zTOx` / `unpack_zTOx` and `pack_xTOz` / `unpack_xTOz` (22-24 lines
  each) -- the alltoall case, one linear buffer index per element.

They are much closer to each other than `FFT`/`IFT` were, but **not by three
tokens**, and the differences are not all cosmetic:

- **The loop nest order differs on purpose.** Every kernel puts the *output's*
  fastest index innermost (`repack_zTOx_local` runs `iz` inner, its mirror
  runs `ix` inner; same in the pack pair). That is a coalescing choice on the
  device, so a shared body must keep each site's own nest order -- swapping
  one is a performance change, invisible to every correctness check in the
  gate.
- **The linear index differs in which extent is the stride**: `iz + nzB*(ix-1)
  + nzB*nxB*(iy-1)` against `ix + nxB*(iz-1) + nxB*nzB*(iy-1)`.
- **The rank offset moves sides**: `pack_*` offsets the *array* subscript
  (`Vz(dest*nzB + iz, ...)`) while `unpack_*` offsets it on the other array
  (`Vx(ix + src*nxB, ...)`), and the buffer subscript carries `dest`/`src`
  either way.

So the honest question is whether a fypp macro taking each site's own nest
order and index text (the technique in "Techniques worth reusing") reads
better than six kernels that each fit on one screen -- or whether naming the
buffer index alone (one `#:def` used six times) gets most of the value. Both
"unify" and "leave them, because ..." are acceptable; deciding without reading
all six is not.

Practicalities:

- `mpi_transpose.f90` is plain `.f90`, so this means **`.f90` -> `.fypp`** and
  moving it from the `add_library(channel_core ...)` list into
  `channel_fypp_generate(...)` (`CMakeLists.txt:105`). That move is routine
  here -- eight files have made it, most recently `ffts`. fypp and `#ifdef`
  coexist fine; `y_line_solvers.fypp` does both.
- **Aim for byte-identical generated Fortran** (see below); it is what retires
  the performance question for a device kernel without an idle GPU.
- These kernels are on the alltoall path, which **a 2-rank run does not reach**
  -- it autotunes to npxz=1. Force `CHANNEL_NPXZ=2 CHANNEL_NPY=1`.
- The 15 GPU tests that fail environmentally are exactly the ones whose xz
  transpose goes over MPI, i.e. the ones that would exercise this code on the
  device. Plan on the CPU suite plus generated-Fortran identity carrying the
  weight, and say so.

**"These should stay written out, because ..." is a perfectly good answer.**
Several extractions have been abandoned across these sessions and recording
the reason has been worth more than forcing them through.

## Settled recently: `FFT` and `IFT`, and what was left in `ffts.f90`

`ffts.f90` -> `ffts.fypp`, with a `z2z_transform(name, direction)` macro in
the preamble. The two 33-line bodies really did differ in three tokens per
backend, and the plan handles already follow the routine name everywhere
(`cu_p<name>` / `hip_p<name>` / `p<name>`), so the call sites are
`$:z2z_transform('FFT', 'FORWARD')` and `$:z2z_transform('IFT', 'INVERSE')`.
**The generated Fortran is identical to the old file inside both routines on
all three backends** -- which is also what makes the unbuildable HIP arm safe.
Elsewhere the generated file differs only by a removed no-op trailing `;` in
`FFT` and by wrapping the two `fftw_plan_many_dft` calls in `init_fft`, which
were 139 and 140 columns -- past the free-form limit, and fypp would have
folded them anyway.

Left written out, with the reason now a comment above `RFT`:

- **`RFT` and `HFT`.** Being inverses is *not* the reason -- `FFT` and `IFT`
  are inverses too and share a body fine. The reason is that `FFT`/`IFT` are
  one in-place complex-to-complex transform with a direction flag, while these
  two move between layouts: `RFT` takes the packed complex spectrum to the
  padded real field and `HFT` brings it back, so the arguments swap roles as
  well as types (`RFT(x, rx)` against `HFT(rx, x)`) and the calls are `Z2D`
  against `D2Z`. `RFT` also zeroes the padding columns on HIP, with no
  counterpart in `HFT`. A shared body would have to parameterise which dummy
  is the input, plus a block that exists on one side of one backend.
- **The three `free_fft` bodies.** The CUDA and HIP ones are the same eight
  lines up to the handle prefix, but the FFTW one is a different routine
  (a `target exit data`, four `fftw_destroy_plan`, four deallocates, three
  nullifies). Unifying the two device arms means a macro used twice to replace
  eight lines, in the arm nothing builds. Not worth it.
- Noticed while reading, **not fixed**: the FFTW `free_fft` deallocates
  `products` but does not `nullify` it (the `nullify` lists only `VVdz`,
  `VVdx`, `rVVdx`), so a second call would test `associated(products)` on a
  dangling pointer. Nothing calls it twice today. Flagging rather than fixing,
  since it is a behaviour change in a path no test covers.

## Settled recently: the pipeline instrumentation

`ys_solve_pipelined_lu_distributed` used to carry `roctxPush` / `tic` / step /
`toc` / `roctxPop` at twenty sites -- four lines of bookkeeping around one
line of work, with the schedule buried in its own telemetry. Two fypp macros
in the file's preamble, `timed_step(label, bucket)` and
`timed_comm_post(buffer, label=)`, now take the body as a `#:call` /
`#:endcall` block, so it stays ordinary Fortran and the label is written once
and expanded into both the push and the pop -- a mistyped pop label was
silently possible before.

Left in long form on purpose: the `MPI_Test` probe is untimed (the point of
the test is that it may cost nothing), and the outer roctx ranges wrap loops
rather than steps. One deliberate change: the timing report now runs *after*
the halo range is popped, because it barriers across `MPI_COMM_Y` once per
rank and a profiler was charging that to the halo.

**The generated Fortran is byte-identical apart from that one moved
`roctxPop`** -- which is how it was verified, and is worth reaching for
whenever a refactor can be shaped to allow it. `CHANNEL_Y_PIPELINE_TIMING=1`
was also run on CPU and GPU, since no test sets it.

## Settled recently: `-Dbodyforce` is gone

It allocated an `ny x nz x nx x 3` array `F`, mapped it to the device,
computed a wall closure on it and added it to the Orr-Sommerfeld and Squire
right-hand sides -- but nothing ever put a force *into* `F`: no deck
parameter, no file read, no initialiser. It was also broken: `buildrhs_prepare`
indexed `F` over the wrong y range on a y-split rank, so a `-Dbodyforce` build
died at np=2 on an out-of-bounds access while np=1 exited 0. Nothing built it,
which is why that went unseen.

Removing it changed no compiled code (`cpp` over the generated sources before
and after differs only in two comment lines of `header.h`). **If forcing is
wanted back, the missing half is the input path** -- re-adding the terms alone
would recreate exactly what was deleted.

## Settled earlier: the MPI/NCCL pipeline pair

The question was whether `ys_solve_pipelined_lu_distributed` and
`ys_solve_pipelined_lu_nccl_distributed` should be one routine. **They should
not**, and neither should the halo pair. Do not reopen this without new
information.

- **They differ in schedule, not only in transport.** The MPI version
  pre-posts every backward `Irecv` before the forward sweep and then drains
  backward batches *inside* the forward loop via `MPI_Test`, so batch k's back
  substitution overlaps batch k+1's factorisation. The NCCL version runs a
  complete forward loop and then a separate reverse loop. Unifying the two
  bodies either serialises the MPI path or invents an overlap the NCCL path
  does not have -- either way a behaviour change, not a readability change.
- **They are not two peers.** The MPI one is the only entry point; it
  dispatches to the NCCL one near the top and returns.
- **What genuinely was duplicated has been fixed.** Six complex values per
  line carry a batch's factor state and two carry a solve state. Those two
  numbers set the message counts *and* the per-batch strides, and they
  appeared at seven sites across the two routines -- if the transports had
  ever drifted apart the corruption would have been silent. They now appear
  only in `YS_PIPELINE_FACTOR_ELEMS`/`SOLVE_ELEMS`, read by
  `ys_ensure_pipeline_buffers` (sizing) and `ys_pipeline_batch_slice`
  (locating one batch). Each per-batch site went from three-to-five lines of
  arithmetic to one call.
- **The halo pair gets the same verdict for the same reason.**
  `ys_exchange_pipelined_solution_halos` posts both neighbours' Irecv/Isend
  and does one `Waitall`; the `_nccl` one does sendrecv-and-unpack per
  neighbour. MPI overlaps the two neighbours, NCCL does not. What remains in
  common is three lines of `count`/`lower_offset`/`upper_offset` inside a
  40-line routine, both visible on one screen -- naming that would be noise,
  not clarity.

## Looked at and deliberately left alone

Do not redo these without new information; the reasoning is in the commits.

- **The two pipeline pairs above.**
- **`RFT`/`HFT` and the three `free_fft` bodies** in `ffts.fypp` (above).
- **The test-support block in `convvelo`** (`fill_convvelo_synthetic_state_
  for_test` and the two synthetic value functions, ~55 lines, public, shipping
  in the production module). Every shape tried was worse: the filler writes
  three private variables, so a test-side module means making them public, and
  its two `!$omp target update to(...)` calls should stay with the
  declarations (see [[gpu-build-constrains-module-structure]]). The least-bad
  version trades one test-only public routine for another. Worth doing only
  with a shape that shrinks the public surface.
- **`mpi_autotune`'s copy of the blocking transpose.** Converting it would put
  two extra roctx markers inside the loop the autotuner times, and that
  measurement feeds the decomposition choice.
- **The halo unpack pair in `y_line_solvers`** -- `ys_unpack_lower_solution_
  halo` and `ys_unpack_upper_solution_halo`, which differ by destination
  array, not by an index; the reason is a comment at the site. The *pack* pair
  was unified.
- **Reconstruct (step 3 of the wall closure)** stays written twice: the host
  spells it `a - sum(b)`, the batched one `a - b1 - b2 - b3`, and those do not
  round alike.

## Open items that are the user's call, not yours

- **`FFTW_PATIENT`** (`ffts.fypp`) makes CPU runs non-reproducible at ULP
  level. Decided: keep it, and keep patching a scratch tree to
  `FFTW_ESTIMATE` for bit-identity checks. Do not re-litigate.
- **`-Dhalfchannel` and `-DphiNeumann`** do take effect (both change the
  output), but nothing checks their numbers against a reference.
- **The autotuner tries NCCL even when ranks share one GPU.** Arguably it
  should skip NCCL when it sees a shared device -- but that is a behaviour
  change, so ask.
- **15 GPU tests fail for environmental reasons on both boxes** (below).
  Worth raising with whoever owns the hpcx install.

## How to work

Small, independently verified commits -- one logical change each, so a
regression bisects to one peel. The gate, none of which is optional:

1. CPU `ctest` (30 tests) -- necessary, never sufficient. Add explicit runs of
   whatever branch you touched; see "ask which branch your run took" below.
2. An NVHPC GPU build, and the GPU suite compared against **HEAD built in a
   scratch tree**, not against an expectation. The environment moves
   underneath both boxes, so establish the failing set at HEAD every time.
3. **Run the `channel` binary and check `$?`.** `ctest` never runs it -- every
   test drives its own `test_*` binary -- which is how it segfaulted in
   `MPI_Finalize` on every run with 30 green tests, and how the `-Dbodyforce`
   crash surfaced. There is no longer a `-Dbodyforce` leg; `-Dhalfchannel` and
   `-DphiNeumann` still exist and nothing builds them, so compile one of those
   if you touch the code they guard.
4. Bit-identity against the parent commit, both trees built with
   `FFTW_PATIENT` patched to `FFTW_ESTIMATE`. The working recipe: `git archive
   HEAD` into a scratch dir, `sed` the plan type, build, then run
   `mpirun -np 2 .../channel` in a directory holding
   `tests/convvelo/dns_runtime.in` as `dns.in` and
   `tests/data/start_field_scalar.out` as `Dati.cart.out`, and `cmp` the
   `convvelo.*.bin`, `Dati.cart*.out` and `Runtimedata*` outputs. That deck
   drives convvelo, pressure_output and the solver in one run. `run.log`
   always differs -- autotune timings and rank interleaving.

### The GPU environment, as of this session

- **Export `UCX_MEMTYPE_CACHE=n`.** Without it, GPU-aware MPI point-to-point
  segfaults (`Caught signal 11 (invalid permissions for mapped object)` inside
  `ucp_eager_bcopy`) and 17 of 30 tests fail. With it, the two pipelined-LU
  regression tests pass and the failing set drops to 15. This is the UCX CUDA
  memtype-cache bug, not a code defect.
- **The remaining 15 are the tests whose xz transpose goes over MPI with
  device pointers.** That crash survives `UCX_TLS=tcp`/`sm`,
  `--mca coll ^hcoll`, `--mca pml ucx`, and `--mca pml ob1 --mca btl
  self,vader`. `CHANNEL_COMM=nccl` does not rescue them either -- it swaps two
  failures for two others, since NCCL wants one GPU per rank and the 3-rank
  tests share two. The same 15 fail at HEAD on both boxes.
- **Run the GPU suite serially.** `ctest -j2` on two GPUs adds failures the
  serial run does not have.
- **Gate leg 3 on GPU: `channel` itself exits 139 at np=2** on cetus, and does
  so at HEAD too (verified this session by running both trees). Same
  environmental crash as the 15. At np=1 it exits 0, so run leg 3 at np=1 on
  GPU and at np=2 on CPU.
- **Both boxes share the NVHPC tree over NFS**, so they share its hpcx/UCX
  bug. The old "30/30 on istmcetus" note is stale.
- Fresh configure recipe (`module load toolkits/nvhpc/25.9` silently does
  nothing on the workstation):
  `NV=/opt/Nvidia/nvhpc/Linux_x86_64/25.9`,
  `PATH=$NV/compilers/bin:$NV/comm_libs/12.9/hpcx/latest/ompi/bin:$PATH`,
  `-DCMAKE_CUDA_ARCHITECTURES=86 -DCMAKE_Fortran_COMPILER=$NV/compilers/bin/nvfortran`.
- **`$HOME` is shared between the workstation and `istmcetus`**, so a scratch
  tree staged under `~` locally is already visible there -- no copying. `/tmp`
  is not shared.

### Verifying an NCCL change

No test sets `CHANNEL_COMM`, and NCCL cannot make a communicator with two
ranks on one device, so the local single-GPU box cannot exercise it at all.
On `istmcetus`, multi-rank runs over a plain non-interactive `ssh` **do** work
(2 and 3 ranks, verified this session; the old `opal_init` note was about 6
ranks, i.e. oversubscription). The recipe, from the repo root:

```
CHANNEL_COMM=nccl CHANNEL_NPY=2 CHANNEL_Y_SOLVER=pipelined_lu \
CHANNEL_Y_PIPELINE_BATCHES=3 mpiexec -n 2 build-cetus/test_regression \
tests/data/dns_test_npy2.in
```

It checks against the committed reference, so a PASS is real evidence. And
`CHANNEL_COMM=nccl` on a build without NCCL `error stop`s, so a passing run
also proves the NCCL branch was actually taken.

### Techniques worth reusing, in order of how much they save

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
  alltoall path is never entered. Likewise a 1-batch pipelined solve makes
  every per-batch offset collapse to its first term, so it proves nothing
  about the stride arithmetic -- this session's check swept
  `CHANNEL_Y_PIPELINE_BATCHES=1..5` at npy=3 for exactly that reason. Same
  lesson as `CHANNEL_OVERLAPPING=1` and `CHANNEL_YS_FORCE_CUSTOM_GPSV`, which
  no test sets and which was hiding a GPU crash: "30/30 green" is a statement
  about the 30.
- **Keep the bit-identity trees pristine.** A scratch tree you also used for a
  flag experiment is no longer the thing you meant to compare: toggling
  `header.h` in the work tree and rebuilding it silently turned a
  bit-identity run into a comparison of two different builds. One rebuild
  lost. Do side experiments in their own copy.
- **Check that a flag you are relying on actually took effect.** Two claims
  went wrong this way. `CHANNEL_OVERLAPPING=1` is bit-identical by design, so
  output cannot confirm it -- the memory line does (buffers double). And a
  comparison meant to isolate one flag had a second flag varying with it.
- **Compare against the checked-in references too**, not only base-vs-work.
  A reference comparison says "correct"; base-vs-work only says "unchanged".
- **Re-check the obstacle before working around it.** This handoff once said
  a shared transform module would need `VVdx`/`VVdz` behind an `#if`. They had
  not for some time -- one `grep` retired the whole complication. Notes about
  *why something is hard* age faster than notes about what was done. The
  `opal_init` and "30/30 on cetus" notes were both stale this session too.
- `fprettify` corrupts `.fypp` (it eats the space in `${macro}$ - x`); the
  pre-commit hook is restricted to `\.f90$`. `pre-commit` is not on PATH on
  this box -- run `fprettify -i 2 -w 2` from a throwaway venv instead.

Report honestly when something is reverted and why. Several extractions have
been abandoned across these sessions for good reasons, and recording that has
been more useful than forcing them through.
