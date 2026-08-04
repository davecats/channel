# Continuing on the HPC cluster: multi-GPU validation and performance

Copy the block at the end as the opening message of the cluster session.

Branch: **`main_GPU`** (`git@github.com:davecats/channel.git`), at `385552e`.
Read this file, then `src/physics/README.md`, then `NEXT_SESSION.md` (whose
"Do this first" is a *different*, still-open subject — the pack kernels).

---

## Why this session exists

**Multi-GPU has never been validated.** Every multi-rank GPU test fails on both
development boxes for environmental reasons, so the entire distributed GPU path
— the xz transpose with device pointers, the y-Schur solver across ranks, the
pipelined-LU pipeline, NCCL — has only ever been exercised at **np=1** here, or
on the CPU. The code is written and self-consistent; it is not *verified*
distributed on a GPU.

That is the first job. The second is performance, which cannot even be
discussed until the first is settled.

### What "environmental" means, precisely

On `ISTM-io2` / `istmcetus` (shared NVHPC tree over NFS, so shared UCX bug),
with `UCX_MEMTYPE_CACHE=n` exported and run serially, `ctest` gives
**16 failures out of 31**, always these:

```
 2 regression_test_4rank                    17 regression_test_scalar_4rank
 3 regression_test_autotune_default_2rank   18 regression_test_scalar_6rank_np4
 4 regression_test_autotune_report_2rank    19 regression_test_scalar_6rank_np4_schur
 5 regression_test_manual_decomp_2rank      20 regression_test_scalar_pipelined_lu_npy3
 9 pressure_dpdy_2rank                      25 reduced_ghost_backend_2rank
10 pressure_dpdy_npy2_2rank                 26 reduced_ghost_backend_2rank_schur
11 pressure_dpdy_npy2_2rank_schur           27 reduced_ghost_backend_2rank_pipelined_lu
12 pressure_dpdy_npy3_3rank_schur           31 seeded_start_field
```

Every one is multi-rank. The signature is a segfault inside UCX, not in our
code:

```
Caught signal 11 (Segmentation fault: invalid permissions for mapped object)
... ucp_eager_bcopy ...
```

It survives `UCX_TLS=tcp/sm`, `--mca coll ^hcoll`, `--mca pml ucx`,
`--mca pml ob1 --mca btl self,vader`, and `CHANNEL_COMM=nccl`. Without
`UCX_MEMTYPE_CACHE=n` it is 17–18 failures instead of 16.

**On a healthy cluster these should pass.** If they do, that confirms the
diagnosis and unblocks everything. If some still fail, those are real bugs and
this is the first time anyone has seen them — treat that as the finding, not as
a nuisance.

## Step 1 — build

Match `CMAKE_CUDA_ARCHITECTURES` to the cluster's GPUs (`90` for H100, `80` for
A100, `100`/`120` for newer; the build defaults to a wide list if you pass
something below 75).

```bash
cmake -S . -B build-gpu \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DCMAKE_C_COMPILER=nvc \
  -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build-gpu -j
```

Check the configure output for three lines that decide what you can test:

- `NCCL y-pipeline backend enabled: version X, native AlltoAll=ON/OFF` — if
  NCCL is not found, `CHANNEL_COMM=nccl` silently does nothing. Set
  `NCCL_ROOT_DIR` if needed.
- `NVTX profiling ranges enabled` — without it every `roctxPush`/`Pop` compiles
  out and Nsight will show no ranges. Needed for step 4.
- `Physics option:` lines — should be absent; you want the default build.

Also build the CPU reference (`-DCMAKE_Fortran_COMPILER=mpifort`) on the same
machine. Several checks below compare against it.

## Step 2 — establish the failing set, then validate multi-GPU

```bash
export UCX_MEMTYPE_CACHE=n          # harmless if the cluster does not need it
ctest --test-dir build-gpu -j1      # serially: the tests share a work directory
```

Record the failing set **before reading anything into any other run**. The
suite is 31 tests on a GPU build (32 on CPU; `physics_options_compile` is
GNU-only).

Then, in order:

1. **`seeded_start_field` is the cheapest real multi-rank check.** It runs the
   `channel` binary at np=1 and np=2 and at `CHANNEL_NPY=2` and demands
   byte-identical output. It is new, it is the only test that runs the
   production binary multi-rank, and it fails here purely environmentally.
2. **The three y-solver paths, each multi-rank.** They are independent code:
   ```bash
   CHANNEL_Y_SOLVER=schur        CHANNEL_NPY=4 mpirun -np 4 ./channel
   CHANNEL_Y_SOLVER=pipelined_lu CHANNEL_NPY=4 mpirun -np 4 ./channel
   CHANNEL_NPXZ=4 CHANNEL_NPY=1  mpirun -np 4 ./channel
   ```
3. **Both transports**, since NCCL and MPI differ in *schedule*, not just
   transport: `CHANNEL_COMM=mpi` and `CHANNEL_COMM=nccl`. NCCL needs **one GPU
   per rank** — it cannot make a communicator with two ranks on one device, and
   the autotuner will still try it (a known wart; ask before changing that).
4. **Batch sweep for the pipelined LU**, which needs at least two batches and
   `npy >= 2`: `CHANNEL_Y_PIPELINE_BATCHES=1..5` at `CHANNEL_NPY=3`. Note
   `y_batches=` in the decomposition line is the *autotuner's* field and reads
   0 even when the env var is set — confirm the count took with
   `CHANNEL_Y_PIPELINE_TIMING=1`, whose report prints `nbatches` per rank.

### The check that makes this easy — use it

`[velocity] seed = <n>` (new) makes a restartless run's start field a **pure
function of the deck**: bit-for-bit identical at any rank count and any
`CHANNEL_NPY`. So you can validate a decomposition without shipping a restart
file:

```ini
[velocity]
seed = 4242
```

Run the same deck at npy=1, 2, 4 and at different npxz and compare
`Dati.cart.out`. With `nstep = 0` the time loop never runs and the file *is*
the start field, so any difference is in the decomposition, not the solver. Then
raise `nstep` and the same comparison tests the solver instead.

Caveat: across *compilers* the field agrees to ~4 ULP, not bit-for-bit — the
mixer is bit-identical but `EXP(dcmplx(0,phase))` and the `tanh` mesh come from
each compiler's libm. So compare GPU-to-GPU and CPU-to-CPU exactly, and
GPU-to-CPU with a tolerance.

Beyond that, `tests/compare_complex_fields.py` compares fields, and
`tests/data/` holds committed reference `end_field*.out` for the small decks.

**`ctest` green is necessary and not sufficient.** Run the `channel` binary
yourself and check `$?` — that is how it once segfaulted in `MPI_Finalize` on
every run with a fully green suite.

## Step 3 — decomposition rules and the deck to scale up

Manual decomposition must satisfy

```
nproc = CHANNEL_NPXZ * CHANNEL_NPY
CHANNEL_NPXZ divides nx + 1
CHANNEL_NPXZ divides nzd            (nzd = 3*nz, adjusted by the FFT-fit logic)
```

The bundled decks are tiny — the test decks are `nx=5 ny=16 nz=8`, which will
tell you nothing about performance. `examples/channel_kmm1987/dns.in`
(`nx=95 ny=128 nz=80`, Re_tau ~ 180) is the smallest realistic case and the
right starting point; scale it up from there for strong/weak scaling. Choose
`nx+1` and `nzd` with plenty of divisors — `nx=95` gives `nx+1 = 96` and
`nzd = 240`, which is deliberately friendly.

For benchmarking:

```bash
export CHANNEL_DISABLE_RESTART_WRITE=1   # skip the final write
```

and set `nstep` to something small and fixed rather than using `t_max`.

## Step 4 — performance

`#define chron` is on by default, so every step prints
`STEP n/N TIME PER TIMESTEP x ELAPSED RUN TIME y`. **Discard the first few
steps**: FFT plan creation, the autotuner and first-touch allocation all land
there.

Knobs, roughly in order of expected effect:

| Knob | What to try |
| --- | --- |
| `CHANNEL_NPY` / `CHANNEL_NPXZ` | the main decomposition trade-off at fixed `nproc` |
| `CHANNEL_MPI_AUTOTUNE` | default on; `report` prints a recommendation without applying it, `off` disables. Compare its pick against your own sweep |
| `CHANNEL_COMM` | `mpi` vs `nccl` — expect this to matter most at large npy |
| `CHANNEL_Y_SOLVER` | `schur` vs `pipelined_lu`; they scale differently in npy |
| `CHANNEL_Y_PIPELINE_BATCHES` | more batches = more overlap, more messages |
| `CHANNEL_Y_SCHUR_PASSES` | e.g. `2x4`; product must equal `CHANNEL_NPY` |
| `CHANNEL_Y_SCHUR_EXCHANGE` | `auto` / `alltoall` / `allgather` |
| `CHANNEL_OVERLAPPING` | overlap comms with computation. **Bit-identical by design**, so output cannot confirm it took — the memory report does (buffers double) |
| `CHANNEL_YS_BATCH_MAX_COMPLEX`, `CHANNEL_YS_CHUNK_NX` | y-solve workspace and chunking caps |
| `CHANNEL_YS_FORCE_CUSTOM_GPSV` | built-in pentadiagonal solve instead of the vendor batched one |

`CHANNEL_EXIT_AFTER_MPI_AUTOTUNE=1` exits after parsing and decomposition —
useful to validate a deck and see the autotuner's choice cheaply.

Two instrumented views:

- `CHANNEL_Y_PIPELINE_TIMING=1` — per-rank breakdown of the pipelined LU,
  including `nbatches`. Note the report barriers across `MPI_COMM_Y` once per
  rank, so it costs something; it is deliberately outside the halo roctx range.
- **Nsight Systems** — the code is covered in NVTX ranges (`timestep`,
  `rk_substep`, `transform_to_physical`, `buildrhs_prepare`, `linsolve_*`,
  `outstats`, plus fine-grained ones in the solver), but only if the configure
  line said `NVTX profiling ranges enabled`.

```bash
nsys profile -t cuda,nvtx,mpi -o prof_%q{OMPI_COMM_WORLD_RANK} ./channel
```

Also print-worthy: startup prints an estimated per-rank memory breakdown before
allocating, which is the quickest way to see whether a case will fit.

## Traps this project has already paid for

- **Establish the failing set at HEAD every time.** The environment moves under
  these builds; a note saying "30/30 on cetus" was stale twice.
- **`module load toolkits/nvhpc/25.9` silently did nothing** on the dev boxes;
  the existing build worked only because its CMake cache held an absolute
  compiler path. On the cluster, check `which nvfortran` actually changed.
- **A flag you rely on may not have taken effect.** Check, do not assume — the
  `y_batches=0` case above is exactly this.
- **`gfortran -O2` miscompiles 64-bit wraparound multiply** (found this week
  writing the seed mixer: correct at `-O0`, constant output at `-O2`). If you
  write any hash or PRNG here, test it at `-O2`, not at the project's `-O0`.
- **CPU runs are not bit-reproducible** unless `FFTW_PATIENT` is patched to
  `FFTW_ESTIMATE` in `src/fft/ffts.fypp`. Irrelevant on the GPU path, relevant
  the moment you compare against the CPU build.

## Open questions worth answering while you are there

- Do the 16 environmental failures pass on real hardware? (The whole premise.)
- Should the autotuner skip NCCL when it sees ranks sharing a device? It
  currently tries and fails. That is a behaviour change — **ask before doing it**.
- `-DCHANNEL_HALF_CHANNEL=ON` and `-DCHANNEL_PHI_NEUMANN=ON` compile and run,
  but **nothing checks their numbers**. A cluster with a known validated case
  is the place to fix that.
- Does the autotuner's decomposition pick actually win on this machine?

---

## Opening message for the cluster session

> Continue `channel` on the HPC cluster, on branch `main_GPU`. Read
> `HPC_SESSION.md` first, then `src/physics/README.md`, and the memory notes
> `channel-src-layout` and `channel-verification-gates`.
>
> Goal, in order: (1) build with NVHPC for this cluster's GPU architecture and
> establish the `ctest` failing set at HEAD; (2) find out whether the 16
> multi-rank failures listed in `HPC_SESSION.md` are the known UCX environment
> bug or real — this is the first time the distributed GPU path will run on
> healthy hardware, and nothing beyond np=1 has ever been validated; (3) once
> multi-GPU is trusted, scale `examples/channel_kmm1987/dns.in` up and measure
> strong and weak scaling across `CHANNEL_NPY`/`CHANNEL_NPXZ`,
> `CHANNEL_COMM=mpi` vs `nccl`, and the two y solvers.
>
> Use `[velocity] seed` with `nstep = 0` to compare decompositions without a
> restart file. Report what the autotuner picks against what actually wins.
> Do not treat a failure as environmental without reproducing the UCX signature
> quoted in `HPC_SESSION.md`.
