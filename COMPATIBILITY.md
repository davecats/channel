# COMPATIBILITY.md — keeping one source tree working on NVIDIA and AMD

`channel` runs on two GPU platforms that accept different subsets of Fortran +
OpenMP offload:

| | NVIDIA (HoreKa A100 / H100) | AMD (Hunter MI300A) |
|---|---|---|
| compiler | `nvfortran` (NVHPC 25.x) | Cray `ftn` (CCE 20/21) or `amdflang` (TheRock) |
| libraries | cuFFT, cuSPARSE, NCCL (`cudafor`, `cufft` modules) | hipFFT, hipSPARSE, RCCL (hipfort) |
| memory | discrete GPU, 40–94 GB, separate host memory | APU, 128 GB HBM shared by CPU and GPU (`HSA_XNACK=1`) |
| `YS_COEF` (`build_options.h`) | **complex** | real |
| macro | `HAVE_CUDA` | `HAVE_HIP` (`_CRAYFTN` on CCE) |

Every item below broke one platform after a PR tested only on the other, and
was fixed afterwards (PR #25, PR #26). Each entry gives the symptom, the cause
and the form to write instead. **Read it before opening a PR, and run the
checklist at the end.**

## 1. No `merge()` with a complex value in device code (nvfortran)

- **Symptom:** the build "succeeds" (CMake exits 0) but the log contains
  `Accelerator restriction: unsupported operation: CDSELECT`; at run time the
  kernel is missing: `Failed to find device function nvkernel_…`, for every
  solver.
- **Cause:** nvfortran cannot offload `merge(x, y, c)` when `x`, `y` are
  complex. `YS_COEF` / `YS_CZERO` are complex on CUDA and real on HIP, so code
  that is fine on AMD becomes a complex merge on NVIDIA.
- **Write instead** (valid for real and complex `YS_COEF` alike):

  ```fortran
  a = YS_CZERO; if (i >= 2) a = ds(p)          ! scalar target
  if (i >= 2) then                              ! array element: no double store
    ys_batch_ds(p) = c_m2
  else
    ys_batch_ds(p) = (0.0d0, 0.0d0)
  end if
  ```

  `cmplx(merge(real(x), …), merge(aimag(x), …))` also offloads, but breaks the
  real `YS_COEF` of HIP — do not use it. Fixed in 7243565 (PR #25, 20 sites)
  and in the factor-ahead band macro of PR #26.

## 2. `map(...)` is not a data-sharing clause (nvfortran **and** CCE)

- **Symptom:** nvfortran `NVFORTRAN-S-0155-loc must appear in a SHARED or
  PRIVATE clause`; CCE `ftn-1510 … Variable "LOC" must have its data scope
  explicitly declared because DEFAULT(NONE) was specified`. amdflang accepts it.
- **Cause:** on a combined `target teams distribute parallel do default(none)`
  a variable named only in `map(from: x)` still needs a data-sharing attribute.
- **Write instead:** `map(from: loc) shared(loc, …)`. (PR #26,
  `ys_dvdy_build_tables`: the PR head did not compile with CCE 21 either.)

## 3. nvfortran leaks a module's private imports (`default(none)` false errors)

- **Symptom:** `NVFORTRAN-S-0155-<name> must appear in a SHARED or PRIVATE
  clause` for a name that **is** in the clause, often in a file the PR did not
  touch (PR #26: `der` in `compact_component_solve`, `k2`, `ni`, `nx0` in
  `channel_equations`). gfortran, CCE and amdflang accept the code.
- **Cause:** module `A` imports `x` at module level (`use channel_grid, only:
  …, k2` in `y_line_solvers`, private). A scope that uses `A` — even with
  `use A, only: something_else` — and also sees `x` through its own import
  then has two paths to `x`, and nvfortran's `default(none)` check rejects the
  kernel. Importing the same module both at module level and again inside a
  routine (`use channel_grid, only: k2` in the routine, `only: ny0, …` at the
  top) triggers it too.
- **Write instead:**
  - In a module that others `use` (`y_line_solvers`, `y_schur_solver`, …),
    import variables that kernels elsewhere also name (`der`, `k2`, `ni`,
    grid sizes) **inside the routines that need them**, not at module level —
    unless the module already imports that module at module level, in which
    case add the name to that existing list.
  - In a module full of kernels (`channel_equations`), reach such a module
    through a small kernel-free routine with its own local `use`
    (`register_matrix_free` in `channel_equations`, from the PR #26 nvfortran
    fixes), never with a module-level `use`.
- Related, older rules (Aug 2026): a module must not both `use` another module
  and declare its own `!$omp declare target` variables; device-callable
  procedures cannot cross a module boundary (`nvlink` undefined references;
  `-gpu=rdc` does not help).

## 4. Never forward an assumed-shape dummy into a `declare target` call (nvfortran)

- **Symptom:** `CUDA_ERROR_ILLEGAL_ADDRESS` at launch; `compute-sanitizer`
  shows no invalid access.
- **Write instead:** name the module arrays at the call site; when one kernel
  must serve two array sets, generate it per set with a fypp macro
  (`interleaved_penta_solver`), not through a dummy argument.

## 5. The APU hides mapping bugs; discrete GPUs do not

- On MI300A with `HSA_XNACK=1` host and device share one memory: a missing
  `map`, a host array read inside a kernel, or an unbalanced
  `target enter/exit data` still gives the right answer at full speed. On
  NVIDIA the same code crashes, or silently frees and re-uploads data.
- Example: an unmatched `target exit data map(release:)` of the Schur views
  (2bcdc5a) dropped the refcount of the whole byte workspace to zero; every
  y-solve then re-uploaded 2.5 GB+ — 2× step time, +6 to +25 GB per GPU on
  HoreKa, invisible on Hunter. Fixed in 6b08bfa.
- **Rule:** every `target exit data map(release: v)` must be paired with a
  `target enter data` of `v` on the same code path; keep the guard that
  decides both in one place. Prefer views into the byte workspace that are
  never entered or released on their own.

## 6. Memory is not free on NVIDIA

- An A100 has 40 GB, an H100 94 GB, against 128 GB per MI300A APU. A default
  that costs memory for a few percent of speed on AMD can leave no headroom on
  NVIDIA: stacking the y-solves took the 2-node case from 15 to 37 GB per A100
  for no gain. CUDA therefore defaults to `CHANNEL_Y_STACK_MAX=1`
  (`ys_stack_limit` in `y_line_solvers`); keep platform-specific defaults
  under `#if defined(HAVE_CUDA)` rather than removing them. (With PR #26 the
  stack became cheap — +0.05 to +1 GB per H100 for STACK_MAX 4 / 8, 2–3 %
  faster — so the default is worth re-measuring, not deleting.)
- `CHANNEL_MPI_AUTOTUNE` trials hold their candidates' workspaces at once:
  peak GPU memory during autotune was 28 GB per H100 for a run whose
  estimate and steady state are 12–16 GB. Check autotune against the memory
  of the smallest NVIDIA card, not the APU.
- Any new workspace must be counted in `get_solver_workspace_estimate` /
  `ys_workspace_layout`, which size the reservation at startup.

## 7. Every backend call needs both spellings

- A feature written only for hipFFT / hipfort / RCCL does not exist on
  NVIDIA. `CHANNEL_SINGLE` arrived HIP-only in PR #26 (`#error` on CUDA); in
  the CUDA port that goes with it the cuFFT calls follow the FFT kind through `CUFFT_ZEXEC`, `CUFFT_ZTYPE`, `CUFFT_C2REXEC`,
  `CUFFT_R2CEXEC`, `CUFFT_C2RTYPE`, `CUFFT_R2CTYPE` in `ffts.fypp`, as the
  hipFFT calls do with their `HIPFFT_*` macros. Add both, or make the CMake
  option fail on the other backend.
- NCCL and RCCL differ in what blocks: a receive posted early sits in its
  lane's stream ahead of sends the partner waits for. The y pipeline keeps
  NCCL receives one batch ahead (f43456e) and the meet swap blocking under
  NCCL (edb06bf). Test both `CHANNEL_COMM=mpi` and `nccl`.
- `dreal` / `dimag` take double complex only; use the generic `real(x, kind)` /
  `aimag(x)` so the code works for `RHS_COMPLEX` in single precision.
- hipfort is not versioned with the code: `HIPFFT_INVERSE` exists in the
  hipfort `rocm-7.0.2` tag and is gone on hipfort master. Build hipfort from
  the tag that matches the loaded ROCm.

## 8. CPU (gfortran) build

- `c_loc(x)` needs `x` to have `TARGET` (or `POINTER`); nvfortran does not
  check this, gfortran rejects it (c39e144, `mpi_transpose`).
- gfortran at `-O2` miscompiles 64-bit wraparound arithmetic; never rely on
  integer overflow.

## 9. Branches compiled by one compiler only

- `#if defined(_CRAYFTN)` (CCE), `HAVE_HIP`, `HAVE_CUDA` branches are compiled
  by one toolchain only. A change inside one of them, or a new kernel that
  is "off on CCE", needs a build with that toolchain. Several PR #26 commits
  say "the CCE branch is edited alike but was not compiled".

## 10. Build traps

- **CMake exits 0 even when nvfortran drops a kernel.** Grep the build log:
  `grep "Accelerator restriction" build.log` must be empty.
- Building on a CPU-only node, nvfortran targets a list of GPUs or the build
  host's; pass the target explicitly: `-DCMAKE_Fortran_FLAGS="-gpu=cc80"`
  (A100) / `-gpu=cc90` (H100). `CMAKE_CUDA_ARCHITECTURES` does not reach the
  Fortran flags.
- `make` stops at the first failing file: build with `-k` (`cmake --build …
  -- -k`) to see every error of a port at once.
- CCE compiles slowly (the three trees of PR #26 took ~50 min on Hunter with
  32 build jobs); start that build first.

## 11. A branch in front of a collective must be the same on every rank

Not a vendor difference, but it hides the same way: a GPU benchmark at
production size never takes the branch that breaks.

- **Symptom:** a run with an uneven y split hangs at 100 % CPU; one rank sits
  in `MPI_Allgather`, the others in `MPI_Wait` of the pipeline. PR #26:
  `regression_test_pipelined_lu_npy3_3rank` and
  `regression_test_scalar_pipelined_lu_npy3` never finish.
- **Cause:** `ys_dvdy_register` chose the dv/dy table solve when *this rank*
  had at least 6 rows; the table is built collectively on `MPI_COMM_Y`. Once
  the wall ranks got one row fewer (2ce7ff1), ny = 16 on 3 y ranks splits 4,
  6, 5 and only the middle rank enters the collective. It hangs on CPU and
  on H100 (same deck, 3 ranks); Hunter's CCE build does not hang only because
  the table solve is compiled out there (`#if !defined(_CRAYFTN)`, see §9).
- **Write instead:** decide from data every rank has alike — here the fewest
  rows of any y rank, computed from `y_split_first` — or reduce the decision
  (`MPI_Allreduce` of the flag) once at setup, never per rank and per call.
  Any change to the y split must be followed by every code that assumed the
  old one (`mpi_autotune`'s `valid_decomposition` still checked the even
  split).

## Checklist for every GPU PR

1. NVIDIA: build with nvfortran (`-gpu=cc80` or `cc90`), `grep "Accelerator
   restriction"` on the log must be empty, and run at least one step on 1 and
   on 2 nodes, `CHANNEL_COMM=nccl` and `mpi`.
2. AMD: build with **both** CCE (`tools/build_hunter_cray.sh`) and amdflang
   (`tools/build_hunter_amd.sh`) when a `_CRAYFTN` branch or a kernel with
   `default(none)` changed; run on Hunter.
3. CPU: gfortran build and the **full** `ctest` (with `--timeout 300`: the
   odd-rank tests `*_npy3*` are the ones that catch uneven-split bugs and
   deadlocks; a hung test otherwise runs forever).
4. Results: bit-identical to the previous commit where the PR claims so, on
   both platforms; peak GPU memory against the startup estimate on NVIDIA.
5. New code: no complex `merge` in kernels, `map` variables also in
   `shared()`, no module-level imports of kernel variables in modules others
   use, every `exit data` paired with an `enter data`, no rank-local test
   in front of a collective.
