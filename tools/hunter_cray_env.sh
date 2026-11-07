# Module environment for building and running channel on HLRS Hunter
# with the HPE Cray compiler (CCE / crayftn).  Source it, do not execute it:
#
#     source tools/hunter_cray_env.sh
#
# Used by tools/build_hunter_cray.sh; source it as well in interactive shells
# and job scripts before running ctest / mpiexec, so the runtime libraries
# match the ones the binaries were linked against.
#
# Versions follow the HPE CPE 26.03 release (HLRS/APU/testing-2026.2 stack):
# CCE 21 supports AMD GPU offload only with ROCm >= 7.0, so rocm/7.0.2.
# The HLRS stack module itself is not loaded: it also sets HSA_ENABLE_SDMA=0,
# GPU_*_ALLOC_PERCENT, NCCL_* ... which tools/hunter_env.sh (amdflang) does not,
# and the two builds should run in comparable environments.
# Override individual versions via HUNTER_CCE_VERSION etc. if needed.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "hunter_cray_env.sh must be sourced: source ${BASH_SOURCE[0]}" >&2
    exit 1
fi
hunter_prgenv_version="${HUNTER_PRGENV_VERSION:-8.7.0}"
hunter_cce_version="${HUNTER_CCE_VERSION:-21.0.2}"
hunter_craype_version="${HUNTER_CRAYPE_VERSION:-2.7.36}"
hunter_rocm_version="${HUNTER_ROCM_VERSION:-7.0.2}"
hunter_mpich_version="${HUNTER_MPICH_VERSION:-9.1.0}"
hunter_pals_version="${HUNTER_PALS_VERSION:-1.7.1}"
if ! type module >/dev/null 2>&1; then
    echo "module command not available; use a Hunter login shell." >&2
    return 1
fi
module purge &&
module load "PrgEnv-cray/${hunter_prgenv_version}" &&
module load "cce/${hunter_cce_version}" "craype/${hunter_craype_version}" &&
module load craype-x86-genoa craype-network-ofi xpmem &&
module load "rocm/${hunter_rocm_version}" craype-accel-amd-gfx942 &&
module load "cray-mpich/${hunter_mpich_version}" "cray-pals/${hunter_pals_version}" || {
    echo "hunter_cray_env.sh: loading modules failed." >&2
    return 1
}
# Without this the binaries pick up the *default* CPE runtime libraries from
# /opt/cray/pe/lib64 instead of the CCE / MPICH versions loaded above.
export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
# Cray MPICH must know buffers may be device memory (GPU-aware MPI); without
# it, passing device buffers to MPI SIGBUSes.  Also makes the compiler
# wrappers link libmpi_gtl_hsa.
export MPICH_GPU_SUPPORT_ENABLED=1
# Kept for a like-for-like comparison with the amdflang build, where results
# were wrong without it.  Override with HUNTER_HSA_XNACK=0 before sourcing,
# or simply re-export HSA_XNACK afterwards.
export HSA_XNACK="${HUNTER_HSA_XNACK:-1}"
# ROCm 7.0.2's RCCL segfaults in ncclIbInit (ibv_free_device_list) on these
# nodes, which have libibverbs but no InfiniBand: the first RCCL communicator
# of every multi-rank run dies.  Intra-node RCCL does not need IB anyway.
# See tools/cray_bugs/rccl-ib-init-segv.
export NCCL_IB_DISABLE="${NCCL_IB_DISABLE:-1}"
