
# Module environment for building and running channel on HLRS Hunter
# (TheRock/amdflang toolchain). Source it, do not execute it:
#
#     source tools/hunter_env.sh
#
# Used by tools/build_hunter.sh; source it as well in interactive shells and
# job scripts before running ctest / mpiexec, so the runtime libraries match
# the ones the binaries were linked against.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "hunter_env.sh must be sourced: source ${BASH_SOURCE[0]}" >&2
    exit 1
fi

hunter_therock_version="afar-24.3.0-multiarch-10.1.0-592954c"
hunter_module_root="${THEROCK_MODULE_ROOT:-/opt/hlrs/testing/unsupported/modulefiles/therock-${hunter_therock_version}}"

if ! type module >/dev/null 2>&1; then
    echo "module command not available; use a Hunter login shell." >&2
    return 1
fi

module purge &&
module use "$hunter_module_root" &&
module load "therock/${hunter_therock_version}" &&
module load mpich/4.3.0b1 &&
module load rocm-libs/4.0.0 &&
# The mpich module is a Fortran shim over Cray MPICH: runs need the PALS
# launcher and libfabric, which `module purge` removed.
module load cray-pals libfabric || {
    echo "hunter_env.sh: loading modules failed." >&2
    return 1
}


# Cray MPICH must know buffers may be device memory (GPU-aware MPI).
export MPICH_GPU_SUPPORT_ENABLED=1
# Required for correct results; see debugging notes (GPU touches host memory).
export HSA_XNACK=1
