cmake_minimum_required(VERSION 3.16)
project(ChannelSolver LANGUAGES Fortran C)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

# --------------------------------------------------------------------
# CMake module path (for custom FindFFTW)
# --------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

function(read_nccl_header_version INCLUDE_DIR OUT_VAR)
    set(_NCCL_HEADER "")
    if(EXISTS "${INCLUDE_DIR}/nccl.h")
        set(_NCCL_HEADER "${INCLUDE_DIR}/nccl.h")
    elseif(EXISTS "${INCLUDE_DIR}/rccl/rccl.h")
        set(_NCCL_HEADER "${INCLUDE_DIR}/rccl/rccl.h")
    endif()

    set(_NCCL_VERSION "unknown")
    if(_NCCL_HEADER)
        set(_NCCL_MAJOR "")
        set(_NCCL_MINOR "")
        set(_NCCL_PATCH "")
        file(STRINGS "${_NCCL_HEADER}" _NCCL_MAJOR_LINE REGEX "^#define[ \t]+NCCL_MAJOR[ \t]+[0-9]+")
        file(STRINGS "${_NCCL_HEADER}" _NCCL_MINOR_LINE REGEX "^#define[ \t]+NCCL_MINOR[ \t]+[0-9]+")
        file(STRINGS "${_NCCL_HEADER}" _NCCL_PATCH_LINE REGEX "^#define[ \t]+NCCL_PATCH[ \t]+[0-9]+")
        if(_NCCL_MAJOR_LINE MATCHES "^#define[ \t]+NCCL_MAJOR[ \t]+([0-9]+)")
            set(_NCCL_MAJOR "${CMAKE_MATCH_1}")
        endif()
        if(_NCCL_MINOR_LINE MATCHES "^#define[ \t]+NCCL_MINOR[ \t]+([0-9]+)")
            set(_NCCL_MINOR "${CMAKE_MATCH_1}")
        endif()
        if(_NCCL_PATCH_LINE MATCHES "^#define[ \t]+NCCL_PATCH[ \t]+([0-9]+)")
            set(_NCCL_PATCH "${CMAKE_MATCH_1}")
        endif()
        if(_NCCL_MAJOR AND _NCCL_MINOR AND _NCCL_PATCH)
            set(_NCCL_VERSION "${_NCCL_MAJOR}.${_NCCL_MINOR}.${_NCCL_PATCH}")
        endif()
    endif()

    set(${OUT_VAR} "${_NCCL_VERSION}" PARENT_SCOPE)
endfunction()

# --------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------
add_library(channel_core
    config.f90
    rbmat.f90
    roctx.f90
    byte_workspace.f90
    y_pipeline_nccl.f90
    y_schur_solver.f90
    y_line_solvers.f90
    compact_line_solvers.f90
    mpi_transpose.f90
    ffts.f90
    mpi_autotune.f90
    dnsdata.f90
    convvelo.f90
    pressure_output.f90
    driver.f90
)

# --------------------------------------------------------------------
# MPI
# --------------------------------------------------------------------
option(NOMPI "Build without MPI support (no reading / writing of files)" OFF)

if(NOT NOMPI)
    find_package(MPI REQUIRED Fortran)
    target_compile_definitions(channel_core PUBLIC HAVE_MPI)
    target_link_libraries(channel_core PUBLIC MPI::MPI_Fortran)
endif()

# --------------------------------------------------------------------
# Executable
# --------------------------------------------------------------------
add_executable(channel ${SOURCES} channel.f90)
target_link_libraries(channel PRIVATE channel_core)

add_executable(post_pressure post_pressure.f90)
target_link_libraries(post_pressure PRIVATE channel_core)

add_executable(post_convvelo post_convvelo.f90)
target_link_libraries(post_convvelo PRIVATE channel_core)

# --------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------

set(TEST_NPROCS 1 CACHE STRING "Number of MPI ranks for regression test")

add_executable(test_regression test_regression.f90)
target_link_libraries(test_regression PRIVATE channel_core)

add_executable(test_regression_scalar test_regression_scalar.f90)
target_link_libraries(test_regression_scalar PRIVATE channel_core)

add_executable(test_pressure_dpdy test_pressure_dpdy.f90)
target_link_libraries(test_pressure_dpdy PRIVATE channel_core)

add_executable(test_convvelo_stats test_convvelo_utils.f90 test_convvelo_stats.f90)
target_link_libraries(test_convvelo_stats PRIVATE channel_core)

add_executable(test_convvelo_runtime test_convvelo_utils.f90 test_convvelo_runtime.f90)
target_link_libraries(test_convvelo_runtime PRIVATE channel_core)

add_executable(test_post_convvelo test_convvelo_utils.f90 test_post_convvelo.f90)
target_link_libraries(test_post_convvelo PRIVATE channel_core)

add_executable(test_y_pencil_transpose test_y_pencil_transpose.f90)
target_link_libraries(test_y_pencil_transpose PRIVATE channel_core)

add_executable(test_mean_correction test_mean_correction.f90)
target_link_libraries(test_mean_correction PRIVATE channel_core)

add_executable(test_mpi_autotune_candidates test_mpi_autotune_candidates.f90)
target_link_libraries(test_mpi_autotune_candidates PRIVATE channel_core)

# Add a CTest target
enable_testing()
include(cmake/tests.cmake)

# --------------------------------------------------------------------
# Compiler flags
# --------------------------------------------------------------------



if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O0 -malign-double -fall-intrinsics -ffree-line-length-none -g -fcheck=all -fbacktrace -Wall")

    target_compile_definitions(channel_core PRIVATE HAVE_FFTW)
    find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)
    target_link_libraries(channel_core PUBLIC FFTW::Double)


elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    find_package(CUDAToolkit REQUIRED)

    if(CMAKE_CUDA_ARCHITECTURES LESS 75)
        set(CMAKE_CUDA_ARCHITECTURES 75 80 86 90 100 120)
    endif()

    target_include_directories(channel_core PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_link_libraries(channel_core PUBLIC CUDA::cufft CUDA::cudart CUDA::cusparse)

    set(_NCCL_ROOT_HINTS)
    foreach(_NCCL_ENV_PATH IN ITEMS
        "$ENV{CPATH}"
        "$ENV{C_INCLUDE_PATH}"
        "$ENV{LD_LIBRARY_PATH}"
        "$ENV{LIBRARY_PATH}"
    )
        if(_NCCL_ENV_PATH)
            string(REPLACE ":" ";" _NCCL_ENV_HINTS "${_NCCL_ENV_PATH}")
            foreach(_NCCL_ENV_HINT IN LISTS _NCCL_ENV_HINTS)
                list(APPEND _NCCL_ROOT_HINTS "${_NCCL_ENV_HINT}")
                get_filename_component(_NCCL_ENV_ROOT "${_NCCL_ENV_HINT}" DIRECTORY)
                list(APPEND _NCCL_ROOT_HINTS "${_NCCL_ENV_ROOT}")
            endforeach()
        endif()
    endforeach()
    foreach(_NCCL_ROOT IN ITEMS
        "$ENV{NCCL_HOME}"
        "$ENV{NVHPC_ROOT}"
        "$ENV{NVCOMPILER}"
        "$ENV{NVHPC}"
    )
        if(_NCCL_ROOT)
            list(APPEND _NCCL_ROOT_HINTS "${_NCCL_ROOT}")
            file(GLOB _NCCL_VERSIONED_HINTS
                "${_NCCL_ROOT}/comm_libs/${CUDAToolkit_VERSION}/nccl*"
                "${_NCCL_ROOT}/comm_libs/${CUDAToolkit_VERSION}/openmpi*/latest/nccl"
                "${_NCCL_ROOT}/Linux_x86_64/*/comm_libs/${CUDAToolkit_VERSION}/nccl*"
                "${_NCCL_ROOT}/Linux_x86_64/*/comm_libs/${CUDAToolkit_VERSION}/openmpi*/latest/nccl"
            )
            if(_NCCL_VERSIONED_HINTS)
                list(SORT _NCCL_VERSIONED_HINTS)
                list(REVERSE _NCCL_VERSIONED_HINTS)
            endif()
            list(APPEND _NCCL_ROOT_HINTS ${_NCCL_VERSIONED_HINTS})
            file(GLOB _NCCL_GLOB_HINTS
                "${_NCCL_ROOT}/comm_libs/*/nccl*"
                "${_NCCL_ROOT}/comm_libs/*/openmpi*/latest/nccl"
                "${_NCCL_ROOT}/Linux_x86_64/*/comm_libs/*/nccl*"
                "${_NCCL_ROOT}/Linux_x86_64/*/comm_libs/*/openmpi*/latest/nccl"
                "${_NCCL_ROOT}/REDIST/comm_libs/*/nccl*"
                "${_NCCL_ROOT}/Linux_x86_64/*/REDIST/comm_libs/*/nccl*"
            )
            if(_NCCL_GLOB_HINTS)
                list(SORT _NCCL_GLOB_HINTS)
                list(REVERSE _NCCL_GLOB_HINTS)
            endif()
            list(APPEND _NCCL_ROOT_HINTS ${_NCCL_GLOB_HINTS})
        endif()
    endforeach()

    list(REMOVE_DUPLICATES _NCCL_ROOT_HINTS)
    set(_NCCL_CANDIDATE_ROOTS)
    foreach(_NCCL_CANDIDATE IN LISTS _NCCL_ROOT_HINTS)
        if(EXISTS "${_NCCL_CANDIDATE}/include/nccl.h"
           AND EXISTS "${_NCCL_CANDIDATE}/lib/libnccl.so")
            list(APPEND _NCCL_CANDIDATE_ROOTS "${_NCCL_CANDIDATE}")
        endif()
    endforeach()
    if(_NCCL_CANDIDATE_ROOTS)
        list(REMOVE_DUPLICATES _NCCL_CANDIDATE_ROOTS)
        message(STATUS "NCCL candidates found:")
        foreach(_NCCL_CANDIDATE IN LISTS _NCCL_CANDIDATE_ROOTS)
            read_nccl_header_version("${_NCCL_CANDIDATE}/include" _NCCL_CANDIDATE_VERSION)
            message(STATUS "  NCCL ${_NCCL_CANDIDATE_VERSION}: ${_NCCL_CANDIDATE}")
        endforeach()
    endif()

    foreach(_NCCL_CANDIDATE IN LISTS _NCCL_ROOT_HINTS)
        if(NOT NCCL_ROOT_DIR
           AND EXISTS "${_NCCL_CANDIDATE}/include/nccl.h"
           AND EXISTS "${_NCCL_CANDIDATE}/lib/libnccl.so")
            set(NCCL_ROOT_DIR "${_NCCL_CANDIDATE}" CACHE PATH "NCCL installation root" FORCE)
        endif()
    endforeach()

    if(NCCL_ROOT_DIR)
        set(NCCL_INCLUDE_DIR "${NCCL_ROOT_DIR}/include" CACHE PATH "NCCL include directory" FORCE)
        set(NCCL_LIBRARY "${NCCL_ROOT_DIR}/lib/libnccl.so" CACHE FILEPATH "NCCL library" FORCE)
    else()
        find_path(NCCL_INCLUDE_DIR nccl.h
            HINTS ${_NCCL_ROOT_HINTS}
            PATH_SUFFIXES include
        )
        find_library(NCCL_LIBRARY nccl
            HINTS ${_NCCL_ROOT_HINTS}
            PATH_SUFFIXES lib lib64
        )
    endif()
    if(NCCL_INCLUDE_DIR AND NCCL_LIBRARY)
        target_sources(channel_core PRIVATE channel_nccl_p2p.c)
        target_include_directories(channel_core PRIVATE ${NCCL_INCLUDE_DIR})
        target_link_libraries(channel_core PUBLIC ${NCCL_LIBRARY})
        target_compile_definitions(channel_core PRIVATE HAVE_NCCL)
        read_nccl_header_version("${NCCL_INCLUDE_DIR}" NCCL_VERSION)
        message(STATUS "NCCL y-pipeline backend enabled: version ${NCCL_VERSION}, ${NCCL_LIBRARY}")
    else()
        message(STATUS "NCCL not found; y-pipeline NCCL backend disabled")
    endif()
    if(TARGET CUDA::nvToolsExt)
        target_link_libraries(channel_core PUBLIC CUDA::nvToolsExt)
    elseif(CUDA_nvToolsExt_LIBRARY)
        target_link_libraries(channel_core PUBLIC ${CUDA_nvToolsExt_LIBRARY})
    endif()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O3 -mp -mp=gpu -Minfo=mp -Mcache_align -cuda") #
    target_compile_definitions(channel_core PRIVATE HAVE_CUDA)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -O0 -fp-model=precise")
    target_compile_definitions(channel_core PRIVATE HAVE_FFTW)
    find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)
    target_link_libraries(channel_core PUBLIC FFTW::Double)

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -e Z -O2 -hfp0 -homp -lrocprofiler-sdk-roctx")

    find_package(hipfort REQUIRED COMPONENTS hip hipfft hipsparse)
    target_link_libraries(channel_core PRIVATE hipfort::hip hipfort::hipfft hipfort::hipsparse)
    find_path(RCCL_INCLUDE_DIR
        NAMES rccl/rccl.h nccl.h
        HINTS ENV RCCL_HOME ENV ROCM_PATH /opt/rocm
        PATH_SUFFIXES include
    )
    find_library(RCCL_LIBRARY
        NAMES rccl
        HINTS ENV RCCL_HOME ENV ROCM_PATH /opt/rocm
        PATH_SUFFIXES lib lib64
    )
    if(RCCL_INCLUDE_DIR AND RCCL_LIBRARY)
        target_sources(channel_core PRIVATE channel_nccl_p2p.c)
        target_include_directories(channel_core PRIVATE ${RCCL_INCLUDE_DIR})
        target_link_libraries(channel_core PRIVATE ${RCCL_LIBRARY})
        target_compile_definitions(channel_core PRIVATE HAVE_NCCL HAVE_RCCL)
        read_nccl_header_version("${RCCL_INCLUDE_DIR}" RCCL_VERSION)
        message(STATUS "RCCL communication backend enabled: version ${RCCL_VERSION}, ${RCCL_LIBRARY}")
    else()
        message(STATUS "RCCL not found; CHANNEL_COMM=nccl backend disabled for HIP")
    endif()
    target_compile_definitions(channel_core PRIVATE HAVE_HIP)
    message(STATUS "hipfort_DIR = ${hipfort_DIR}")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O2 -fopenmp -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa --offload-arch=gfx942")

    find_package(hipfort REQUIRED COMPONENTS hip hipfft hipsparse)
    target_link_libraries(channel_core PUBLIC hipfort::hip hipfort::hipfft hipfort::hipsparse)
    find_path(RCCL_INCLUDE_DIR
        NAMES rccl/rccl.h nccl.h
        HINTS ENV RCCL_HOME ENV ROCM_PATH /opt/rocm
        PATH_SUFFIXES include
    )
    find_library(RCCL_LIBRARY
        NAMES rccl
        HINTS ENV RCCL_HOME ENV ROCM_PATH /opt/rocm
        PATH_SUFFIXES lib lib64
    )
    if(RCCL_INCLUDE_DIR AND RCCL_LIBRARY)
        target_sources(channel_core PRIVATE channel_nccl_p2p.c)
        target_include_directories(channel_core PRIVATE ${RCCL_INCLUDE_DIR})
        target_link_libraries(channel_core PUBLIC ${RCCL_LIBRARY})
        target_compile_definitions(channel_core PRIVATE HAVE_NCCL HAVE_RCCL)
        read_nccl_header_version("${RCCL_INCLUDE_DIR}" RCCL_VERSION)
        message(STATUS "RCCL communication backend enabled: version ${RCCL_VERSION}, ${RCCL_LIBRARY}")
    else()
        message(STATUS "RCCL not found; CHANNEL_COMM=nccl backend disabled for HIP")
    endif()
    target_compile_definitions(channel_core PRIVATE HAVE_HIP)
    message(STATUS "hipfort_DIR = ${hipfort_DIR}")

endif()
