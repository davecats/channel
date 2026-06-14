cmake_minimum_required(VERSION 3.16)
project(ChannelSolver LANGUAGES Fortran C)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

# --------------------------------------------------------------------
# CMake module path (for custom FindFFTW)
# --------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


# --------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------
add_library(channel_core
    config.f90
    rbmat.f90
    roctx.f90
    y_schur_solver.f90
    y_line_solvers.f90
    compact_line_solvers.f90
    mpi_transpose.f90
    ffts.f90
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
    target_link_libraries(channel_core PRIVATE hipfort::hipfft hipfort::hipsparse)
    target_compile_definitions(channel_core PRIVATE HAVE_HIP)
    message(STATUS "hipfort_DIR = ${hipfort_DIR}")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O2 -fopenmp -fopenmp -fopenmp-targets=amdgcn-amd-amdhsa --offload-arch=gfx942")

    find_package(hipfort REQUIRED COMPONENTS hip hipfft hipsparse)
    target_link_libraries(channel_core PUBLIC hipfort::hipfft hipfort::hipsparse)
    target_compile_definitions(channel_core PRIVATE HAVE_HIP)
    message(STATUS "hipfort_DIR = ${hipfort_DIR}")

endif()
