cmake_minimum_required(VERSION 3.16)
project(ChannelSolver LANGUAGES Fortran)

# --------------------------------------------------------------------
# CMake module path (for custom FindFFTW)
# --------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# --------------------------------------------------------------------
# MPI
# --------------------------------------------------------------------
find_package(MPI REQUIRED Fortran)

# --------------------------------------------------------------------
# FFTW
# --------------------------------------------------------------------
# Looks for FFTW with the FindFFTW.cmake module
find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB)

# --------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------
add_library(channel_core
    rbmat.f90
    mpi_transpose.f90
    ffts.f90
    dnsdata.f90
    driver.f90
)

target_link_libraries(channel_core PUBLIC MPI::MPI_Fortran FFTW::Double)

# --------------------------------------------------------------------
# Executable
# --------------------------------------------------------------------
add_executable(channel ${SOURCES} channel.f90)
target_link_libraries(channel PRIVATE channel_core)

# --------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------

set(TEST_NPROCS 1 CACHE STRING "Number of MPI ranks for regression test")

add_executable(test_regression test_regression.f90)
target_link_libraries(test_regression PRIVATE channel_core)

# Add a CTest target
enable_testing()

# Register the regression test
# You can pass input filenames here if needed

# 1-rank regression test
add_test(NAME regression_test_1rank
    COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} 1
            ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_regression>
            ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(regression_test_1rank PROPERTIES PROCESSORS 1)

# 4-rank regression test
add_test(NAME regression_test_4rank
    COMMAND ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test_regression>
            ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(regression_test_4rank PROPERTIES PROCESSORS 4)
# --------------------------------------------------------------------
# Compiler flags
# --------------------------------------------------------------------
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O0 -malign-double -fall-intrinsics -ffree-line-length-none")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O0 -Wall -Kieee -Mnofma -Mbounds -traceback")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -O0 -fp-model=precise")
endif()