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
set(SOURCES
    rbmat.f90
    mpi_transpose.f90
    ffts.f90
    dnsdata.f90
    channel.f90
)

# --------------------------------------------------------------------
# Executable
# --------------------------------------------------------------------
add_executable(channel ${SOURCES})

target_link_libraries(channel PRIVATE MPI::MPI_Fortran)
target_link_libraries(channel PRIVATE FFTW::Double)

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