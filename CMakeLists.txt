cmake_minimum_required(VERSION 3.23)

project(ChannelKokkosPort LANGUAGES CXX)

option(CHANNEL_USE_VENDORED_KOKKOS "Build the bundled third_party/kokkos checkout" ON)
option(CHANNEL_USE_VENDORED_UMPIRE "Build the bundled third_party/umpire checkout" ON)
option(CHANNEL_USE_VENDORED_KOKKOSFFT "Build the bundled third_party/kokkos-fft checkout" ON)
option(CHANNEL_BUILD_TESTS "Build C++ port tests" ON)
option(CHANNEL_ENABLE_CUDA "Enable the Kokkos CUDA backend for NVIDIA GPUs" OFF)
option(CHANNEL_ENABLE_HIP "Enable the Kokkos HIP backend for AMD GPUs" OFF)

if(CHANNEL_ENABLE_CUDA AND CHANNEL_ENABLE_HIP)
  message(FATAL_ERROR "CHANNEL_ENABLE_CUDA and CHANNEL_ENABLE_HIP are mutually exclusive")
endif()

function(channel_apply_source_patch patch_file work_dir description)
  if(NOT EXISTS "${work_dir}")
    message(FATAL_ERROR "${description} patch target is missing: ${work_dir}")
  endif()
  if(NOT EXISTS "${patch_file}")
    message(FATAL_ERROR "${description} patch file is missing: ${patch_file}")
  endif()

  execute_process(
    COMMAND git -C "${work_dir}" apply --check "${patch_file}"
    RESULT_VARIABLE CHANNEL_PATCH_CAN_APPLY
    OUTPUT_QUIET
    ERROR_QUIET
  )
  if(CHANNEL_PATCH_CAN_APPLY EQUAL 0)
    execute_process(
      COMMAND git -C "${work_dir}" apply "${patch_file}"
      RESULT_VARIABLE CHANNEL_PATCH_APPLY_RESULT
      OUTPUT_VARIABLE CHANNEL_PATCH_APPLY_OUTPUT
      ERROR_VARIABLE CHANNEL_PATCH_APPLY_ERROR
    )
    if(NOT CHANNEL_PATCH_APPLY_RESULT EQUAL 0)
      message(FATAL_ERROR "Failed to apply ${description} patch:\n${CHANNEL_PATCH_APPLY_OUTPUT}\n${CHANNEL_PATCH_APPLY_ERROR}")
    endif()
    message(STATUS "Applied ${description} patch")
    return()
  endif()

  execute_process(
    COMMAND git -C "${work_dir}" apply --reverse --check "${patch_file}"
    RESULT_VARIABLE CHANNEL_PATCH_ALREADY_APPLIED
    OUTPUT_QUIET
    ERROR_QUIET
  )
  if(CHANNEL_PATCH_ALREADY_APPLIED EQUAL 0)
    message(STATUS "${description} patch already applied")
  else()
    message(FATAL_ERROR "${description} patch cannot be applied cleanly and is not already applied")
  endif()
endfunction()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(channel_core
  src/compact_setup.cpp
  src/decomposition.cpp
  src/disabled_hooks.cpp
  src/dnsdata.cpp
  src/dns_solver.cpp
  src/ffts.cpp
  src/input.cpp
  src/mean_correction.cpp
  src/memory.cpp
  src/mpi_transpose.cpp
  src/mpi.cpp
  src/nonlinear.cpp
  src/restart.cpp
  src/runtime.cpp
  src/runge_kutta.cpp
  src/schur.cpp
  src/velocity_recovery.cpp
  src/workspace.cpp
  src/yline.cpp
)

target_include_directories(channel_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(channel_core PUBLIC cxx_std_20)

find_package(MPI REQUIRED COMPONENTS CXX)
target_compile_definitions(channel_core PUBLIC CHANNEL_HAS_MPI=1)
target_link_libraries(channel_core PUBLIC MPI::MPI_CXX)

if(CHANNEL_USE_VENDORED_KOKKOS)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/kokkos/CMakeLists.txt")
    set(Kokkos_ENABLE_TESTS OFF CACHE BOOL "Disable vendored Kokkos tests" FORCE)
    set(Kokkos_ENABLE_EXAMPLES OFF CACHE BOOL "Disable vendored Kokkos examples" FORCE)
    set(Kokkos_ENABLE_BENCHMARKS OFF CACHE BOOL "Disable vendored Kokkos benchmarks" FORCE)
    if(CHANNEL_ENABLE_CUDA)
      set(Kokkos_ENABLE_CUDA ON CACHE BOOL "Enable the Kokkos CUDA backend" FORCE)
      set(Kokkos_ENABLE_HIP OFF CACHE BOOL "Disable the Kokkos HIP backend" FORCE)
    endif()
    if(CHANNEL_ENABLE_HIP)
      set(Kokkos_ENABLE_HIP ON CACHE BOOL "Enable the Kokkos HIP backend" FORCE)
      set(Kokkos_ENABLE_CUDA OFF CACHE BOOL "Disable the Kokkos CUDA backend" FORCE)
    endif()
    if(NOT DEFINED Kokkos_ENABLE_SERIAL)
      set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "Enable the Kokkos Serial backend")
    endif()
    if(NOT DEFINED Kokkos_ENABLE_OPENMP)
      set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "Enable the Kokkos OpenMP backend")
    endif()
    add_subdirectory(third_party/kokkos EXCLUDE_FROM_ALL)
  else()
    message(FATAL_ERROR "CHANNEL_USE_VENDORED_KOKKOS=ON but third_party/kokkos is missing")
  endif()
else()
  find_package(Kokkos REQUIRED)
endif()
if(CHANNEL_ENABLE_CUDA AND NOT Kokkos_ENABLE_CUDA)
  message(FATAL_ERROR "CHANNEL_ENABLE_CUDA=ON requires a Kokkos build with Kokkos_ENABLE_CUDA=ON")
endif()
if(CHANNEL_ENABLE_HIP AND NOT Kokkos_ENABLE_HIP)
  message(FATAL_ERROR "CHANNEL_ENABLE_HIP=ON requires a Kokkos build with Kokkos_ENABLE_HIP=ON")
endif()
target_compile_definitions(channel_core PUBLIC CHANNEL_HAS_KOKKOS=1)
target_link_libraries(channel_core PUBLIC Kokkos::kokkos)

if(CHANNEL_USE_VENDORED_KOKKOSFFT)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/kokkos-fft/CMakeLists.txt")
    set(KokkosFFT_ENABLE_INTERNAL_KOKKOS OFF CACHE BOOL "Use channel's Kokkos target" FORCE)
    set(KokkosFFT_ENABLE_EXAMPLES OFF CACHE BOOL "Disable vendored KokkosFFT examples" FORCE)
    set(KokkosFFT_ENABLE_TESTS OFF CACHE BOOL "Disable vendored KokkosFFT tests" FORCE)
    set(KokkosFFT_ENABLE_BENCHMARK OFF CACHE BOOL "Disable vendored KokkosFFT benchmarks" FORCE)
    set(KokkosFFT_ENABLE_DOCS OFF CACHE BOOL "Disable vendored KokkosFFT docs" FORCE)
    set(KokkosFFT_ENABLE_DISTRIBUTED OFF CACHE BOOL "Disable vendored KokkosFFT distributed layer for now" FORCE)
    add_subdirectory(third_party/kokkos-fft EXCLUDE_FROM_ALL)
  else()
    message(FATAL_ERROR "CHANNEL_USE_VENDORED_KOKKOSFFT=ON but third_party/kokkos-fft is missing")
  endif()
else()
  find_package(KokkosFFT REQUIRED)
endif()
target_compile_definitions(channel_core PUBLIC CHANNEL_HAS_KOKKOSFFT=1)
target_link_libraries(channel_core PUBLIC KokkosFFT::fft)

set(CHANNEL_PATCH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cmake/patches")
if(CHANNEL_USE_VENDORED_UMPIRE)
  channel_apply_source_patch(
    "${CHANNEL_PATCH_DIR}/umpire-fmt-nvhpc-bit-cast.patch"
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/umpire/src/tpl/umpire/fmt"
    "vendored Umpire fmt NVHPC bit_cast"
  )
endif()
channel_apply_source_patch(
  "${CHANNEL_PATCH_DIR}/umpire-kokkos-missing-standard-includes.patch"
  "${CMAKE_CURRENT_SOURCE_DIR}/third_party/umpire-kokkos"
  "umpire-kokkos missing standard includes"
)

if(CHANNEL_USE_VENDORED_UMPIRE)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/umpire/CMakeLists.txt")
    set(ENABLE_TESTS OFF CACHE BOOL "Disable vendored Umpire tests" FORCE)
    set(ENABLE_BENCHMARKS OFF CACHE BOOL "Disable vendored Umpire benchmarks" FORCE)
    set(ENABLE_EXAMPLES OFF CACHE BOOL "Disable vendored Umpire examples" FORCE)
    set(ENABLE_FORTRAN OFF CACHE BOOL "Disable vendored Umpire Fortran" FORCE)
    set(ENABLE_OPENMP ${Kokkos_ENABLE_OPENMP} CACHE BOOL "Match Umpire OpenMP to Kokkos" FORCE)
    set(ENABLE_CUDA ${Kokkos_ENABLE_CUDA} CACHE BOOL "Match Umpire CUDA to Kokkos" FORCE)
    set(ENABLE_HIP ${Kokkos_ENABLE_HIP} CACHE BOOL "Match Umpire HIP to Kokkos" FORCE)
    set(ENABLE_MPI OFF CACHE BOOL "Umpire does not need MPI for channel workspace pooling" FORCE)
    set(UMPIRE_ENABLE_TESTS OFF CACHE BOOL "Disable vendored Umpire tests" FORCE)
    set(UMPIRE_ENABLE_BENCHMARKS OFF CACHE BOOL "Disable vendored Umpire benchmarks" FORCE)
    set(UMPIRE_ENABLE_EXAMPLES OFF CACHE BOOL "Disable vendored Umpire examples" FORCE)
    set(UMPIRE_ENABLE_FORTRAN OFF CACHE BOOL "Disable vendored Umpire Fortran" FORCE)
    add_subdirectory(third_party/umpire EXCLUDE_FROM_ALL)
  else()
    message(FATAL_ERROR "CHANNEL_USE_VENDORED_UMPIRE=ON but third_party/umpire is missing")
  endif()
else()
  find_package(umpire REQUIRED CONFIG)
endif()

find_path(UMPIRE_KOKKOS_INCLUDE_DIR UmpireSpace.hpp
          PATHS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/umpire-kokkos"
          NO_DEFAULT_PATH)
if(NOT UMPIRE_KOKKOS_INCLUDE_DIR)
  find_path(UMPIRE_KOKKOS_INCLUDE_DIR UmpireSpace.hpp REQUIRED)
endif()
target_compile_definitions(channel_core PUBLIC CHANNEL_HAS_UMPIRE_SPACE=1)
target_include_directories(channel_core PUBLIC ${UMPIRE_KOKKOS_INCLUDE_DIR})
if(TARGET umpire)
  target_link_libraries(channel_core PUBLIC umpire)
elseif(TARGET umpire::umpire)
  target_link_libraries(channel_core PUBLIC umpire::umpire)
else()
  message(FATAL_ERROR "Umpire target was not created")
endif()
message(STATUS "Using UmpireSpace at ${UMPIRE_KOKKOS_INCLUDE_DIR}")

set(CHANNEL_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/channel-bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CHANNEL_RUNTIME_OUTPUT_DIRECTORY}")

add_executable(channel src/main.cpp)
target_link_libraries(channel PRIVATE channel_core)

if(CHANNEL_BUILD_TESTS)
  enable_testing()
  set(CHANNEL_TEST_MPI_RANKS 1 CACHE STRING "MPI ranks used by focused C++ tests")

  function(add_channel_cpp_test name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE channel_core)
    add_test(NAME ${name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${CHANNEL_TEST_MPI_RANKS}
                     $<TARGET_FILE:${name}>)
  endfunction()

  add_channel_cpp_test(test_workspace_cpp tests/cpp/test_workspace.cpp)
  add_channel_cpp_test(test_input_cpp tests/cpp/test_input.cpp)
  add_channel_cpp_test(test_restart_io_cpp tests/cpp/test_restart_io.cpp)
  add_channel_cpp_test(test_dnsdata_cpp tests/cpp/test_dnsdata.cpp)
  add_channel_cpp_test(test_fft_cpp tests/cpp/test_fft.cpp)
  add_channel_cpp_test(test_mpi_transpose_cpp tests/cpp/test_mpi_transpose.cpp)
  add_test(NAME test_mpi_transpose_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_mpi_transpose_cpp>)
  add_channel_cpp_test(test_mean_correction_cpp tests/cpp/test_mean_correction.cpp)
  add_channel_cpp_test(test_nonlinear_product_transform_cpp tests/cpp/test_nonlinear_product_transform.cpp)
  add_test(NAME test_nonlinear_product_transform_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_nonlinear_product_transform_cpp>)
  add_channel_cpp_test(test_nonlinear_velocity_cpp tests/cpp/test_nonlinear_velocity.cpp)
  add_channel_cpp_test(test_velocity_recovery_cpp tests/cpp/test_velocity_recovery.cpp)
  add_channel_cpp_test(test_runge_kutta_cpp tests/cpp/test_runge_kutta.cpp)
  add_test(NAME test_mean_correction_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_mean_correction_cpp>)
  add_executable(test_dns_fft_transpose_cpp tests/cpp/test_dns_fft_transpose.cpp)
  target_link_libraries(test_dns_fft_transpose_cpp PRIVATE channel_core)
  add_test(NAME test_dns_fft_transpose_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${CHANNEL_TEST_MPI_RANKS}
                   $<TARGET_FILE:test_dns_fft_transpose_cpp>)
  add_test(NAME test_dns_fft_transpose_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_dns_fft_transpose_cpp>)
  add_executable(test_dns_linear_step_cpp tests/cpp/test_dns_linear_step.cpp)
  target_link_libraries(test_dns_linear_step_cpp PRIVATE channel_core)
  add_test(NAME test_dns_linear_step_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${CHANNEL_TEST_MPI_RANKS}
                   $<TARGET_FILE:test_dns_linear_step_cpp>)
  add_test(NAME test_dns_linear_step_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_dns_linear_step_cpp>)
  add_test(NAME test_channel_linear_fixture_2rank_alltoall
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:channel> --linear-fixture alltoall)
  add_test(NAME test_channel_linear_fixture_2rank_allgather
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:channel> --linear-fixture allgather)
  add_test(NAME test_channel_tiny_timestep_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                   $<TARGET_FILE:channel> --run ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/tiny_dns.in --steps 1)
  add_test(NAME test_channel_tiny_timestep_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:channel> --run ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/tiny_dns.in --steps 2)
  add_channel_cpp_test(test_pentadiagonal_cpp tests/cpp/test_pentadiagonal.cpp)
  add_channel_cpp_test(test_schur_level_cpp tests/cpp/test_schur_level.cpp)
  add_channel_cpp_test(test_schur_numeric_cpp tests/cpp/test_schur_numeric.cpp)
  add_executable(test_schur_solver_multirank_cpp tests/cpp/test_schur_solver_multirank.cpp)
  target_link_libraries(test_schur_solver_multirank_cpp PRIVATE channel_core)
  add_test(NAME test_schur_solver_multirank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_schur_solver_multirank_cpp>)
  add_executable(test_schur_solver_multilevel_cpp tests/cpp/test_schur_solver_multilevel.cpp)
  target_link_libraries(test_schur_solver_multilevel_cpp PRIVATE channel_core)
  add_test(NAME test_schur_solver_multilevel_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                   $<TARGET_FILE:test_schur_solver_multilevel_cpp>)
  add_test(NAME test_schur_solver_multilevel_6rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6
                   $<TARGET_FILE:test_schur_solver_multilevel_cpp>)
  add_test(NAME test_schur_solver_multilevel_6rank_root_alltoall_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6
                   $<TARGET_FILE:test_schur_solver_multilevel_cpp> root-alltoall)
  add_executable(test_endpoint_schur_yline_cpp tests/cpp/test_endpoint_schur_yline.cpp)
  target_link_libraries(test_endpoint_schur_yline_cpp PRIVATE channel_core)
  add_test(NAME test_endpoint_schur_yline_1rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                   $<TARGET_FILE:test_endpoint_schur_yline_cpp>)
  add_test(NAME test_endpoint_schur_yline_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_endpoint_schur_yline_cpp>)
  add_test(NAME test_endpoint_schur_yline_3rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3
                   $<TARGET_FILE:test_endpoint_schur_yline_cpp>)
  add_test(NAME test_endpoint_schur_yline_4rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                   $<TARGET_FILE:test_endpoint_schur_yline_cpp>)
  add_test(NAME test_endpoint_schur_yline_6rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6
                   $<TARGET_FILE:test_endpoint_schur_yline_cpp>)
  add_test(NAME test_endpoint_schur_yline_6rank_root_alltoall_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6
                   $<TARGET_FILE:test_endpoint_schur_yline_cpp> pass23)
  add_executable(test_dns_yline_component_cpp tests/cpp/test_dns_yline_component.cpp)
  target_link_libraries(test_dns_yline_component_cpp PRIVATE channel_core)
  add_test(NAME test_dns_yline_component_1rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                   $<TARGET_FILE:test_dns_yline_component_cpp>)
  add_test(NAME test_dns_yline_component_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_dns_yline_component_cpp>)
  add_test(NAME test_dns_yline_component_4rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                   $<TARGET_FILE:test_dns_yline_component_cpp>)
  add_test(NAME test_dns_yline_component_6rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6
                   $<TARGET_FILE:test_dns_yline_component_cpp>)
  add_test(NAME test_dns_yline_component_6rank_root_alltoall_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6
                   $<TARGET_FILE:test_dns_yline_component_cpp> pass23)
  add_executable(test_dns_implicit_yline_step_cpp tests/cpp/test_dns_implicit_yline_step.cpp)
  target_link_libraries(test_dns_implicit_yline_step_cpp PRIVATE channel_core)
  add_test(NAME test_dns_implicit_yline_step_1rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                   $<TARGET_FILE:test_dns_implicit_yline_step_cpp>)
  add_test(NAME test_dns_implicit_yline_step_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_dns_implicit_yline_step_cpp>)
  add_test(NAME test_dns_implicit_yline_step_4rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
                   $<TARGET_FILE:test_dns_implicit_yline_step_cpp>)
  add_test(NAME test_dns_implicit_yline_step_6rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6
                   $<TARGET_FILE:test_dns_implicit_yline_step_cpp>)
  add_test(NAME test_dns_implicit_yline_step_6rank_root_alltoall_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 6
                   $<TARGET_FILE:test_dns_implicit_yline_step_cpp> pass23)
  add_executable(test_regression_cpp tests/cpp/test_regression_cpp.cpp)
  find_library(CHANNEL_FFTW3_LIBRARY NAMES fftw3 REQUIRED)
  target_link_libraries(test_regression_cpp PRIVATE channel_core)
  target_link_libraries(test_regression_cpp PRIVATE ${CHANNEL_FFTW3_LIBRARY})
  add_test(NAME test_regression_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                   $<TARGET_FILE:test_regression_cpp>
                   --regression
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dns_test.in
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/start_field.out
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/end_field.out)
  add_test(NAME test_regression_2rank_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_regression_cpp>
                   --regression
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dns_test_npy2.in
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/start_field.out
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/end_field.out)
  add_test(NAME test_regression_2rank_npy1_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                   $<TARGET_FILE:test_regression_cpp>
                   --regression
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dns_test_npy1.in
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/start_field.out
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/end_field.out)
  add_executable(test_regression_scalar_cpp tests/cpp/test_regression_scalar_cpp.cpp)
  target_link_libraries(test_regression_scalar_cpp PRIVATE channel_core)
  add_test(NAME test_regression_scalar_cpp
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
                   $<TARGET_FILE:test_regression_scalar_cpp>
                   --regression-scalar
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/dns_test_scalar.in
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/start_field_scalar.out
                   ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/end_field_scalar.out)
  set_tests_properties(test_regression_scalar_cpp PROPERTIES SKIP_RETURN_CODE 77)
  add_channel_cpp_test(test_decomposition_cpp tests/cpp/test_decomposition.cpp)
  add_channel_cpp_test(test_mpi_autotune_candidates_cpp tests/cpp/test_mpi_autotune_candidates.cpp)
endif()
