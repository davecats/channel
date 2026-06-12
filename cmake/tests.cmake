function(add_channel_mpi_test)
    set(options USE_YSLAB_MATRIX)
    set(one_value_args NAME TARGET NPROCS PROCESSORS)
    set(multi_value_args ARGS ENVIRONMENT)
    cmake_parse_arguments(TEST "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    if(NOT TEST_NAME OR NOT TEST_TARGET OR NOT TEST_NPROCS OR NOT TEST_PROCESSORS)
        message(FATAL_ERROR "add_channel_mpi_test requires NAME, TARGET, NPROCS, and PROCESSORS")
    endif()

    set(test_command
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} ${TEST_NPROCS}
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:${TEST_TARGET}>
        ${TEST_ARGS}
        ${MPIEXEC_POSTFLAGS}
    )

    if(TEST_USE_YSLAB_MATRIX)
        foreach(env_value IN ITEMS 0 1)
            set(test_env ${TEST_ENVIRONMENT})
            if(env_value EQUAL 0)
                set(test_suffix schur)
            else()
                set(test_suffix yslab)
            endif()
            list(APPEND test_env "CHANNEL_USE_YSLAB_LINSOLVE=${env_value}")
            set(test_name "${TEST_NAME}_${test_suffix}")
            add_test(NAME ${test_name}
                COMMAND ${test_command}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            )
            set_tests_properties(${test_name} PROPERTIES
                PROCESSORS ${TEST_PROCESSORS}
                ENVIRONMENT "${test_env}"
            )
        endforeach()
        return()
    endif()

    add_test(NAME ${TEST_NAME}
        COMMAND ${test_command}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
    set_tests_properties(${TEST_NAME} PROPERTIES
        PROCESSORS ${TEST_PROCESSORS}
        ENVIRONMENT "${TEST_ENVIRONMENT}"
    )
endfunction()

add_channel_mpi_test(
    NAME regression_test_1rank
    TARGET test_regression
    NPROCS 1
    PROCESSORS 1
)

add_channel_mpi_test(
    NAME regression_test_4rank
    TARGET test_regression
    NPROCS 2
    PROCESSORS 4
)

add_channel_mpi_test(
    NAME regression_test_npy2_4rank
    TARGET test_regression
    NPROCS 4
    PROCESSORS 4
    USE_YSLAB_MATRIX
    ARGS tests/data/dns_test_npy2.in
)

add_channel_mpi_test(
    NAME pressure_dpdy_1rank
    TARGET test_pressure_dpdy
    NPROCS 1
    PROCESSORS 1
)

add_channel_mpi_test(
    NAME pressure_dpdy_2rank
    TARGET test_pressure_dpdy
    NPROCS 2
    PROCESSORS 2
)

add_channel_mpi_test(
    NAME pressure_dpdy_npy2_2rank
    TARGET test_pressure_dpdy
    NPROCS 2
    PROCESSORS 2
    USE_YSLAB_MATRIX
    ARGS tests/data/dns_test_npy2.in
)

add_channel_mpi_test(
    NAME pressure_dpdy_npy2_4rank
    TARGET test_pressure_dpdy
    NPROCS 4
    PROCESSORS 4
    USE_YSLAB_MATRIX
    ARGS tests/data/dns_test_npy2.in
)

add_channel_mpi_test(
    NAME convvelo_stats_1rank
    TARGET test_convvelo_stats
    NPROCS 1
    PROCESSORS 1
)

add_channel_mpi_test(
    NAME convvelo_runtime_1rank
    TARGET test_convvelo_runtime
    NPROCS 1
    PROCESSORS 1
    ARGS full
)

add_channel_mpi_test(
    NAME convvelo_runtime_minimal_1rank
    TARGET test_convvelo_runtime
    NPROCS 1
    PROCESSORS 1
    ARGS minimal
)

add_channel_mpi_test(
    NAME regression_test_scalar_1rank
    TARGET test_regression_scalar
    NPROCS 1
    PROCESSORS 1
    ARGS tests/data/dns_test_scalar.in
)

add_channel_mpi_test(
    NAME regression_test_scalar_4rank
    TARGET test_regression_scalar
    NPROCS 2
    PROCESSORS 4
    ARGS tests/data/dns_test_scalar.in
)

add_channel_mpi_test(
    NAME regression_test_scalar_6rank_np4
    TARGET test_regression_scalar
    NPROCS 6
    PROCESSORS 8
    USE_YSLAB_MATRIX
    ARGS tests/data/dns_test_scalar_npy3.in
)

add_channel_mpi_test(
    NAME regression_test_scalar_npy4_1rank
    TARGET test_regression_scalar
    NPROCS 1
    PROCESSORS 1
    ARGS tests/data/dns_test_scalar_npy4.in tests/data/start_field_scalar_npy4.out tests/data/end_field_scalar_npy4.out
    ENVIRONMENT CHANNEL_NPY=1
)

add_channel_mpi_test(
    NAME regression_test_scalar_npy4_4rank
    TARGET test_regression_scalar
    NPROCS 4
    PROCESSORS 4
    ARGS tests/data/dns_test_scalar_npy4.in tests/data/start_field_scalar_npy4.out tests/data/end_field_scalar_npy4.out
    ENVIRONMENT CHANNEL_USE_YSLAB_LINSOLVE=0
)

add_channel_mpi_test(
    NAME regression_test_scalar_npy4_8rank
    TARGET test_regression_scalar
    NPROCS 8
    PROCESSORS 8
    ARGS tests/data/dns_test_scalar_npy4.in tests/data/start_field_scalar_npy4.out tests/data/end_field_scalar_npy4.out
    ENVIRONMENT CHANNEL_USE_YSLAB_LINSOLVE=0
)

set(POST_PRESSURE_NPROCS 1)
set(POST_PRESSURE_EXE ${CMAKE_CURRENT_BINARY_DIR}/post_pressure)
configure_file(cmake/run_post_pressure_test.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/run_post_pressure_test.cmake @ONLY)

add_test(NAME post_pressure_1rank
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/run_post_pressure_test.cmake
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(post_pressure_1rank PROPERTIES PROCESSORS 1)

set(POST_CONVVELO_NPROCS 1)
set(POST_CONVVELO_EXE ${CMAKE_CURRENT_BINARY_DIR}/post_convvelo)
set(TEST_POST_CONVVELO_EXE ${CMAKE_CURRENT_BINARY_DIR}/test_post_convvelo)
configure_file(cmake/run_post_convvelo_test.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/run_post_convvelo_test.cmake @ONLY)

add_test(NAME post_convvelo_1rank
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/run_post_convvelo_test.cmake
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(post_convvelo_1rank PROPERTIES PROCESSORS 1)

add_channel_mpi_test(
    NAME reduced_ghost_backend_2rank
    TARGET test_y_pencil_transpose
    NPROCS 2
    PROCESSORS 2
    USE_YSLAB_MATRIX
)

add_channel_mpi_test(
    NAME mean_correction_npy2_2rank
    TARGET test_mean_correction
    NPROCS 2
    PROCESSORS 2
    USE_YSLAB_MATRIX
)
