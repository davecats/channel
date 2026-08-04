# Syntax-checks the code guarded by the physics options.
#
# CHANNEL_HALF_CHANNEL and CHANNEL_PHI_NEUMANN select cpp-guarded branches in
# channel_bcs and case_setup.  Nothing in the default build compiles those
# branches, so a rename or a removal elsewhere could break them silently and
# stay broken for a long time -- which is exactly what happened to the
# bodyforce option before it was deleted.
#
# This compiles the two guarded files with both switches defined, syntax-only,
# reusing the .mod files the normal build already produced.  It is a second or
# two, and it catches the failure mode that matters: guarded code referring to
# something that no longer exists.  It does not check that the numbers are
# right -- nothing does, and that is a separate open question.
#
# Invoked with -DFC, -DFFLAGS, -DBUILD_DIR and -DSRC_DIR.

set(scratch "${BUILD_DIR}/physics_options_check")
file(REMOVE_RECURSE "${scratch}")
file(MAKE_DIRECTORY "${scratch}")

separate_arguments(flags NATIVE_COMMAND "${FFLAGS}")

set(sources
    "${SRC_DIR}/physics/channel_bcs.f90"
    "${BUILD_DIR}/generated/case_setup.f90")

foreach(source IN LISTS sources)
    if(NOT EXISTS "${source}")
        message(FATAL_ERROR "physics option check: missing source ${source}")
    endif()
    execute_process(
        COMMAND "${FC}" ${flags} -cpp -fsyntax-only
                -Dhalfchannel -DphiNeumann
                -I "${SRC_DIR}/util" -I "${BUILD_DIR}"
                -J "${scratch}"
                "${source}"
        RESULT_VARIABLE status
        OUTPUT_VARIABLE out
        ERROR_VARIABLE err)
    if(NOT status EQUAL 0)
        message(FATAL_ERROR
            "physics option check failed for ${source}\n"
            "with -Dhalfchannel -DphiNeumann\n${out}${err}")
    endif()
endforeach()

message(STATUS "physics options compile: OK")
