# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# Get DOLFIN configuration data (dolfin-config.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(dolfin)

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Add special DOLFIN compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DOLFIN_CXX_FLAGS}")

# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS} ${DOLFIN_3RD_PARTY_INCLUDE_DIRS} ${PROJECT_BINARY_DIR})

file(
    WRITE ${PROJECT_BINARY_DIR}/hgid.h.in
    "#define __HG_ID__ \"@HG_ID@\"\n"
    )

file(
    WRITE ${PROJECT_BINARY_DIR}/hgid.cmake
    "find_program(HAVE_HG hg)
    if(HAVE_HG)
        execute_process (
           COMMAND hg id -i 
           WORKING_DIRECTORY \${DIR}
           OUTPUT_VARIABLE HG_ID
           ERROR_VARIABLE HG_ERR
           OUTPUT_STRIP_TRAILING_WHITESPACE
           ERROR_STRIP_TRAILING_WHITESPACE
           )
        if(HG_ERR)
           set(HG_ID "hg_id_not_found")
        endif(HG_ERR)
    else(HAVE_HG)
        set(HG_ID "hg_not_found")
    endif(HAVE_HG)
    configure_file(\${SRC} \${DST} @ONLY)
    "
    )

add_custom_target(
    hgid
    ${CMAKE_COMMAND} -D SRC="${PROJECT_BINARY_DIR}/hgid.h.in"
                     -D DST="${PROJECT_BINARY_DIR}/hgid.h"
                     -D DIR="${PROJECT_SOURCE_DIR}"
                     -P "${PROJECT_BINARY_DIR}/hgid.cmake"
    )

add_library(BucketTools StatFile.cpp DiagnosticsFile.cpp)
add_dependencies(BucketTools hgid)


