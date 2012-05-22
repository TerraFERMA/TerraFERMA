project(tests)

# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# Get DOLFIN configuration data (dolfin-config.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(dolfin)

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

if (NOT LOGLEVEL)
  set(LOGLEVEL "INFO" CACHE STRING
       "The log level that should be used in run targets, options are: DBG, TRACE, PROGRESS, INFO, WARNING, ERROR, CRITICAL." FORCE)
endif()

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Add special DOLFIN compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DOLFIN_CXX_FLAGS}")

# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS} ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})

set(BUCKETTOOLS_DIR "${PROJECT_SOURCE_DIR}/../buckettools" CACHE PATH
     "Location of buckettools." FORCE)

find_package(PythonLibs)

if(NOT ${PYTHONLIBS_FOUND})
  message(FATAL_ERROR "Could not find python libraries.")
endif(NOT ${PYTHONLIBS_FOUND})

include_directories(${PYTHON_INCLUDE_PATH})

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/../cmake/modules")

find_package(Spud REQUIRED)

if(NOT ${SPUD_FOUND})
  message(FATAL_ERROR "Could not find Spud.")
endif(NOT ${SPUD_FOUND})

include_directories(${SPUD_INCLUDE_DIRS})

add_subdirectory("${BUCKETTOOLS_DIR}/cpp" "buckettools_cpp")

include_directories("${BUCKETTOOLS_DIR}/include")

add_subdirectory(diamond_validation)
add_subdirectory(courant_1d)
add_subdirectory(alt_courant_1d)
add_subdirectory(adaptive_dt_1d)
add_subdirectory(alt_adaptive_dt_1d)
add_subdirectory(poisson)
add_subdirectory(neumann_poisson)
add_subdirectory(projection_mixedelement)
add_subdirectory(projection_timedependent)
add_subdirectory(projection_timedependent_checkpoint)
add_subdirectory(projection_cpp)
add_subdirectory(rbconvection_steadystate_split)
add_subdirectory(rbconvection_steadystate)
add_subdirectory(rbconvection1a)
add_subdirectory(rbconvection1a_refpt)
add_subdirectory(rbconvection2a)
add_subdirectory(region_ids)
add_subdirectory(subductionbenchmark1a)
add_subdirectory(subductionbenchmark1b)
add_subdirectory(subductionbenchmark1c)
add_subdirectory(subductionbenchmark1c_split)
add_subdirectory(subductionbenchmark2a)
add_subdirectory(subductionbenchmark2a_split)
add_subdirectory(subductionbenchmark2b)
add_subdirectory(dg_advection_1d_vi)

add_custom_target(run)
add_dependencies(run poisson_run)

add_custom_target(test)
add_dependencies(test diamond_validation_test courant_1d_test alt_courant_1d adaptive_dt_1d_test alt_adaptive_dt_1d_test
projection_mixedelement_test projection_timedependent_test projection_timedependent_checkpoint_test projection_cpp_test
rbconvection_steadystate_split_test rbconvection_steadystate_test rbconvection1a_test rbconvection1a_refpt_test region_ids_test 
subductionbenchmark1a_test subductionbenchmark1b_test subductionbenchmark1c_test subductionbenchmark1c_split_test
subductionbenchmark2a_test subductionbenchmark2a_split_test subductionbenchmark2b_test neumann_poisson_test)

add_custom_target(test3p2)
add_dependencies(test3p2 test dg_advection_1d_vi_test)

add_custom_target(mediumtest)
add_dependencies(mediumtest rbconvection2a_test)

