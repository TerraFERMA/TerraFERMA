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

add_subdirectory(poisson)
add_subdirectory(projection_mixedelement)
add_subdirectory(rbconvection_steadystate_split)
add_subdirectory(rbconvection_steadystate)
add_subdirectory(rbconvection1a)
add_subdirectory(rbconvection2a)

add_custom_target(run)
add_dependencies(run poissonrun projectionrun rbconvectionsplitrun rbconvectioncoupledrun rbconvection1arun rbconvection2arun)

add_custom_target(test)
add_dependencies(test projectiontest rbconvectionsplittest rbconvectioncoupledtest rbconvection1atest rbconvection2atest)

