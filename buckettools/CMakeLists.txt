# Require CMake 2.8
cmake_minimum_required(VERSION 2.8)

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)

# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# The project name for buckettools (sets up environment variables for binary and source directories)
project(buckettools)

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# Get DOLFIN configuration data (dolfin-config.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(dolfin)

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Add special DOLFIN compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DOLFIN_CXX_FLAGS}")

# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS} ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})

find_package(PythonLibs)

if(NOT ${PYTHONLIBS_FOUND})
  message(FATAL_ERROR "Could not find python libraries.")
endif(NOT ${PYTHONLIBS_FOUND})

include_directories(${PYTHON_INCLUDE_PATH})

find_package(PythonInterp)

if(NOT ${PYTHONINTERP_FOUND})
  message(FATAL_ERROR "Could not find python interpreter.")
endif(NOT ${PYTHONINTERP_FOUND})

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

find_package(Spud REQUIRED)

if(NOT ${SPUD_FOUND})
  message(FATAL_ERROR "Could not find Spud.")
endif(NOT ${SPUD_FOUND})

include_directories(${SPUD_INCLUDE_DIRS})

add_subdirectory(cpp)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

if (NOT DEFINED BUCKETTOOLS_PYTHON_INSTALL_DIR)
  # Get path for pure Python modules
  execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys, distutils.sysconfig; sys.stdout.write(distutils.sysconfig.get_python_lib(plat_specific=False, prefix='${CMAKE_INSTALL_PREFIX}'))"
    OUTPUT_VARIABLE BUCKETTOOLS_PYTHON_INSTALL_DIR
    )
  # Strip off CMAKE_INSTALL_PREFIX (is added later by CMake)
  string(REGEX REPLACE "${CMAKE_INSTALL_PREFIX}(/|\\\\)([^ ]*)" "\\2"
    BUCKETTOOLS_PYTHON_INSTALL_DIR "${BUCKETTOOLS_PYTHON_INSTALL_DIR}")
  set(BUCKETTOOLS_PYTHON_INSTALL_DIR ${BUCKETTOOLS_PYTHON_INSTALL_DIR}
    CACHE PATH "Python module installation directory.")
endif()

if (BUCKETTOOLS_PYTHON_INSTALL_DIR)
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/python/buckettools DESTINATION ${BUCKETTOOLS_PYTHON_INSTALL_DIR})
endif()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/schemas DESTINATION share/buckettools)

# Write a schemata file for bmls that can be installed
file(
    WRITE ${PROJECT_BINARY_DIR}/bml
    "BucketTools Markup Language\ndefault=${CMAKE_INSTALL_PREFIX}/share/buckettools/schemas/buckettools.rng"
    )

install(FILES ${PROJECT_BINARY_DIR}/bml DESTINATION share/buckettools/diamond/schemata)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/cmake DESTINATION share/buckettools)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/ufc   DESTINATION share/buckettools)

file(GLOB SCRIPTS ${PROJECT_SOURCE_DIR}/scripts/*)
install(FILES ${SCRIPTS} DESTINATION bin 
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)

