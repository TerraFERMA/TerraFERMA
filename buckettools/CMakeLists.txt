# Copyright (C) 2013 Columbia University in the City of New York and others.
#
# Please see the AUTHORS file in the main source directory for a full list
# of contributors.
#
# This file is part of TerraFERMA.
#
# TerraFERMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TerraFERMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with TerraFERMA. If not, see <http://www.gnu.org/licenses/>.

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
find_package(DOLFIN)

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Add special DOLFIN compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DOLFIN_CXX_FLAGS}")

# Add Mac specific loader flags
if (APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-undefined,dynamic_lookup" CACHE STRING "should only be set if macosx" FORCE)
endif()

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

# Write a schemata file for tfmls that can be installed
file(
    WRITE ${PROJECT_BINARY_DIR}/tfml
    "TerraFERMA Markup Language\ndefault=${CMAKE_INSTALL_PREFIX}/share/buckettools/schemas/terraferma.rng"
    )

install(FILES ${PROJECT_BINARY_DIR}/tfml DESTINATION share/buckettools/diamond/schemata)

# Write a schemata file for tfmls that can be installed
file(
    WRITE ${PROJECT_BINARY_DIR}/shml
    "TerraFERMA Simulation Harness Markup Language\ndefault=${CMAKE_INSTALL_PREFIX}/share/buckettools/schemas/simulations.rng"
    )

install(FILES ${PROJECT_BINARY_DIR}/shml DESTINATION share/buckettools/diamond/schemata)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/cmake DESTINATION share/buckettools)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/ufc   DESTINATION share/buckettools)

file(GLOB SCRIPTS ${PROJECT_SOURCE_DIR}/scripts/*)
file(GLOB EXCLUDESCRIPTS ${PROJECT_SOURCE_DIR}/scripts/.*)
list(REMOVE_ITEM SCRIPTS ${EXCLUDESCRIPTS})
install(FILES ${SCRIPTS} DESTINATION bin 
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)

configure_file(${PROJECT_SOURCE_DIR}/cmake/templates/BUCKETTOOLSConfig.cmake.in
  ${PROJECT_BINARY_DIR}/BUCKETTOOLSConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/BUCKETTOOLSConfig.cmake DESTINATION share/buckettools/cmake)

